# Script for Running Final Model.
# list.of.packages <- c("tidyverse", "tidymodels")

options(scipen = 999)
library(tidymodels)
source("data_prep.R")
source("hic_code.R")

# function for running the entire set of models
# Input: response of interest, which structure to use (wweia for non-hierarchical), and food_value measure
# Output: model training and testing results are written to model registry for all feature selection algorithms
#       and all values of k.
run_all <- function(response, structure, food_value,model = "nn"){

  # get the necessary data written to global environment
  get_data(response, structure,food_value)
  
  print(str_c("Got ", response, " data"))
  
  mod_formula <<- formula(str_c(response," ~ ", " . "))
  # Declare covariates to be used.
    if(response=="waist_circ_risk"){
      covariates <- c("ridageyr", "riagendr", "education", "race_ethnicity", "smoking","indfmpir", "physical_activity") 
    } else{
      covariates <- c("ridageyr", "riagendr", "education",'bmxbmi', "race_ethnicity", "smoking","indfmpir", "physical_activity") 
    }
  # vector of feature selection algorithms to consider for the structure
  if(structure=="fndds"){
    fs_vec <- c("hic_mod_coef_raw")
  } else {
    fs_vec <- c("hic_mod_coef_raw","mi","cv","su","gr","test","cor", "relieff")
  }
  # create dataframe of feature selection/number of features, k, to pass into prep_run_models function
  fs_grid <- expand.grid(fs = fs_vec,top_k = c(5,10,50,100,200,500), stringsAsFactors = FALSE)
  results_list <- pmap(fs_grid, .f = ~prep_run_models(fs = ..1, top_k = ..2,
                                                      model=model,response=response,
                                                      structure=structure, food_value=food_value,
                                                      covariates=covariates))
  
}

# Function that runs the model and outputs test and training set features.
# Input: feature selection algorithm (fs), number of food features (top_k), model type (model), 
#       response, structure, food_value
# Output: test and train model results for combination of inputs.
prep_run_models <- function(fs = "hic_mod_coef_raw", top_k =10, model = "lr", response, structure, food_value,covariates){

  fs_summary <- get_fs_summary(response, structure,food_value, relieff=TRUE,resampling)
  fs_column <- str_glue(structure,response,food_value,fs,"rank",.sep ="_")
  predictors <- fs_summary %>% select(fdcd,drxfcsd,as.name(fs_column)) %>% filter(!!sym(fs_column)<=top_k) %>% pull(fdcd)
  
  # Get train data with all predictors
  train_data_all <-
    train_data %>%
    left_join(seqn_fdcd_df %>% select(seqn,all_of(predictors)), by = "seqn") %>%
    select(seqn,sub_weights,all_of(response), all_of(predictors), all_of(covariates)) 
 
  # Get test data with all predictors
  test_data_all <- 
    test_data %>%
    left_join(seqn_fdcd_df %>% select(seqn,all_of(predictors)), by = "seqn")  %>%
    select(seqn,sub_weights,all_of(response), all_of(predictors), all_of(covariates))
      
  # conduct 6-fold cross validation
  set.seed(123)
  folds <- vfold_cv(train_data_all, v=6)
  
  # run the model with train and test data.
  model_out <- run_model(train_data_all, test_data_all, folds, response,structure,fs,top_k,model,food_value)
  
  # register model results in registry
  register_model(model_out %>% 
                   mutate(food_value= food_value))
  return(model_out )

}

# function that runs whichever model is requested on train and test data.
run_model <- function(train_data_all, test_data_all, folds, response,structure,fs,top_k,model,food_value){
  
  # establish recipe and preprocessing of data
  mod_recipe<-
    train_data_all %>% 
    recipe(mod_formula, data=.) %>% 
    step_impute_mode(all_nominal_predictors()) %>%
    step_impute_median(all_numeric_predictors()) %>%
    step_naomit(all_predictors(), skip = F) %>%
    step_dummy(all_nominal_predictors()) %>% 
    {if(model == "nn") step_normalize(.,all_predictors()) else .}

  # establish model engine and parameter tuning
  if(model == "rf"){
    tune_spec <-
      rand_forest(
        mtry = tune(),
        trees =  500,
        min_n =  tune()) %>%
      set_engine("ranger", importance = "impurity") %>%
      set_mode("classification") 
  } else if(model == "lr"){
    tune_spec <-
      logistic_reg(
        penalty=tune()
      ) %>%
      set_engine("glmnet", importance = "impurity") %>%
      set_mode("classification")
  } else if(model == "lasso"){
    tune_spec <-
      logistic_reg(
        penalty=tune(),
        mixture = 1
      ) %>%
      set_engine("glmnet", importance = "impurity") %>%
      set_mode("classification")
  } else if(model == "nn"){
    tune_spec <-
      mlp(
        penalty = tune(),
        hidden_units=tune(),
        epochs = 500
      ) %>%
      set_engine("nnet", MaxNWts=31000) %>% 
      set_mode("classification")
  }

  # Initialize model workflow
  tune_wf <-
    workflow() %>%
    add_model(tune_spec) %>%
    add_recipe(mod_recipe)
  
  # tune the model on crossvalidation folds
  if(model == "rf"){
    tune_output <- tune_grid(
     tune_wf, # workflow
     resamples = folds, # cv folds
     metrics = metric_set(roc_auc, pr_auc,accuracy, precision, recall, f_meas,average_precision), # metrics to track
     control = control_grid(save_pred = TRUE, event_level = 'second',pkgs = c("tidyverse")),
     grid = 10)
  } else if(model %in% c("lr","lasso")){
    lr_reg_grid <- tibble(penalty = 10^seq(-4, -1, length.out = 10))
    
    tune_output <- tune_grid(
      tune_wf, # workflow
      resamples = folds, # cv folds
      metrics = metric_set(roc_auc, pr_auc,accuracy, precision, recall, f_meas, average_precision),
      control = control_grid(save_pred = TRUE, event_level = 'second', pkgs = c("tidyverse")),
      grid = lr_reg_grid # penalty grid defined above
    )
  } else if(model== "nn"){
    nn_grid <- expand.grid(penalty = c(0,.1,.01),hidden_units = seq(1,3, by=1))
    
    tune_output <- tune_grid(
      tune_wf, # workflow
      resamples = folds, # cv folds
      metrics = metric_set(roc_auc, pr_auc,accuracy, precision, recall, f_meas, average_precision),
      control = control_grid(save_pred = TRUE, event_level = 'second', pkgs = c("tidyverse")),
      grid = 10 # penalty grid defined above
    )
  }
  
  # select model with best auroc
  best_auc <- select_best(tune_output, "roc_auc")

  results_df <- bind_rows(tune_output$.metrics) %>%
    {if(model=="lr") filter(., penalty==best_auc$penalty) else .} %>% 
    {if(model=="rf") filter(.,mtry==best_auc$mtry & min_n == best_auc$min_n) else .} %>% 
    {if(model=="nn") filter(.,penalty==best_auc$penalty & hidden_units == best_auc$hidden_units) else .} 
  # gather training metrics
  train_results <- get_model_metrics(results_df, "train",structure,response,fs,top_k,model,resampling, food_value, which_cov) 
 
  # create final model based on best auroc
  final_mod <- finalize_model(tune_spec, best_auc)
      
  final_wf <-
    workflow() %>%
    add_model(final_mod) %>%
    add_recipe(mod_recipe)

  # fit final model on all training data
  set.seed(456)
  final_rf_fit <- 
    final_wf %>% 
    fit(data = train_data_all)

  rf_testing_pred <-
    predict(final_rf_fit, test_data_all) %>% 
    bind_cols(predict(final_rf_fit, test_data_all, type = "prob")) %>% 
    bind_cols(test_data_all %>% select(as.name(response)))
  
  rf_testing_pred %>% count(.pred_class)
  
  test_results <- get_model_metrics(rf_testing_pred, "test",structure,response,fs,top_k,model,resampling, food_value,which_cov)
  
  return(bind_rows(test_results, train_results))
  register_model(bind_rows(test_results, train_results)%>% 
                   mutate(food_value= food_value,
                          covariates_used= which_cov))
  # write_csv(bind_rows(test_results, train_results), "ults/model_registry.csv")
  # roc = pROC::roc(response = rf_testing_pred %>% pull(!!sym(response)),
  #                  predictor = rf_testing_pred$.pred_1)
  # pROC::ci.auc(roc, method = "bootstrap")
  # pROC::ci.auc(roc)
  yardstick::roc_auc(test_data,ada_class,truth = 'diabetic',event_level = "second")
  yardstick::average_precision(test_data,ada_class,truth = 'diabetic',event_level = "second")
  # pROC::auc(roc)
#### Logistic Regression / Random Forest / NN End
}

#### Lasso Begin
run_lasso <- function(response){
  # Logistic LASSO Regression Model Spec
  logistic_lasso_spec_tune <- logistic_reg() %>%
    set_engine('glmnet') %>%
    set_args(mixture = 1, penalty = tune()) %>%
    set_mode('classification')
  
  # Recipe
  logistic_rec <- 
    recipe(mod_formula, data = train_data_all) %>%
    update_role(seqn, new_role = "ID") %>% 
    # step_naomit(all_predictors(), skip = T) %>% 
    step_dummy(all_nominal_predictors())
      # step_normalize(all_numeric_predictors()) %>% 
    # step_dummy(all_nominal_predictors())
  
  # Workflow (Recipe + Model)
  log_lasso_wf <- workflow() %>% 
    add_recipe(logistic_rec) %>%
    add_model(logistic_lasso_spec_tune)     
  
  # Tune Model (trying a variety of values of Lambda penalty)
  penalty_grid <- grid_regular(
    penalty(range = c(-3, 0)), #log10 transformed
    levels = 10)
  
  tune_output <- tune_grid( 
    log_lasso_wf, # workflow
    resamples = folds, # cv folds
    metrics = metric_set(roc_auc, accuracy, precision, recall, f_meas),
    control = control_resamples(save_pred = TRUE, event_level = 'second'),
    grid = penalty_grid # penalty grid defined above
  )
  
  tune_output$.metrics
  
  # Visualize Model Evaluation Metrics from Tuning
  # autoplot(tune_output) + theme_classic()
  
  # Select Penalty
  if(response=="diabetic"){
    best_se_penalty <- select_best(tune_output, metric = 'average_precision', desc(penalty)) # choose penalty value based on the largest penalty within 1 se of the lowest CV roc_auc
  } else {
    best_se_penalty <- select_best(tune_output, metric = 'roc_auc', desc(penalty)) # choose penalty value based on the largest penalty within 1 se of the lowest CV roc_auc
  }
  # best_se_penalty <- show_best(tune_output, metric = 'roc_auc') # choose penalty value based on the largest penalty within 1 se of the lowest CV roc_auc
  best_se_penalty
  
  fold_average <- bind_rows(tune_output$.metrics) %>% 
    filter(penalty == best_se_penalty$penalty) %>% 
    group_by(.metric) %>% 
    summarize(metric_average = mean(.estimate))
  
  # Fit Final Model
  final_fit_se <- finalize_workflow(log_lasso_wf, best_se_penalty) %>% # incorporates penalty value to workflow 
    fit(data = train_data)
  
  lasso_output <- final_fit_se %>% tidy() %>% left_join(food_codes_all %>% select(category_number,drxfcsd) %>% mutate(category_number=as.character(category_number)),by= c("term"="category_number"))
  lasso_output %>% write_rds("02_results/lasso_wweia_bp_coeff.rds")
  
  lasso_aug <- 
    augment(final_fit_se, test_data) %>% 
    select(seqn,as.name(response),.pred_class,.pred_0,.pred_1)
  lasso_aug %>% write_rds("02_results/lasso_wweia_bp_results.rds")
  
  # The data look like: 
  lasso_aug %>%
    select(as.name(response), .pred_class, .pred_1)
  
  lasso_aug %>% 
    roc_curve(truth = !!sym(as.name(response)), .pred_0) %>% 
    autoplot()
  
  test_metrics <- bind_rows(accuracy(lasso_aug, truth = !!sym(as.name(response)), .pred_class),
                            pr_auc(lasso_aug, truth = !!sym(as.name(response)), .pred_class, event_level="second"),
                            precision(lasso_aug, truth = !!sym(as.name(response)), .pred_class, event_level="second"),
                            recall(lasso_aug, truth = !!sym(as.name(response)), .pred_class, event_level="second"),
                            f_meas(lasso_aug, truth = !!sym(as.name(response)), .pred_class, event_level="second"),
                            average_precision(lasso_aug, truth = !!sym(as.name(response)), .pred_class, event_level="second"),
                            roc_auc(lasso_aug,truth = !!sym(as.name(response)), .pred_1, event_level="second"))
}  
####### End Lasso 

get_model_metrics <- function(predictions = rf_testing_pred, data_split="test",structure,response,fs,top_k,model,resampling=resampling,food_value= food_value,which_cov){
  model_meta <- tibble_row(date_time= lubridate::now(),
                           structure = structure, 
                           response = response, 
                           resampling = resampling, 
                           feature_selection = ifelse(fs=="hic_test",str_c("hic_",structure),fs), 
                           top_k = top_k, 
                           model = model,
                           data= data_split,
                           food_value = food_value,
                           covariates_used = which_cov)
  if(data_split=="train" & model !='glm'){
    return(bind_cols(model_meta,
                     predictions %>% 
                       group_by(.metric) %>% 
                       summarize(mean_fold_metric = mean(.estimate, na.rm=TRUE)) %>% 
                       pivot_wider(names_from = .metric, values_from = mean_fold_metric) %>% 
                       mutate(data_prev = read_csv("01_data/prevalences.csv") %>% filter(response==!!response) %>%  pull(train_prev))) %>% 
             relocate(pr_auc,.after=roc_auc))
  } else {
    conf_mat <- conf_mat(predictions, truth = !!sym(as.name(response)), estimate= .pred_class)
    test_metrics <- bind_rows(accuracy(predictions, truth = !!sym(as.name(response)),estimate=  .pred_class),
                              precision(predictions, truth = !!sym(as.name(response)),estimate=  .pred_class,event_level='second'),
                              recall(predictions, truth = !!sym(as.name(response)),estimate=  .pred_class,event_level='second'),
                              f_meas(predictions, truth = !!sym(as.name(response)),estimate=  .pred_class,event_level='second'),
                              roc_auc(predictions,truth = !!sym(as.name(response)),.pred_1,event_level='second'),
                              pr_auc(predictions, truth = !!sym(as.name(response)), .pred_1, event_level='second'),
                              average_precision(predictions, truth = !!sym(as.name(response)), .pred_1, event_level='second')) %>% 
      select(-.estimator) %>% 
      pivot_wider(names_from = .metric, values_from = .estimate) %>% 
      mutate(true_pos = conf_mat$table[2,2], false_pos = conf_mat$table[2,1], true_neg = conf_mat$table[1,1], false_neg = conf_mat$table[1,2],
             pred_prev = (true_pos + false_pos)/(true_pos + false_neg + false_pos + true_neg),
             data_prev = read_csv("01_data/prevalences.csv") %>% filter(response==!!response) %>%  pull(test_prev))
    return(bind_cols(model_meta, test_metrics))
  }
}

### Register Model
register_model <- function(results){
  # old_registry <- read_csv("02_results/model_registry.csv")
  old_registry <- read_csv("02_results/model_registry.csv",guess_max = 100000)# %>% 
    # mutate(structure = ifelse(feature_selection == "hic_fndds", "fndds",structure)) %>% 
    # write_csv("02_results/model_registry.csv")
  write_csv(bind_rows(old_registry,results), "02_results/model_registry.csv")
}
