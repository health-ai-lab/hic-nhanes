options(java.parameters = "-Xmx6g")
update_data=FALSE
source("data_prep.R")

# response options: triglycerides_risk, reduced_HDL_risk, elevated_blood_pressure_risk,
#                   elevated_fasting_glucose_risk, waist_circ_risk, at_risk
# structure options: fndds (for HIC FNDDS only), wweia (HIC WWEIA and all non-hierarchical methods)
# food_value options: food_consumed, cal_consumed, prop_cal

library(foreach)  
library(FSelector)
run_hic <- function(response = "diabetic", structure = "fndds",food_value = "food_consumed"){
  
  get_data(response, structure,food_value)
  
  # Covariates differ for waist_circ_risk than the other outcomes    
  if(response=="waist_circ_risk"){
    covariates <- c("ridageyr", "riagendr", "education", "race_ethnicity", "smoking","indfmpir", "physical_activity") 
  } else{
    covariates <- c("ridageyr", "riagendr", "education",'bmxbmi', "race_ethnicity", "smoking","indfmpir", "physical_activity") 
  }
  
  # Join food consumption data with the training data for the participants
  feature_sel_df <-
    train_data %>% 
    select(seqn, as.name(response), all_of(covariates)) %>% 
    inner_join(seqn_fdcd_df)
  

    if(food_value != "food_consumed"){
      food_consumed_results <<- read_csv(str_c("01_data/hic_results/hic_",structure,"_",response,"_food_consumed.csv"), col_types=cols()) %>% 
        mutate(fdcd=trimws(as.character(format(fdcd, scientific=FALSE))))
    } 
  
  # create global dictionary of necessary proportions for each food/response combination
  prop_dictionary <<- make_prop_dictionary(feature_sel_df %>% 
                                             select(-all_of(covariates)) %>% 
                                             mutate(!!response := as.numeric(as.character(feature_sel_df %>% pull(as.name(response))))),
                                           response, food_value,structure)
 
  fdcds <- 
    prop_dictionary %>% 
    filter(n_fdcd>0  & fdcd!=response) %>% 
    filter(!is.na(n_fdcd)) %>% 
    arrange(fdcd) %>% 
    pull(fdcd)
  print("get_hic_terms")

  # Get all terms for each food by response and food_value 
  terms_sub <<- foreach(fdcd = fdcds, .combine = bind_rows, .packages = "tidyverse") %do% get_hic_terms_raw(fdcd, feature_sel_df, food_value, response,covariates)
  terms <- bind_rows(terms_sub)
 
  # calculate alpha value used to determine the branch and tree statistical significance weights.
  alpha <-
    prop_dictionary %>% 
    select(fdcd, n_fdcd, p_x, n) %>% 
    filter(fdcd != response) %>% 
    left_join(food_codes_all %>% 
                select(category_number, max_level_skips) %>% 
                mutate(category_number=trimws(as.character(format(category_number,scientific = FALSE)))), 
              by = c("fdcd"='category_number')) %>% 
    mutate(numerator = ifelse(max_level_skips==1,p_x*logb(n_fdcd,max_level_skips+1),p_x*logb(n_fdcd,max_level_skips))) %>% 
    replace_na(list(numerator = 0)) %>% #view(.)
    summarize(alpha_num = sum(numerator), 
              alpha_den=(logb(max(n),max(max_level_skips))*logb(length(.$fdcd),max(max_level_skips))),
              alpha = alpha_num/alpha_den) %>% 
    pull(alpha)
  
  # Calculate weights for branch and tree statistical significance
  terms_weights <- terms %>% 
    mutate(w_b = 1/(alpha + 1),
              w_t = alpha/(alpha+1))
  
  # calculate the final term values for each of the 5 hic terms and the final hic score along with all other ranks.
  hic_df_temp <- terms_weights %>% 
    mutate(mod_coef = sigmoid(abs(mod_coef_raw)),
           lvl_term = sigmoid(lvl_term_raw),
           branch_stat = w_b*max_branch_pvalue,
           tree_stat = w_t* max_tree_pvalue,
           sample_size_term = logb(n_fdcd+1,length(feature_sel_df$seqn)),
           p_y = ifelse(food_value=="food_consumed",prop_dictionary$p_x[prop_dictionary$fdcd==response],prop_dictionary$avg_fdcd[prop_dictionary$fdcd==response])) %>% 
    left_join(food_codes_all %>% 
                select(category_number,drxfcsd,ancestors) %>% 
                mutate(fdcd=trimws(as.character(format(category_number,scientific = FALSE)))) %>% 
                select(-category_number), by = "fdcd") 
  
  # determine which food codes are leaf nodes (and thus used for all non-hierarchical calculations)
  if(structure == "wweia"){
    leaf_nodes <- 
       food_codes_all %>% 
       filter(level==max_level_by_num) %>% 
       mutate(fdcd=trimws(as.character(format(category_number,scientific = FALSE)))) %>% 
       pull(fdcd)
    
  non_hic_ranks <- 
      hic_df_temp %>% 
      filter(fdcd %in% leaf_nodes) %>% 
      arrange(desc(mi_raw)) %>% 
      mutate(mi_rank = row_number()) %>% 
      arrange(desc(gr_raw)) %>% 
      mutate(gr_rank = row_number()) %>% 
      arrange(desc(su_raw)) %>% 
      mutate(su_rank = row_number()) %>% 
      arrange(desc(cor_raw)) %>% 
      mutate(cor_rank = row_number()) %>% 
      arrange(desc(test_raw)) %>% 
      mutate(test_rank = row_number()) %>% 
      arrange(desc(abs(mod_coef_raw))) %>% 
      mutate(mod_coef_rank = row_number()) %>% 
      select(fdcd, contains("rank"))
  
  hic_df_new <-
    hic_df_temp %>% 
    left_join(non_hic_ranks, by = "fdcd") %>% 
    mutate(hic_mod_coef_raw = abs(mod_coef_raw) - lvl_term - tree_stat - branch_stat + sample_size_term)
  
  } else {
    
  hic_df_new <-    
    hic_df_temp %>% 
    mutate(hic_mod_coef_raw = abs(mod_coef_raw) - lvl_term - tree_stat - branch_stat + sample_size_term)
  }
  
  hic_cols=names(hic_df_new)[str_detect(names(hic_df_new), "hic")]
  # get the ranks for hic scores.
  hic_ranks <- map(.x = hic_cols,.f = ~get_hic_rank(col=.x, hic_df = hic_df_new))
  hic_ranks_join <- reduce(hic_ranks,full_join, by = "fdcd")
  
  hic_df <-
    left_join(hic_df_new,hic_ranks_join, by = "fdcd")

  # writes the csv file for the hic values. 
  # fndds structure will only have hic_fndds. wweia will have all else.
  write_csv(hic_df, str_c("01_data/hic_results/hic_",structure,"_",response,"_",food_value,".csv"))
  
  return(hic_df)
}

# returns sigmoid of input value
sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

# Receives: food code, dataset, food_value measure being used, response of interest, and the vector of covariates
# Returns: a data frame row with all feature selection terms as well as the hic term values.
get_hic_terms_raw <- function(fdcd="95342000",feature_sel_df,food_value,response,covariates){
    print(fdcd)
  
  # Calculate feature selection score amount
  mi_raw <- information.gain(formula(str_c(response,"~`",fdcd,"`")), data =  feature_sel_df, unit = "log2")$attr_importance
  gr_raw <- gain.ratio(formula(str_c(response,"~`",fdcd,"`")), data =  feature_sel_df, unit = "log2")$attr_importance
  su_raw <- symmetrical.uncertainty(formula(str_c(response,"~`",fdcd,"`")), data =  feature_sel_df, unit = "log2")$attr_importance
  cor_raw <- ifelse(food_value=="food_consumed" & structure != "fped",
                    chi.squared(formula(str_c(response,"~`",fdcd,"`")), data =  feature_sel_df)$attr_importance,
                    abs(ltm::biserial.cor(feature_sel_df %>% pull(!!sym(fdcd)),feature_sel_df %>% pull(!!sym(response)), use = "all.obs",level=2)))
  test_raw <- ifelse(food_value=="food_consumed",get_chi(fdcd,response, feature_sel_df),get_ttest(fdcd,response,feature_sel_df))
  
  mod_formula <<- formula(str_c(response," ~ ", " . "))
  mod_lm <- glm(mod_formula, data = feature_sel_df %>% select(all_of(covariates),!!sym(fdcd),all_of(response)), family = "binomial")
  mod_coef_raw<-summary(mod_lm)$coefficients[length(summary(mod_lm)$coefficients[,1]),3]
  
  # Calculate Level Term
  if(food_value == "food_consumed"){
    lvl_term_raw <- 
    food_codes_all %>% 
      summarize(lvl_term_raw = ifelse(max_level_skips==1,1.5,logb(level+max_level_skips,max_level_skips))) %>% 
      pull(lvl_term_raw)

  prop_dict_sort <-
    prop_dictionary %>% 
    filter(fdcd!=response) %>% 
    filter(n_fdcd>0) %>% 
    arrange(fdcd) %>% 
    select(fdcd,n_fdcd,n_fdcd_y,p_xy)
  }
  
  pre_df <- 
    tibble(fdcd = fdcd,
         mi_raw = mi_raw,
         gr_raw = gr_raw,
         su_raw = su_raw,
         cor_raw = cor_raw,
         test_raw = abs(test_raw),
         mod_coef_raw = mod_coef_raw,
         n_fdcd = prop_dictionary$n_fdcd[prop_dictionary$fdcd==fdcd],
         p_x = prop_dictionary$p_x[prop_dictionary$fdcd==fdcd],
         p_y_given_x = prop_dictionary$p_y_given_x[prop_dictionary$fdcd==fdcd])
  
  # food_consumed should be done first so the proportions only need to be calculated once.
  if(food_value == "food_consumed"){
    # Calculate z-scores for each the food code passed into function with all other foods.
    pvals <- get_z_scores(fdcd,response,prop_dict_sort, food_value=food_value, feature_sel_df=feature_sel_df)
    df <- pre_df %>% 
      mutate(lvl_term_raw = lvl_term_raw,
             max_tree_pvalue=pvals$max_tree,
             max_branch_pvalue=pvals$max_branch)
  } else {
    df <- 
      pre_df %>% 
      left_join(food_consumed_results %>% 
                  select(fdcd,lvl_term_raw,max_tree_pvalue,max_branch_pvalue), by = "fdcd")
  }
  return(df)
}

# Checks to see if two foods are in the same branch.
check_branch <- function(a="10",b="1008"){
  if(as.numeric(a) > 100000){
    a_lvl <- food_codes_all$level[food_codes_all$category_number == as.numeric(a)]
    b_lvl <- food_codes_all$level[food_codes_all$category_number == as.numeric(b)]
    min <- ifelse(a_lvl<b_lvl,a_lvl,b_lvl)
    return(case_when(a_lvl==b_lvl~FALSE,
                     str_sub(a,1,min)==str_sub(b,1,min) ~ TRUE,
                     TRUE ~ FALSE))
  } else {
  b_ancestors <- as.numeric(str_split(food_codes_all$ancestors[food_codes_all$category_number==as.numeric(b)],",",simplify=TRUE)[1,])
  a_ancestors <- as.numeric(str_split(food_codes_all$ancestors[food_codes_all$category_number==as.numeric(a)],",",simplify=TRUE)[1,])
  return((as.numeric(a) %in% b_ancestors) | (as.numeric(b) %in% a_ancestors))
  }
}

# Returns the maximum p-value for a food compared to other foods in the branch and tree.
get_z_scores <- function(fdcd = "89902100", response=response,prop_dict_sort=prop_dict_sort, food_value, feature_sel_df){

  if(str_length(fdcd)<3)  print(str_c("Get Z Score for ", fdcd))
  
  branch = str_split(food_codes_all$branch[food_codes_all$category_number==fdcd],",",simplify=TRUE)[1,]
  branch_clean =branch[branch!="0" & branch != fdcd & branch %in% prop_dict_sort$fdcd]
  
  if(!is_empty(branch_clean)){
    branch_z <- sapply(branch_clean, do_test, i=fdcd, prop_dict_sort = prop_dict_sort,food_value=food_value, feature_sel_df=feature_sel_df, response=response)
    max_branch <- max(branch_z)
  } else {
    max_branch <- 1
    branch_z <- 0
  }

  tree <- prop_dict_sort$fdcd[!(prop_dict_sort$fdcd %in% branch)]
  tree_z <- sapply(tree, do_test, i=fdcd, prop_dict_sort = prop_dict_sort,food_value=food_value, feature_sel_df=feature_sel_df, response=response)
  max_tree <- max(tree_z,branch_z)
  
  return(tibble_row(max_branch = max_branch, max_tree = max_tree))
}

# Calculates statistical test for calculating z-scores for tree and branch comparison.
do_test <- function(j= "11111160",i="11111170", prop_dict_sort = prop_dict_sort, food_value, feature_sel_df, response){
  if(food_value =="food_consumed"){
    if(prop_dict_sort$n_fdcd[prop_dict_sort$fdcd==i]>0 & prop_dict_sort$n_fdcd[prop_dict_sort$fdcd==j] >0){
      test <- prop.test(x=c(prop_dict_sort$n_fdcd_y[prop_dict_sort$fdcd==i],prop_dict_sort$n_fdcd_y[prop_dict_sort$fdcd==j]), n=c(prop_dict_sort$n_fdcd[prop_dict_sort$fdcd==i],prop_dict_sort$n_fdcd[prop_dict_sort$fdcd==j]), correct=FALSE)
      return(if_else(is.na(test$p.value),1,test$p.value))
    }
  } else {
    test <- t.test(feature_sel_df %>% filter(!!sym(response)=='1') %>%  select(any_of(i)),feature_sel_df %>% filter(!!sym(response)=='1')  %>% select(any_of(j)))
    return(if_else(is.na(test$p.value),1,test$p.value))
  }
}

# Make the dictionary of p(x) probabilities.  Need to pass in filtered dataset
make_prop_dictionary <- function(df= feature_sel_df, response = "elevated_blood_pressure_risk", food_value = "food_consumed",structure){
  print("Make Prop Dictionary")
  
  # probabilities are calculated differently for binary vs continuous food measures.
  if(food_value=="food_consumed"){
    dict <-
    df %>%
       adorn_totals() %>% 
       filter(seqn=="Total") %>% 
       pivot_longer(cols = c(2:length(.)),names_to = "fdcd", values_to = "n_fdcd") %>% 
       mutate(n = length(df$seqn),
              p_x = n_fdcd/n) %>% 
       select(fdcd, n_fdcd, n, p_x) %>% 
    left_join(df %>%
                filter(!!sym(response) ==1) %>% 
                adorn_totals() %>% 
                filter(seqn=="Total") %>% 
                pivot_longer(cols = c(2:length(.)),names_to = "fdcd", values_to = "n_fdcd_y") %>% 
                mutate(n_y = n_fdcd_y[fdcd==response]) %>%
                select(fdcd, n_y, n_fdcd_y), by = c("fdcd")) %>% 
   mutate(n_not_fdcd = n-n_fdcd,
           n_fdcd_not_y = n_fdcd-n_fdcd_y,
           n_not_fdcd_not_y = (n-n_y-n_fdcd_not_y),
           n_not_fdcd_y = n_y-n_fdcd_y,
           p_xy = n_fdcd_y/n,
           p_x_not_y = (n_fdcd_not_y/n),
           p_y_not_x = (n_not_fdcd_y/n),
           p_not_x_not_y = n_not_fdcd_not_y/n,        
           p_y_given_x = n_fdcd_y/n_fdcd) 
  } else {
    dict <- 
      df %>%
      adorn_totals() %>% 
      filter(seqn=="Total") %>% 
      pivot_longer(cols = c(2:length(.)),names_to = "fdcd", values_to = "sum_fdcd") %>% 
      mutate(n = length(df$seqn),
             avg_fdcd = sum_fdcd/n) %>% 
      select(fdcd, avg_fdcd, n) %>% 
      left_join(food_consumed_results %>% 
                  select(fdcd,n_fdcd,p_x,p_y_given_x), by = "fdcd")  
    }
  return(dict)
}

# get the mutual information for the food with the response
get_mutual_information <- function(fdcd="8004",response = "waist_circ_risk", base = 2) {
  
  p_y <- prop_dictionary$p_x[prop_dictionary$fdcd==response]
  p_x <- prop_dictionary$p_x[prop_dictionary$fdcd==as.character(fdcd)]

  p_xy <- prop_dictionary$p_xy[prop_dictionary$fdcd==fdcd]
  p_x_not_y = prop_dictionary$p_x_not_y[prop_dictionary$fdcd==fdcd]
  p_y_not_x = prop_dictionary$p_y_not_x[prop_dictionary$fdcd==fdcd]
  p_not_x_not_y = prop_dictionary$p_not_x_not_y[prop_dictionary$fdcd==fdcd]

  sum1 = p_xy*log2(p_xy/(p_x*p_y)) 
  sum2 = p_y_not_x*log2(p_y_not_x/((1-p_x)*p_y)) 
  sum3 = p_x_not_y*log2(p_x_not_y/((p_x)*(1-p_y))) 
  sum4 = p_not_x_not_y*log2(p_not_x_not_y/((1-p_x)*(1-p_y)))
  sigmoid(sum1+sum2+sum3+sum4)
  return(ifelse(sum1 %in% c(Inf, -Inf, NaN),0, sum1) + ifelse(sum2 %in% c(Inf, -Inf, NaN),0, sum2) + 
           ifelse(sum3 %in% c(Inf, -Inf, NaN),0, sum3) + ifelse(sum4 %in% c(Inf, -Inf, NaN),0, sum4))
}

# get gain ratio for the food with the response
get_gain_ratio <- function(fdcd="8004", response= "at_risk", base =2){
  mi <- get_mutual_information(fdcd, response, base)
  p_x <- prop_dictionary$p_x[prop_dictionary$fdcd==as.character(fdcd)]
  p_y_given_x <- prop_dictionary$p_y_given_x[prop_dictionary$fdcd==as.character(fdcd)]
  gain_ratio <- mi/-(p_x*log2(p_x) + (1-p_x)*log2(1-p_x))  
  
  return(gain_ratio)
}

# get gini index for the food with the response
get_gini <- function(cd="8004", response= "at_risk"){
  p_x <- prop_dictionary$p_x[prop_dictionary$fdcd==as.character(cd)]

  gini <- 
    prop_dictionary %>% 
    mutate(gini = p_x*(1-((n_fdcd_y/(n_fdcd_y+n_fdcd_not_y))^2 + (n_fdcd_not_y/(n_fdcd_y + n_fdcd_not_y))^2)) +
             (1-p_x)*(1-((n_not_fdcd_y/(n_not_fdcd_y+n_not_fdcd_not_y))^2 + (n_not_fdcd_not_y/(n_not_fdcd_y + n_not_fdcd_not_y))^2))) %>% 
    filter(fdcd==cd) %>% 
    pull(gini) %>% 
    unique()
  return(gini)
}

# get chi-squared test statistic for the food with the response
get_chi <- function(cd="9", response = "at_risk", feature_sel_df){
  if(sum(feature_sel_df %>% pull(!!sym(cd))=="0")==0) {
    return(0)
  } else {
    test <- chisq.test(feature_sel_df %>% pull(!!sym(cd)), feature_sel_df %>% pull(!!sym(response)))
    return(test$statistic)  
  }
}

# get ttest for the food with the response
get_ttest <- function(cd="28345140", response = "at_risk", feature_sel_df){
  test <- t.test(feature_sel_df %>% filter(!!sym(response)=='1') %>% pull(!!sym(cd)), feature_sel_df %>% filter(!!sym(response)=='0')  %>% pull(!!sym(cd)))
  return(test$statistic)  
}

# calculate the rank for the hic term
get_hic_rank <- function(col ="hic_mod_coef_raw", hic_df=hic_df_new){
  col_name = str_c(col,"_rank")
  temp <- hic_df %>% 
    select(fdcd,!!sym(col)) %>% 
    arrange(desc(!!sym(col))) %>% 
    left_join(food_codes_all %>% 
                # filter(level==max_level_by_num) %>% 
                mutate(fdcd=trimws(as.character(format(category_number,scientific = FALSE)))) %>% 
                select(fdcd,branch), by = "fdcd")
  ranks <- tibble(fdcd=NA)
  while(length(temp$fdcd > 1)){
    ranks <- bind_rows(ranks,temp %>% select(fdcd) %>% slice_head(n=1))
    branch <- str_split(temp$branch[1],",",simplify=TRUE)[1,]
    temp <- temp[!(temp$fdcd %in% branch),]
  }
  return(ranks %>% 
    filter(complete.cases(.)) %>% 
    rownames_to_column(var = col_name) %>% 
    mutate(!!sym(col_name) := as.double(!!sym(col_name))))
}

# pulls all the rankings for each feature selection method.
get_fs_summary <- function(response=response, structure=structure, food_value=food_value, relieff = TRUE, resampling=resampling, which_cov=which_cov){
  hic_raw <- 
    read_csv(str_c("01_data/hic_results/hic_",structure,"_",response,"_",food_value,".csv"))
  
  hic <- hic_raw %>% 
    mutate(.,fdcd = trimws(format(as.numeric(fdcd),scientific=FALSE))) %>% 
           select(fdcd,drxfcsd, contains("raw"), contains("rank"),-lvl_term_raw) %>% 
           rename_with(.fn = ~(str_c(structure,"_", response, "_", food_value,"_",.)), .cols = -c(1:2))
  
  if(structure!="fndds"){
  leaf_nodes <- 
    food_codes_all %>% 
    filter(level==max_level_by_num) %>% 
    mutate(fdcd=trimws(as.character(format(category_number,scientific = FALSE)))) %>% 
    pull(fdcd)
  
  if(relieff & structure!="fndds"){
  relieff_column <- str_c(str_glue(structure,response,food_value,.sep ="_"), "_relieff")
  relieff_rank <- str_c(relieff_column,"_rank")
  relieff_summary <-
    read_csv(str_c("01_data/hic_results/relieff_",structure,"_",food_value,"_results.csv")) %>% 
    mutate(.,fdcd = trimws(format(as.numeric(fdcd),scientific=FALSE))) %>% 
    mutate({{relieff_column}} := !!sym(str_c(response, "_",structure,"_relieff"))) %>% 
    select(fdcd, as.name(relieff_column)) %>% 
    filter(fdcd %in% leaf_nodes) %>% 
    arrange(desc(!!sym(relieff_column))) %>% 
    rownames_to_column(var = str_c(relieff_column,"_rank")) %>% 
    relocate(str_c(relieff_column,"_rank"), .after = last_col()) %>% 
    mutate({{relieff_rank}} := as.numeric(!!sym(relieff_rank)))

    
  return(hic %>% left_join(relieff_summary))
  } else{
    return(hic)
  }
  }else {
    return(hic)
  }
  
}