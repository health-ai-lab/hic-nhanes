library(tidyverse)
library(fst)
library(janitor)

# pulls in all demographic data to useable environment.
pull_dem_data <- function(i){
  
 ## demographics
 
  dem_var_wt_file <- ifelse(i %in% c(2007,2009,2011), 
                            str_c("01_data/", i, "/Demographics/demographic_variables_sample_weights_",i,".fst"),
                            str_c("01_data/", i, "/Demographics/demographic_variables_and_sample_weights_",i,".fst"))
  assign(str_c("demographic_variables_sample_weights_", i),
        read_fst(dem_var_wt_file) %>%
          clean_names()  %>%
          mutate(year=i) %>% 
          as_tibble())
  if(!exists('demographic_variables_sample_weights')){
    assign('demographic_variables_sample_weights', get(str_c("demographic_variables_sample_weights_", i)), envir = globalenv())
  } else {
    demographic_variables_sample_weights <<- bind_rows(demographic_variables_sample_weights, get(str_c("demographic_variables_sample_weights_", i)))
  }
}

# pulls in all health data to useable environment.
pull_health_data <- function(i){
  
  ## blood pressure
  assign(str_c("blood_pressure_", i), 
         read_fst(str_c("01_data/", i, "/Examination/blood_pressure_", i, ".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('blood_pressure')){
    assign('blood_pressure', get(str_c("blood_pressure_", i)), envir = globalenv())
  } else {
    blood_pressure <<- bind_rows(blood_pressure, get(str_c("blood_pressure_", i)))
  }
  
  ## body measures
  assign(str_c("body_measures_", i), 
         read_fst(str_c("01_data/", i, "/Examination/body_measures_", i, ".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('body_measures')){
    assign('body_measures', get(str_c("body_measures_", i)), envir = globalenv())
  } else {
    body_measures <<- bind_rows(body_measures, get(str_c("body_measures_", i)))
  }
  
  ## labratory cholesterol hdl
  hdl_file <- ifelse(i %in% c(2007,2009,2011,2013), 
                     str_c("01_data/", i, "/Laboratory/cholesterol_hdl_",i,".fst"),
                     str_c("01_data/", i, "/Laboratory/cholesterol_high_density_lipoprotein_hdl_",i,".fst"))
  
  assign(str_c("cholesterol_hdl_",i), 
         read_fst(hdl_file) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('cholesterol_hdl')){
    assign('cholesterol_hdl', get(str_c("cholesterol_hdl_", i)), envir = globalenv())
  } else {
    cholesterol_hdl <<- bind_rows(cholesterol_hdl, get(str_c("cholesterol_hdl_", i)))
  }
  
  ## labratory cholesterol ldl
  # i=2017
  ldl_file <- ifelse(i %in% c(2007,2009,2011,2013), 
                     str_c("01_data/", i, "/Laboratory/cholesterol_ldl_triglycerides_",i,".fst"),
                     ifelse(i ==2015,
                        str_c("01_data/", i, "/Laboratory/cholesterol_low_density_lipoprotein_ldl_triglycerides_",i,".fst"),
                        str_c("01_data/", i, "/Laboratory/cholesterol_low_density_lipoproteins_ldl_triglycerides_",i,".fst"))
  )
  assign(str_c("cholesterol_ldl_",i), 
         read_fst(ldl_file) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('cholesterol_ldl')){
    assign('cholesterol_ldl', get(str_c("cholesterol_ldl_", i)), envir = globalenv())
  } else {
    cholesterol_ldl <<- bind_rows(cholesterol_ldl, get(str_c("cholesterol_ldl_", i)))
  }
  
  ## labratory cholesterol total
  assign(str_c("cholesterol_total_",i),
         read_fst(str_c("01_data/",i,"/Laboratory/cholesterol_total_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('cholesterol_total')){
    assign('cholesterol_total', get(str_c("cholesterol_total_", i)), envir = globalenv())
  } else {
    cholesterol_total <<- bind_rows(cholesterol_total, get(str_c("cholesterol_total_", i)))
  }
  
  ## labratory glycohemoglobin
  assign(str_c("glycohemoglobin_",i),
         read_fst(str_c("01_data/",i,"/Laboratory/glycohemoglobin_",i,".fst")) %>%
                    clean_names()  %>%
                    as_tibble())
  if(!exists('glycohemoglobin')){
    assign('glycohemoglobin', get(str_c("glycohemoglobin_", i)), envir = globalenv())
  } else {
    glycohemoglobin <<- bind_rows(glycohemoglobin, get(str_c("glycohemoglobin_", i)))
  }
  
  ## labratory oral glucose (No OGTT given in 2017)
  if(i!=2017) {assign(str_c("oral_glucose_tolerance_test_",i),
        read_fst(str_c("01_data/",i,"/Laboratory/oral_glucose_tolerance_test_",i,".fst")) %>%
          clean_names()  %>%
          mutate(year=i) %>% 
          as_tibble())
    if(!exists('oral_glucose_tolerance_test')){
      assign('oral_glucose_tolerance_test', get(str_c("oral_glucose_tolerance_test_", i)), envir = globalenv())
    } else {
      oral_glucose_tolerance_test <<- bind_rows(oral_glucose_tolerance_test, get(str_c("oral_glucose_tolerance_test_", i)))
    }
  }
   
  ## labratory plasma fasting
  plasma_glucose_file <- ifelse(i %in% c(2007,2009,2011), 
                     str_c("01_data/", i, "/Laboratory/plasma_fasting_glucose_insulin_",i,".fst"),
                     str_c("01_data/", i, "/Laboratory/plasma_fasting_glucose_",i,".fst"))
  
  assign(str_c("plasma_fasting_glucose_insulin_",i),
         read_fst(plasma_glucose_file) %>%
          clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('plasma_fasting_glucose_insulin')){
    assign('plasma_fasting_glucose_insulin', get(str_c("plasma_fasting_glucose_insulin_", i)), envir = globalenv())
  } else {
    plasma_fasting_glucose_insulin <<- bind_rows(plasma_fasting_glucose_insulin, get(str_c("plasma_fasting_glucose_insulin_", i)))
  }
  
  # i=2017
  ## labratory pregnancy test
  pregnancy_test_file <- ifelse(i %in% c(2017), 
                                str_c("01_data/", i, "/Laboratory/urine_pregnancy_test_",i,".fst"),
                                str_c("01_data/", i, "/Laboratory/pregnancy_test_urine_",i,".fst"))
  
  assign(str_c("pregnancy_",i),
         read_fst(pregnancy_test_file) %>%
           clean_names()  %>%
           as_tibble())
  if(!exists('pregnancy')){
    assign('pregnancy', get(str_c("pregnancy_", i)), envir = globalenv())
  } else {
    pregnancy <<- bind_rows(pregnancy, get(str_c("pregnancy_", i)))
  }
  
  ## questionnaire physical activity
  assign(str_c("physical_activity_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/physical_activity_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('physical_activity')){
    assign('physical_activity', get(str_c("physical_activity_", i)), envir = globalenv())
  } else {
    physical_activity <<- bind_rows(physical_activity, get(str_c("physical_activity_", i)))
  }
  
  ## questionnaire blood pressure
  assign(str_c("blood_pressure_cholesterol_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/blood_pressure_cholesterol_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('blood_pressure_cholesterol')){
    assign('blood_pressure_cholesterol', get(str_c("blood_pressure_cholesterol_", i)), envir = globalenv())
  } else {
    blood_pressure_cholesterol <<- bind_rows(blood_pressure_cholesterol, get(str_c("blood_pressure_cholesterol_", i)))
  }
  
  ## questionnaire cardiovascular health
  assign(str_c("cardiovascular_health_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/cardiovascular_health_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('cardiovascular_health')){
    assign('cardiovascular_health', get(str_c("cardiovascular_health_", i)), envir = globalenv())
  } else {
    cardiovascular_health <<- bind_rows(cardiovascular_health, get(str_c("cardiovascular_health_", i)))
  }
  
  ## questionnaire diabetes
  assign(str_c("diabetes_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/diabetes_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('diabetes')){
    assign('diabetes', get(str_c("diabetes_", i)), envir = globalenv())
  } else {
    diabetes <<- bind_rows(diabetes, get(str_c("diabetes_", i)))
  }
  
  ## questionnaire physical functioning
  assign(str_c("physical_functioning_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/physical_functioning_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('physical_functioning')){
    assign('physical_functioning', get(str_c("physical_functioning_", i)), envir = globalenv())
  } else {
    physical_functioning <<- bind_rows(physical_functioning, get(str_c("physical_functioning_", i)))
  }
  
  ## Alcohol questionnaire
  assign(str_c("alcohol_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/alcohol_use_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('alcohol')){
    assign('alcohol', get(str_c("alcohol_", i)), envir = globalenv())
  } else {
    alcohol <<- bind_rows(alcohol, get(str_c("alcohol_", i)))
  }  
  
  ## Smoking questionare
  assign(str_c("smoking_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/smoking_cigarette_use_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('smoking')){
    assign('smoking', get(str_c("smoking_", i)), envir = globalenv())
  } else {
    smoking <<- bind_rows(smoking, get(str_c("smoking_", i)))
  }
  
  ## Reproductive Health questionnaire
  assign(str_c("reproductive_health_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/reproductive_health_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('reproductive_health')){
    assign('reproductive_health', get(str_c("reproductive_health_", i)), envir = globalenv())
  } else {
    reproductive_health <<- bind_rows(reproductive_health, get(str_c("reproductive_health_", i)))
  }

## Income questionnaire
  assign(str_c("income_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/income_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('income')){
    assign('income', get(str_c("income_", i)), envir = globalenv())
  } else {
    income <<- bind_rows(income, get(str_c("income_", i)))
  }
  
  ## Sleep questionnaire
  assign(str_c("sleep_disorders_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/sleep_disorders_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('sleep_disorders')){
    assign('sleep_disorders', get(str_c("sleep_disorders_", i)), envir = globalenv())
  } else {
    sleep_disorders <<- bind_rows(sleep_disorders, get(str_c("sleep_disorders_", i)))
  }
  
  ## Current Health questionnaire
  assign(str_c("current_health_status_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/current_health_status_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('current_health_status')){
    assign('current_health_status', get(str_c("current_health_status_", i)), envir = globalenv())
  } else {
    current_health_status <<- bind_rows(current_health_status, get(str_c("current_health_status_", i)))
  }
  
  ## Food Security questionnaire
  assign(str_c("food_security_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/food_security_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('food_security')){
    assign('food_security', get(str_c("food_security_", i)), envir = globalenv())
  } else {
    food_security <<- bind_rows(food_security, get(str_c("food_security_", i)))
  }
  
  ## Consumer Behavior questionnaire
  assign(str_c("consumer_behavior_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/consumer_behavior_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('consumer_behavior')){
    assign('consumer_behavior', get(str_c("consumer_behavior_", i)), envir = globalenv())
  } else {
    consumer_behavior <<- bind_rows(consumer_behavior, get(str_c("consumer_behavior_", i)))
  }
  
  ## Diet Behavior questionnaire
  assign(str_c("diet_behavior_nutrition_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/diet_behavior_nutrition_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('diet_behavior_nutrition')){
    assign('diet_behavior_nutrition', get(str_c("diet_behavior_nutrition_", i)), envir = globalenv())
  } else {
    diet_behavior_nutrition <<- bind_rows(diet_behavior_nutrition, get(str_c("diet_behavior_nutrition_", i)))
  }
  
  ## Weight History questionnaire
  assign(str_c("weight_history_",i),
         read_fst(str_c("01_data/",i,"/Questionnaire/weight_history_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('weight_history')){
    assign('weight_history', get(str_c("weight_history_", i)), envir = globalenv())
  } else {
    weight_history <<- bind_rows(weight_history, get(str_c("weight_history_", i)))
  }
  
}

# pulls in all food data to useable environment.
pull_food_data <- function(i=2011){
  
  ## dietary individual foods 1
  assign(str_c("food_day1pre_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_interview_individual_foods_first_day_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('food_day1pre')){
    assign('food_day1pre', get(str_c("food_day1pre_", i)), envir = globalenv())
  } else {
    food_day1pre <<- bind_rows(food_day1pre, get(str_c("food_day1pre_", i)))
  }
  
  ## dietary individual foods 2
  assign(str_c("food_day2pre_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_interview_individual_foods_second_day_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('food_day2pre')){
    assign('food_day2pre', get(str_c("food_day2pre_", i)), envir = globalenv())
  } else {
    food_day2pre <<- bind_rows(food_day2pre, get(str_c("food_day2pre_", i)))
  }
  
  ## dietary foods codes
  assign(str_c("food_codes_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_interview_technical_support_file_food_codes_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('food_codes')){
    assign('food_codes', get(str_c("food_codes_", i)), envir = globalenv())
  } else {
    food_codes <<- bind_rows(food_codes, get(str_c("food_codes_", i)))
  }
  
  ## dietary total nutrients day 1
  assign(str_c("total_nutrients_day1_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_interview_total_nutrient_intakes_first_day_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('total_nutrients_day1')){
    assign('total_nutrients_day1', get(str_c("total_nutrients_day1_", i)), envir = globalenv())
  } else {
    total_nutrients_day1 <<- bind_rows(total_nutrients_day1, get(str_c("total_nutrients_day1_", i)))
  }
  
  ## dietary total nutrients day 1
  assign(str_c("total_nutrients_day2_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_interview_total_nutrient_intakes_second_day_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('total_nutrients_day2')){
    assign('total_nutrients_day2', get(str_c("total_nutrients_day2_", i)), envir = globalenv())
  } else {
    total_nutrients_day2 <<- bind_rows(total_nutrients_day2, get(str_c("total_nutrients_day2_", i)))
  }
  
  ## dietary total supplements day 1
  assign(str_c("total_supplements_day1_",i),
         read_fst(str_c("01_data/",i,"/Dietary/dietary_supplement_use_24_hour_total_dietary_supplements_first_day_",i,".fst")) %>%
           clean_names()  %>%
           mutate(year=i) %>% 
           as_tibble())
  if(!exists('total_supplements_day1')){
    assign('total_supplements_day1', get(str_c("total_supplements_day1_", i)), envir = globalenv())
  } else {
    total_supplements_day1 <<- bind_rows(total_supplements_day1, get(str_c("total_supplements_day1_", i)))
  }
  
  
  
  years <- c(7,9,11,13,15,17)
  get_fndds_map <- function(year){
    return(readxl::read_xlsx(case_when(
      year==7 ~ str_c("01_data/",year+2000,"/Dietary/WWEIA0708_foodcat_FNDDS.xlsx"),
      year==9 ~ str_c("01_data/",year+2000,"/Dietary/WWEIA0910_foodcat_FNDDS.xlsx"),
      TRUE ~ str_c("01_data/",year+2000,"/Dietary/WWEIA",year,year+1,"_foodcat_FNDDS.xlsx")), sheet = 1) %>%
        mutate(year = year+2000))
  }
  
  years <- c(7,9,11,13,15,17)
  get_fped_map <- function(year){
    return(readxl::read_xls(case_when(
      year==7 ~ str_c("01_data/",year+2000,"/FPED_0708.xls"),
      year==9 ~ str_c("01_data/",year+2000,"/FPED_0910.xls"),
      TRUE ~ str_c("01_data/",year+2000,"/FPED_",year,year+1,".xls")), sheet = 1) %>%
        mutate(year = year+2000))
  }
  

  most_recent_map <<-
    bind_rows(lapply(years,get_fndds_map)) %>% 
    filter(!is.na(category_number)) %>%
    select(food_code,category_number, category_description, year) %>% 
    group_by(food_code) %>%
    slice_max(year) %>% 
    ungroup() %>% 
    select(-year)
  

  food_day1 <<-
    food_day1pre %>% #count(year)
    left_join(most_recent_map, by = c("dr1ifdcd"="food_code"))
  # remove(food_day1pre)
  
  food_day2 <<-
    food_day2pre %>% 
    left_join(most_recent_map, by = c("dr2ifdcd"="food_code"))
  # remove(food_day2pre)
}

years <- c(7,9,11,13,15,17)


if(update_data) {
  map(c(2007,2009,2011,2013,2015,2017),pull_dem_data)
  map(c(2007,2009,2011,2013,2015,2017),pull_health_data)
  map(c(2007,2009,2011,2013,2015,2017),pull_food_data)
  remove(food_day1pre)
  remove(food_day2pre)
}
