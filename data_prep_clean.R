update_data=FALSE
library(tidyverse)
library(tidymodels)
# source("read_in_data_of_interest.R")
# seqn_fdcd_df <- read_rds("01_data/seqn_fdcd.rds")

####### Responses
# NHANES Data Frames Required: diabetes, glycohemoglobin, plasma_fasting_glucose_insulin, oral_glucose_tolerance_test,
#                       demographic_variables, pregnancy, blood pressure, body measures, cholesterol,
#                       physical activity, cardiovascular health, physical health, smoking, income, reproductive health,
#                       weight history, diet behavior, 

# Joins all the required NHANES data and creates the features necessary
prep_data <- function(){
  if(update_data){
    demographics_of_interest <- 
    demographic_variables_sample_weights %>% 
      select(seqn, 
             ridageyr, # Age in years of the participant at the time of screening
             riagendr, # Sex of participant
             ridexprg,  # pregnancy status
             ridreth1, # race ethnicity
             ridreth3, #race ethnicity
             wtmec2yr, #sample weights
             sdmvstra, #sample strata
             sdmvpsu, #psu cluster
             dmdeduc2, # education level
             indfmpir, # poverty income ratio
             year
      ) %>% 
      mutate(education = factor(case_when(
        dmdeduc2 %in% c(1,2) ~ "less_high_school",
        dmdeduc2 == 3 ~ "high_school",
        dmdeduc2 %in% c(4,5) ~ "more_high_school"
      ), levels = c("less_high_school","high_school","more_high_school"),ordered=T)) %>% 
      mutate(race_ethnicity = case_when(
        ridreth3 == 1 | ridreth1 ==1 ~ "Mexican American",
        ridreth3 == 2 | ridreth1 ==2 ~ "Other Hispanic",
        ridreth3 == 3 | ridreth1 ==3 ~ "Non-Hispanic White",
        ridreth3 == 4 | ridreth1 ==4 ~ "Non-Hispanic Black",
        ridreth3 == 6 ~ "Non-Hispanic Asian",
        ridreth3 == 7 | ridreth1 ==5~ "Other Race"
      )) %>% 
      left_join(
        body_measures %>% 
          select(seqn,
                 bmxbmi, # Body Mass Index (kg/m**2)
                 bmxht, # Standing Height (cm)
                 bmxwt, # Weight (kg)
                 bmxwaist # waist circumference
          )
      ) %>% 
      left_join(
        smoking %>% 
          select(seqn,
                 smq020, #smoked 100 cigarettes in life?
                 smq040  #currently smoking status?
          ) %>% 
          mutate(smoking = case_when(
            smq020 == 1 & smq040 %in% c(1,2) ~ "current",
            smq020 == 1 & smq040 == 3 ~ "former",
            smq020 == 2 ~ "never"
          )) %>% 
          mutate(smoking = factor(smoking, levels =  c("never","former","current"), ordered=T)) %>% 
          mutate(smq020 = case_when(smq020 == 1 ~ "Yes",
                                    smq020 == 2 ~ "No",
                                    TRUE ~ "NA")) %>% 
          mutate(smq020 = ifelse(smq020 == "NA", NA, smq020)) 
      ) %>%
      mutate(riagendr = if_else(riagendr == 1, "male", "female")) %>% 
      mutate(ridexprg = ifelse(ridexprg == 1, "1", "0")) %>% 
      mutate(ridexprg = ifelse(is.na(ridexprg), "0", ridexprg)) 
    
  demographics_with_responses <-
    list(
      cholesterol_ldl %>% select(seqn,
                                 lbxtr, # triglyceride mg/dl
                                 wtsaf2yr), #fasting subsample weights
      cholesterol_hdl %>%   select(seqn,
                                   lbdhdd), # direct hdl cholesterol mg/dl
      blood_pressure %>%  select(seqn, 
                                 contains("bpxdi")), # bpxdi1 ... bpxdi4 up to four diastolic blood pressure readings 
      blood_pressure %>%  select(seqn, 
                                 contains("bpxsy")), # bpxsy1 ... bpxsy4 up to four systolic blood pressure readings
      plasma_fasting_glucose_insulin %>% select(seqn,
                                                # wtsaf2yr, # weights for fasting
                                                lbxglu), # fasting plasma glucose mg/dl
      blood_pressure_cholesterol %>% select(seqn,
                                            bpq020, # Ever told you have high blood pressure
                                            bpq040a, # taking prescription for hypertension
                                            bpq050a, # now taking Blood pressure medication
                                            bpq080, # Doctor told you have high cholesterol
                                            bpq100d) , # now taking high cholesterol medication
      diabetes %>% select(seqn,
                          diq010, # Doctor told you have diabetes
                          diq070, # Taking Diabetic Pills
                          diq170, # told at risk for diabetes
                          diq172, # feel at risk for diabetes
                          diq175a, # family history of diabetes
                          diq175s, # gestational diabetes
                          diq050), # Taking insulin now
      glycohemoglobin %>% select(seqn,
                                 lbxgh), # hba1c %
      oral_glucose_tolerance_test %>% select(seqn,
                                             lbxglt), # OGTT mg/dl
      pregnancy %>% select(seqn,
                           urxpreg), # Pregnancy urine test
      physical_activity %>% 
                          select(seqn, 
                                 paq605, # Vigorous work-related activity (8 points)
                                 pad615, # minutes of Vigorous work-related activity (8 points)
                                 paq620, # Moderate work activity (4 points)
                                 pad630, # minutes of Moderate work activity (4 points)
                                 paq635, # walk or bike for transpo (4 point)
                                 pad645, # minutes of walk or bike for transpo (4 point)
                                 paq650, # vigorous leisure-time activity (8 points)
                                 pad660, # minutes of vigorous leisure-time activity (8 points)
                                 paq665, # moderate leisure-time activity (4 points)
                                 pad675), # minutes of moderate leisure-time activity (4 points)
      alcohol %>% 
                          select(seqn,
                                 alq130), # Num alcoholic drinks/day
      reproductive_health %>% select(seqn,
                              rhq162), # During Pregnancy, told you have diabetes (Yes - 1)
      total_nutrients_day1 %>%  select(seqn,
                            wtdrd1, # Sample weights for day 1 of food logs
                            dr1tnumf, # total number of foods day 1
                            dr1tkcal,# total calories day 1
                            drqsdiet, # on a special diet?
                            dr1_300), # Compare food to normal (1 - Much more; 2 - Usual; 3 - Much less; 7 - Refused; 9 - Unknown)
      total_nutrients_day2 %>% select(seqn,
                                      wtdr2d, # Sample weights for day 2 of food logs
                                      dr2tnumf, # total number of foods day 2
                                      dr2tkcal, # total calories day 2
                                      dr2_300), # Compare food to normal (1 - Much more; 2 - Usual; 3 - Much less; 7 - Refused; 9 - Unknown)
      food_day1 %>%  select(seqn,
                            dr1day) %>% # day of week of food recall 1
        unique(),
      food_day2 %>% select(seqn,
                           drdint,    # Number of days of food intake
                           dr2day) %>%  # day of week of food recall 2
        unique()
    )  %>% 
    reduce(
      left_join, .init = demographics_of_interest, by = "seqn"
    ) %>% 
    mutate(pa_mets = paq605*8 + paq620*4 + paq635*4 + paq650*8 + paq665*4) %>% # establish the physical activity measure
    mutate(bpq050a = ifelse(is.na(bpq050a),0,1), bpq100d = ifelse(is.na(bpq100d),0,1)) %>% 
    mutate(triglycerides_risk = factor(ifelse(lbxtr > 150,1,0), levels = c(0,1))) %>% 
    mutate(reduced_HDL_risk = factor(ifelse((lbdhdd < 40 & riagendr == "male") | (lbdhdd < 50 & riagendr == "female"), 1,0), levels = c(0,1))) %>% 
    mutate(bpxsy = rowMeans(subset(., select=c(bpxsy1, bpxsy2, bpxsy3, bpxsy4)), na.rm=T)) %>% 
    mutate(bpxdi = rowMeans(subset(., select=c(bpxdi1, bpxdi2, bpxdi3, bpxdi4)), na.rm=T)) %>% 
    mutate(elevated_blood_pressure_risk = factor(ifelse(bpxsy >= 130 | bpxdi >= 85,1,0), levels = c(0,1))) %>% 
    mutate(elevated_fasting_glucose_risk = factor(ifelse(lbxglu > 100,1,0), levels = c(0,1))) %>% 
    mutate(waist_circ_risk = factor(ifelse(ifelse(riagendr=='male',bmxwaist >=102,bmxwaist >=88),1,0), levels = c(0,1))) %>% 
    rowwise() %>% 
  
    mutate(diabetic = ifelse((replace(as.numeric(lbxgh >= 6.5),is.na(lbxgh),0) + replace(as.numeric(lbxglu >= 126),is.na(lbxglu),0) + replace(as.numeric(lbxglt >= 200),is.na(lbxglt),0))>0,1,0),
           diabetic = factor(ifelse(as.numeric(!is.na(lbxgh)) + as.numeric(!is.na(lbxglu)) + as.numeric(!is.na(lbxglt))>0,diabetic,NA), levels = c(0,1)),
           prediabetic = ifelse(((replace(as.numeric(between(lbxgh, 5.7, 6.4)),is.na(lbxgh),0) + replace(as.numeric(between(lbxglu, 100, 125)),is.na(lbxglu),0) + replace(as.numeric(between(lbxglt, 140, 199)),is.na(lbxglt),0))>0) & diabetic!=1,1,0),
            prediabetic = factor(ifelse(as.numeric(!is.na(lbxgh)) + as.numeric(!is.na(lbxglu)) + as.numeric(!is.na(lbxglt)) > 0, prediabetic, NA), levels = c(0,1))) %>% 
    mutate(risk_factors = sum(as.numeric(as.character(triglycerides_risk)), as.numeric(as.character(reduced_HDL_risk)), as.numeric(as.character(elevated_blood_pressure_risk)), as.numeric(as.character(elevated_fasting_glucose_risk)),as.numeric(as.character(waist_circ_risk)))) %>% 
    mutate(risk_factors_collapsed = case_when(risk_factors == 0 ~ "None",
                                              risk_factors <= 2 ~ "1 or 2",
                                              risk_factors <= 5 ~ ">= 3")) %>% 
    mutate(risk_factors_collapsed = factor(risk_factors_collapsed, levels=c("None", "1 or 2", ">= 3"), ordered=T)) %>% 
    mutate(alcohol = factor(case_when(alq130 > 2 & riagendr == "male" ~ "high",
                               alq130 > 1 & riagendr == "female" ~ "high",
                               alq130 <= 2 & riagendr=="male" ~ "moderate",
                               alq130 <= 1 & riagendr=="female" ~ "moderate",
                               is.na(alq130) ~ "low"),
                            levels = c("low","moderate","high"), ordered=T)) %>% 
    ungroup() %>% 
    mutate(at_risk = factor(ifelse(risk_factors > 2, 1, 0), levels = c(0,1))) %>% 
    mutate(pregnant = ifelse(ridexprg==1 | urxpreg==1, 1,0),
           pregnant = ifelse(is.na(pregnant),0,pregnant)) %>% # 50,588 total
    mutate(day1_wkend = ifelse(dr1day %in% c(1,6,7),1,0),
           day1_wkday = ifelse(dr1day %in% c(2:5),1,0),
           day2_wkend = ifelse(dr2day %in% c(1,6,7),1,0),
           day2_wkday = ifelse(dr2day %in% c(2:5),1,0),
           num_wkend = day1_wkend + day2_wkend,
           num_wkday = day1_wkday + day2_wkday) %>% 
    mutate(family_history = case_when(
      diq175a == 10 ~ "yes",
      diq172 == 2 ~ "no",
      diq172 == 1 & diq175a != 10 ~ "no",
    TRUE ~ "unknown")) %>% 
    mutate(sub_weights = ifelse(is.na(wtsaf2yr)|wtsaf2yr==0,wtmec2yr/2,wtsaf2yr/2)) %>% 
    select(-day1_wkend,-day1_wkday,-day2_wkend,-day2_wkday)
  
  demographics_with_responses_prefiltered <-
    demographics_with_responses %>% 
      filter(ridageyr >= 20) %>% # 34,770 total
      filter(pregnant == 0) %>% # 34,398 total
      filter(!is.na(wtdrd1)) %>% # 31,948 total
      filter(!is.na(wtdr2d)) %>%
      filter(wtdr2d>0 & dr1tkcal >0 & dr2tkcal>0) %>% #26,511
      mutate(physical_activity = factor(ntile(pa_mets,3),ordered=TRUE, levels = c(1,2,3))) 

  write_rds(demographics_with_responses_prefiltered, "01_data/demographics_with_responses_prefiltered.rds")
  } else {
    demographics_with_responses_prefiltered <- read_rds("01_data/demographics_with_responses_prefiltered.rds")
  }
  return(demographics_with_responses_prefiltered)
}
  
# gathers all the code structures and NHANES data and creates test/train sets
get_data <- function(response = "elevated_blood_pressure_risk", structure = "wweia", food_value="food_consumed"){
  
  # Pull in correct food code structure
  if(structure=="wweia") {
    food_codes_all <<- read_rds("01_data/wweia_code_structure.rds") #%>% 
   } else if(structure=="fndds") {
    food_codes_all <<- read_rds("01_data/fndds_code_structure.rds")
   }  else {
    print("wrong structure name")
  }
  
  # Get correct seqn_fdcd_lists.  these files are structures 

  seqn_fdcd_df <<- read_rds(str_c("01_data/",structure,"_seqn_fdcd_",food_value,".rds"))
  
  demo_pre_filtered <- prep_data() 
  
  # create the test/train sets
    set.seed(123)
  data_split <<- initial_split(demo_pre_filtered, prop=.7)
  train_data <<- filter_data("train",response)
  test_data <<- filter_data("test",response)
  
}  

# Filters test and train sets based on response requirements
filter_data <- function(split= "train",response){
 
if(split=="train") original <- training(data_split)
if(split=="test") original <- testing(data_split)
  
  filtered_data <-
    original %>% 
    filter(case_when(response=="at_risk" ~ !is.na(at_risk),
                     response=="waist_circ_risk" ~ !is.na(waist_circ_risk),
                     response=="elevated_fasting_glucose_risk" ~ !is.na(elevated_fasting_glucose_risk) &
                       (is.na(diq010) | diq010 != 1) &
                       (is.na(diq175s) | diq175s!=28) & 
                       (is.na(rhq162) | rhq162 != 1) &
                       (is.na(diq050) | diq050 != 1) &
                       (is.na(diq070) | diq070 != 1),
                     response=="elevated_blood_pressure_risk" ~ !is.na(elevated_blood_pressure_risk) & 
                       (is.na(bpq050a) | bpq050a != 1),
                     response=="reduced_HDL_risk" ~ !is.na(reduced_HDL_risk) & 
                       (is.na(bpq100d) | bpq100d != 1),
                     response=="triglycerides_risk" ~ !is.na(triglycerides_risk) & 
                       (is.na(bpq100d) | bpq100d != 1)
    ))
  return(filtered_data)
}
### Sample Sizes (train/test/total)
# waist_circ_risk - 17900/7682/25582
# at_risk - 8194/3399/11593
# elevated_fasting_glucose_risk - 7275/3057/10332
# elevated_blood_pressure_risk - 12048/5227/17275
# reduced_HDL_risk - 13006/5647/18653
# triglycerides_risk - 6291/2638/8929
