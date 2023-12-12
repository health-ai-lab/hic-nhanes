# NHANES HIC

This readme file explains the files and how to effectively use them. You can reach me at jpleuss@stevens.edu for any questions.

## Data Structures
demographics_with_responses_prefiltered.csv - This csv contains all the participants with their health outcomes and covariate information.

wweia_code_structure.csv - This file contains the structure for the WWEIA food hierarchy.

fndds_code_structure.csv - This file contains the structure for the FNDDS food hierarchy.

The data with all the food information (called "seqn_fdcd_{food consumption measure}.csv was too large to upload here. It can be obtained from online at NHANES for each cycle desired (https://wwwn.cdc.gov/Nchs/Nhanes/search/datapage.aspx?Component=Dietary&CycleBeginYear=2007 as an example for 2007-2008). The structure ends up an nxm matrix with n rows of participants and m columns of food codes.

## read_in_data_of_interest_clean.R
This is a helper script that pulls in all the NHANES data, but requires all the NHANES data to be saved in 

## data_prep_clean.R
This script imports and joins the required data, filters the participants based on the desired health outcome, and creates the test/train split.

## hic_code_clean.R
This script generates the HIC score and ranking for the desired food consumption measure and health outcome. It also generates the ranking for all the other feature selection rankings.

## run_model_clean.R
This is the main model running script. It uses the other scripts to collect the data, then runs the selected model with the designated number of food features for each for each of the feature selection algorithms. It saves the model classification output in a model registry file if desired.
