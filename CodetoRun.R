### Run this code after downloading all files into the same directory ###

# make sure these packages are installed and loaded #
library(glmnet)
library(glmnetUtils)
library(broom)
library(dplyr)


# This source call will make 5 functions: %!in%, data_prep2, DNAmEstimatorAnyModel, 
# DNAmFitnessEstimators, and FitAgeEstimator. 
# It will also make 2 global variables: DNAmFitness_Xvars and FitAge_Xvars
source("DNAmFitnessandDNAmFitAgeSourceFunctions.R")


### data_prep2 ###
# prepares your data frame for estimation. This data frame must have 
# at least an identification column, Age, and Female. In order to get proper estimates 
# based on the DNAm values, the data frame must also have CpGs as columns. 


### DNAmFitnessEstimators ###
# estimates the DNAm fitness biomarkers. This takes the output from data_prep2 
# Specify the identification column. If you have only males or only females in your 
# dataframe, use DNAmEstimatorAnyModel for one model at a time. 


### FitAgeEstimator ###
# calculates DNAmFitAge and FitAgeAcceleration. This data frame MUST contain 
# ID variable, Female, Age, DNAmGait_noAge, DNAmGrip_noAge, DNAmVO2max, and DNAmGrimAge. 
# To calculate DNAmGrimAge, use the online calculator found at http://dnamage.genetics.ucla.edu/ 


### DNAmEstimatorAnyModel ###
# Function to provide estimates for any 1 DNAm fitness model. Use this if you 
# want to calculate individual DNAmFitness estimates and not all 6. For example, 
# this function will need to be used when you have only males or only females in 
# your dataframe. 


##################### EXAMPLE ######################
# This example uses a small simulated dataset. 
Sample_Sim_Data <- readRDS(file = "Sample_Simulation_Data.rds")

### data_prep2 example with 2 CpG columns missing in simulation data ###
sample_data_prep <- data_prep2(dataset = Sample_Sim_Data, idvariable = "SampleID")
# notice how the function prints that 2 CpG columns were missing and the median 
# values were used. 


### DNAmFitnessEstimators example ###
sample_data_FitnessEst <- DNAmFitnessEstimators(sample_data_prep, IDvar = "SampleID")


### FitAgeEstimator example ###

# First need to merge output from function above with dataframe that has DNAmGrimAge in it. 
sample_data_FitAge_prep <- merge(Sample_Sim_Data[, c("SampleID", "DNAmGrimAge")], 
                                 sample_data_FitnessEst, by = "SampleID")

FitAge_out <- FitAgeEstimator(sample_data_FitAge_prep, IDvar = "SampleID")
# notice that the Mean Absolute Deviation of FitAge to Chronological Age for each sex is printed 

colnames(FitAge_out)
# Each DNAm fitness biomarker, DNAmFitAge, and FitAgeAcceleration are in this output. 



