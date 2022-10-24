### SOURCE FUNCTIONS FOR ESTIMATING DNAM FITNESS MODELS AND FITAGE ###
### Updated 10/24/2022 by Kristen McGreevy


library(glmnet)
library(glmnetUtils)
library(broom)
library(dplyr)

DNAmFitnessModels <- readRDS("DNAmFitnessModelsandFitAge_Oct2022.rds")

`%!in%` = Negate(`%in%`)
DNAmFitness_Xvars <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max", "DNAmGait_wAge", 
                       "DNAmGrip_wAge", "DNAmFEV1_wAge", "DNAmGrimAge")

FitAge_Xvars <- c("DNAmGait_noAge", "DNAmGrip_noAge", "DNAmVO2max", "DNAmGrimAge")



## DATAFRAME NEEDS TO HAVE CPGS IN COLUMNS WITH A COLUMN FOR SAMPLEID, AGE, AND FEMALE ## 
## specify dataset and id variable
data_prep2 <- function(dataset, idvariable){
  
  # keep only the CpG loci you need to calculate the DNAm fitness biomarkers
  extract_these <- colnames(dataset)[which(colnames(dataset) %in% DNAmFitnessModels$AllCpGs)]
  output_data <- dataset[, c(idvariable, "Female", "Age", extract_these)]

  # update output_data with medians if we are missing any CpGs # 
  if(length(extract_these) != length(DNAmFitnessModels$AllCpGs)){
      # separate by sex to impute medians by sex. 
      data_fem <- output_data[output_data$Female == 1, ]
      data_male <- output_data[output_data$Female == 0, ]
  
      cpgs_toadd <- colnames(DNAmFitnessModels$Female_Medians_All) %!in% extract_these
      
      # set missing CpG values to medians from our training data. 
      data_fem <- data.frame(data_fem, DNAmFitnessModels$Female_Medians_All[, cpgs_toadd])
      data_male <- data.frame(data_male, DNAmFitnessModels$Male_Medians_All[, cpgs_toadd])
  
      output_data <- rbind(data_fem, data_male)
      
      print(paste0("Total ", sum(cpgs_toadd), 
                   " Missing CpGs that are assigned median values from training data"))
  }
  return(output_data)
}


# Function to provide estimates for any 1 DNAm fitness models 
# Use this if you want to calculate individual DNAmFitness estimates and not all 6
# TidyModel is a specific model in DNAmFitnessModels list
DNAmEstimatorAnyModel <- function(dataset, TidyModel, IDvar){

  int_length <- nrow(dataset)
  mod_length <- length(TidyModel$term)
  
  Xdat <- dataset[, colnames(dataset) %in% TidyModel$term]
  Xdat <- data.frame("Intercept" = rep(1, int_length), Xdat)
  Xdatnew <- as.matrix(Xdat[, c("Intercept", TidyModel$term[2:mod_length])])
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) == mod_length-1){
     estimate <- Xdatnew %*% TidyModel$estimate
  }
  if(sum(colnames(Xdatnew)[2:mod_length] == TidyModel$term[2:mod_length]) != mod_length-1){
    print("Not All Columns in New Dataframe")
    estimate <- rep(NA, int_length)
  }
 
  # put Estimates into one dataset;
  est_data <- data.frame(DNAmEstimate = estimate, ID = dataset[, {{IDvar}}])

  return(est_data)
}


# Function to calculate all DNAm fitness estimates #
DNAmFitnessEstimators <- function(data, IDvar){
  
  data_fem <- data[data$Female ==1,]
  data_male <- data[data$Female ==0,]

  fem_est1 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_noAge_Females, IDvar) # gait without age
  fem_est2 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_noAge_Females, IDvar) # grip
  fem_est3 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar) # vo2max
  fem_est4 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Gait_wAge_Females, IDvar) # gait w age
  fem_est5 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$Grip_wAge_Females, IDvar) # grip w age
  fem_est6 <- DNAmEstimatorAnyModel(dataset = data_fem, TidyModel = DNAmFitnessModels$FEV1_wAge_Females, IDvar) # fev1 w age

  male_est1 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_noAge_Males, IDvar) # gait
  male_est2 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_noAge_Males, IDvar) # grip
  male_est3 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$VO2maxModel, IDvar) # vo2max
  male_est4 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Gait_wAge_Males, IDvar) # gait
  male_est5 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$Grip_wAge_Males, IDvar) # grip
  male_est6 <- DNAmEstimatorAnyModel(dataset = data_male, TidyModel = DNAmFitnessModels$FEV1_wAge_Males, IDvar) # fev1

  # just remerging male and female estimates into one dataset. 
  fem_est1 <- rbind(fem_est1, male_est1)
  fem_est2 <- rbind(fem_est2, male_est2)
  fem_est3 <- rbind(fem_est3, male_est3)
  fem_est4 <- rbind(fem_est4, male_est4)
  fem_est5 <- rbind(fem_est5, male_est5)
  fem_est6 <- rbind(fem_est6, male_est6)

  # reassigning values to proper names 
  fem_est1[, DNAmFitness_Xvars[1]] <- fem_est1$DNAmEstimate
  fem_est2[, DNAmFitness_Xvars[2]] <- fem_est2$DNAmEstimate
  fem_est3[, DNAmFitness_Xvars[3]] <- fem_est3$DNAmEstimate
  fem_est4[, DNAmFitness_Xvars[4]] <- fem_est4$DNAmEstimate
  fem_est5[, DNAmFitness_Xvars[5]] <- fem_est5$DNAmEstimate
  fem_est6[, DNAmFitness_Xvars[6]] <- fem_est6$DNAmEstimate

  all_ests <- Reduce(function(x,y) merge(x = x, y = y, by = "ID", all.x = TRUE, all.y = TRUE), 
                                        list(fem_est1[!is.na(fem_est1$ID), 2:3],
                                             fem_est2[!is.na(fem_est2$ID), 2:3],
                                             fem_est3[!is.na(fem_est3$ID), 2:3],
                                             fem_est4[!is.na(fem_est4$ID), 2:3],
                                             fem_est5[!is.na(fem_est5$ID), 2:3],
                                             fem_est6[!is.na(fem_est6$ID), 2:3]))
  # combine with original dataframe
  match1 <- match(data[, IDvar], all_ests$ID)
  data_and_est <- data.frame(data, all_ests[match1, ])

  return(data_and_est)
}



# Function to calculate DNAmFitAge and FitAgeAcceleration after getting 
# DNAm fitness estimates from functions above. 
# This REQUIRES Female, Age, DNAmVO2max, DNAmGrip_noAge, DNAmGait_noAge, DNAmGrimAge.
# must merge with output from above with other dataframe that has DNAmGrimAge
FitAgeEstimator <- function(data, IDvar){
  
  # prep dataset for fitage- dataset by sex
  fem <- data[data$Female == 1, c("Age", FitAge_Xvars, IDvar)]
  male <- data[data$Female == 0, c("Age", FitAge_Xvars, IDvar)]

  # can only estimate FitAge if all variables are present, remove those without them # 
  fem_comcase <- fem[complete.cases(fem), ]
  male_comcase <- male[complete.cases(male), ]
  
  # Female FitAge 
  female_fitest <- 0.1044232 * ((fem_comcase$DNAmVO2max - 46.825091) / (-0.13620215)) +
                0.1742083 * ((fem_comcase$DNAmGrip_noAge - 39.857718) / (-0.22074456)) +
                0.2278776 * ((fem_comcase$DNAmGait_noAge - 2.508547) / (-0.01245682))  +
                0.4934908 * ((fem_comcase$DNAmGrimAge - 7.978487) / (0.80928530)) 
  
  # Male FitAge
  male_fitest <- 0.1390346 * ((male_comcase$DNAmVO2max - 49.836389) / (-0.141862925)) +
                0.1787371 * ((male_comcase$DNAmGrip_noAge - 57.514016) / (-0.253179827)) +
                0.1593873 * ((male_comcase$DNAmGait_noAge - 2.349080) / (-0.009380061))  +
                0.5228411 * ((male_comcase$DNAmGrimAge - 9.549733) / (0.835120557)) 
  
  fem_X <- data.frame(fem_comcase, DNAmFitAge = female_fitest)
  male_X <- data.frame(male_comcase, DNAmFitAge = male_fitest)
  
  returned_X <- rbind(fem_X, male_X)
  returned_X$FitAgeAccel <- residuals(lm(DNAmFitAge ~ Age, data = returned_X))
  
  # Calculate Mean Absolute Deviation (MAD) of DNAmFitAge from Age
  Female_meanAbsDev <- mean(abs(fem_comcase$Age - female_fitest), na.rm=TRUE)
  print(paste0("Female MAD: ", round(Female_meanAbsDev, 2)))
  
  Male_meanAbsDev <- mean(abs(male_comcase$Age - male_fitest), na.rm=TRUE)
  print(paste0("Male MAD: ", round(Male_meanAbsDev, 2)))

  return(returned_X)
}




