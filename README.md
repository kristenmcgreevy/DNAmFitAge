# DNAmFitAge
This repository houses code and model information needed to calculate estimates for DNAm fitness biomarkers, DNAmFitAge, and FitAge Acceleration. 

This repository accompanies DNAmFitAge: Biological Age Indicator Incorporating Physical Fitness publication in Aging Albany. 

You can calculate DNAmGaitspeed, DNAmGripmax, DNAmFEV1, DNAmFitAge, and FitAgeAcceleration using blood DNA methylation. 


### DNAmFitnessandDNAmFitAgeSourceFunctions.R 
This includes 5 functions that will calculate the DNAm fitness estimates. 

### CodetoRun.R
Is the code users can run. This includes an example based on simulated data. The simulated data is kept in Sample_Simulation_Data.rds

### DNAmFitnessModelsandFitAge_Oct2022.rds
Houses the models used for each fitness biomarker and median values of CpGs. Median values of CpGs are used when your dataframe does not measure them. 

