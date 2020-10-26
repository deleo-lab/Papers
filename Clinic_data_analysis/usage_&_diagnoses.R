# Download required R libraries
require(here)
require(tidyverse)
require(MASS) 
require(lme4) 
require(broom) 

# optional, suppress warnings
warn_on <- getOption("warn")
options(warn = -1)

########## Statistical analyses on clinic usage ########## 
############ ########## ########## ########## ############ 

#### Statistics on proportion of district populations accessing clinic ####
# Test for effect of MOU on the proportion of a district's population that accessed the clinic at least once, controlling for each district's distance to the clinic ##
# Load data
pop.usage.data <- read_csv(paste(getwd(),"/Clinic_data_analysis/ClinicPopData.csv", sep=""))
pop.usage.data$MOUYN <- as.factor(pop.usage.data$MOUYN)

# Apply model
distance.usage <- glmer(UniquePatients ~ scale(GIS_time) + MOUYN + offset(log(Desa_pop_2018)) + (1|Desa), family="poisson", data=pop.usage.data[pop.usage.data$Desa!="Far",])
summary(distance.usage)

#### Statistics on frequency of unique patient visits ####
# Test for effect of MOU on the frequency that a unique patient visited the clinic, controlling for patient's distance to the clinic #

# Download clinic database records from 2008 to 2018 
clinic_data <- read_csv(paste(getwd(),"/Clinic_data_analysis/ClinicPatientData.csv", sep=""))

# Apply model
patient.freq <- glmer.nb(NumMonths ~ scale(GIS_time) + factor(MOUYN) + (1|Desa), data=patient.freq.data)
summary(patient.freq) 


###### Statistical analyses on clinic disease diagnosis trends ###### 
############ ############# ############# ############## #############

# Download data for diagnoses assigned during general doctor visits; ICD10 codes (available in supplementary materials) have been previously categorized and here are summarized as village-level frequencies for each diagnosis category 

diagnosis_data <- read_csv(paste(getwd(), "/Clinic_data_analysis/ClinicDiagnosisData.csv", sep=""))

# set global glmerControl options
opts <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6), check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol))

#### Test for change in diagnosis trends over time, for MOU-signing patients and non-MOU-signing patients, for ICD10 diagnosis categories of interest #### 

# Apply model to each disease category for MOU-signing districts

# Statistical analysis excludes liver disease due to lack of observations
diagnosis_data2 <- diagnosis_data[diagnosis_data$DiseasestoPlot!="Liver Disease",]

change_disease_MOU <- expand.grid(DiseasestoPlot=unique(diagnosis_data2$DiseasestoPlot))

for(i in 1:length(change_disease_MOU$DiseasestoPlot)) {
  test<-glmer(cbind(NumDiagnoses, (TotalYearlyPatients-NumDiagnoses)) ~ time + scale(GIS_time) + (1|Desa),  family="binomial", opts, data=diagnosis_data2[diagnosis_data2$MOUYN == 1 & diagnosis_data2$DiseasestoPlot==change_disease_MOU$DiseasestoPlot[i],])
  change_disease_MOU$MOUYN <- 1
  change_disease_MOU$estimate_late[i]<-tidy(test)[2,2]
  change_disease_MOU$OR[i]<-exp(tidy(test)[2,2])
  change_disease_MOU$OR_min[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,6]) 
  change_disease_MOU$OR_max[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,7]) 
  change_disease_MOU$z_statistic[i]<-tidy(test)[2,4]
  change_disease_MOU$pvalue_late[i]<-tidy(test)[2,5]
}
View(change_disease_MOU)

change_disease_noMOU <- expand.grid(DiseasestoPlot=unique(diagnosis_data2$DiseasestoPlot))

for(i in 1:length(change_disease_noMOU$DiseasestoPlot)) {
  test<-glmer(cbind(NumDiagnoses, (TotalYearlyPatients-NumDiagnoses)) ~ time + scale(GIS_time) + (1|Desa),  family="binomial", opts, data=diagnosis_data2[diagnosis_data2$MOUYN == 0 & diagnosis_data2$DiseasestoPlot==change_disease_noMOU$DiseasestoPlot[i],])
  change_disease_noMOU$MOUYN <- 0
  change_disease_noMOU$estimate_late[i]<-tidy(test)[2,2]
  change_disease_noMOU$OR[i]<-exp(tidy(test)[2,2])
  change_disease_noMOU$OR_min[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,6]) 
  change_disease_noMOU$OR_max[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,7]) 
  change_disease_noMOU$z_statistic[i]<-tidy(test)[2,4]
  change_disease_noMOU$pvalue_late[i]<-tidy(test)[2,5]
}
View(change_disease_noMOU)

# Bind data together
ModelOutput <- unnest(rbind(change_disease_MOU, change_disease_noMOU))
# Add disease category to table for plotting
ModelOutput$Disease_cat <- diagnosis_data2$Disease_cat[match(ModelOutput$DiseasestoPlot, diagnosis_data2$DiseasestoPlot)]

# Make Odds ratio plot
dodger = position_dodge(width = 0.9)
ggplot(ModelOutput, aes(y = OR, x = reorder(DiseasestoPlot, -OR), color = DiseasestoPlot, shape = factor(MOUYN))) +
  geom_pointrange(stat="identity", aes(ymin = OR_min, ymax = OR_max),position = dodger, size = 1) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.01, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10, 20)) +
  labs(y = "Odds ratio", x = "Disease Category") +
  coord_flip(ylim = c(0.001, 20)) +
  facet_grid(Disease_cat~., drop = TRUE, scale="free_y") +
  theme_classic()

# Repeat, not controlling for distance, and including liver disease, to show that MOU effects are consistent whether or not controlling for distance
change_disease_MOU <- expand.grid(DiseasestoPlot=unique(diagnosis_data$DiseasestoPlot))

for(i in 1:length(change_disease_MOU$DiseasestoPlot)) {
  test<-glmer(cbind(NumDiagnoses, (TotalYearlyPatients-NumDiagnoses)) ~ time + (1|Desa),  family="binomial", opts, data=diagnosis_data[diagnosis_data$MOUYN == 1 & diagnosis_data$DiseasestoPlot==change_disease_MOU$DiseasestoPlot[i],])
  change_disease_MOU$MOUYN <- 1
  change_disease_MOU$estimate_late[i]<-tidy(test)[2,2]
  change_disease_MOU$OR[i]<-exp(tidy(test)[2,2])
  change_disease_MOU$OR_min[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,6]) 
  change_disease_MOU$OR_max[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,7]) 
  change_disease_MOU$z_statistic[i]<-tidy(test)[2,4]
  change_disease_MOU$pvalue_late[i]<-tidy(test)[2,5]
}
View(change_disease_MOU)

change_disease_noMOU <- expand.grid(DiseasestoPlot=unique(diagnosis_data$DiseasestoPlot))

for(i in 1:length(change_disease_noMOU$DiseasestoPlot)) {
  test<-glmer(cbind(NumDiagnoses, (TotalYearlyPatients-NumDiagnoses)) ~ time + (1|Desa),  family="binomial", opts, data=diagnosis_data[diagnosis_data$MOUYN == 0 & diagnosis_data$DiseasestoPlot==change_disease_noMOU$DiseasestoPlot[i],])
  change_disease_noMOU$MOUYN <- 0
  change_disease_noMOU$estimate_late[i]<-tidy(test)[2,2]
  change_disease_noMOU$OR[i]<-exp(tidy(test)[2,2])
  change_disease_noMOU$OR_min[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,6]) 
  change_disease_noMOU$OR_max[i]<-exp(tidy(test, conf.int=TRUE, effects="fixed")[2,7]) 
  change_disease_noMOU$z_statistic[i]<-tidy(test)[2,4]
  change_disease_noMOU$pvalue_late[i]<-tidy(test)[2,5]
}
View(change_disease_noMOU)

# Bind data together
ModelOutput <- unnest(rbind(change_disease_MOU, change_disease_noMOU))
# Add disease category to table for plotting
ModelOutput$Disease_cat <- diagnosis_data$Disease_cat[match(ModelOutput$DiseasestoPlot, diagnosis_data$DiseasestoPlot)]

dodger = position_dodge(width = 0.9)
ggplot(ModelOutput, aes(y = OR, x = reorder(DiseasestoPlot, -OR), color = DiseasestoPlot, shape = factor(MOUYN))) +
  geom_pointrange(stat="identity", aes(ymin = OR_min, ymax = OR_max),position = dodger, size = 1) +
  geom_hline(yintercept = 1.0, linetype = "dotted", size = 0.5) +
  scale_y_log10(breaks = c(0.01, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10, 20)) +
  labs(y = "Odds ratio", x = "Disease Category") +
  coord_flip(ylim = c(0.001, 20)) +
  facet_grid(Disease_cat~., drop = TRUE, scale="free_y") +
  theme_classic()


#### Test for effect of participation in the intervention (signing of MOU) on diagnosis trends for each disease category #### 
# Here, we use an interaction term between MOU-status (MOUYN) and time (Beginning of intervention vs. End of intervention) to estimate the intervention impact on contact incidence for each disease category, while controlling for desa-level distance to the clinic

interactions <- expand.grid(DiseasestoPlot=unique(diagnosis_data$DiseasestoPlot))
interactions$Disease_cat <- diagnosis_data$Disease_cat[match(interactions$DiseasestoPlot, diagnosis_data$DiseasestoPlot)]
  
model.list <- list()

for(i in 1:length(interactions$DiseasestoPlot)) {
  test <- glmer(cbind(NumDiagnoses, TotalYearlyPatients-NumDiagnoses) ~ time*factor(MOUYN) + scale(GIS_time) + (1|Desa), data=diagnosis_data[diagnosis_data$DiseasestoPlot==interactions$DiseasestoPlot[i],], family="binomial", opts)
  model.list[i] <- test
  interactions$estimate[i]<-unlist(tidy(test)[5,2])
  interactions$OR[i]<-unlist(exp(tidy(test)[5,2]))
  interactions$OR_min[i]<-unlist(exp(tidy(test, conf.int=TRUE, effects="fixed")[5,6]))
  interactions$OR_max[i]<-unlist(exp(tidy(test, conf.int=TRUE, effects="fixed")[5,7]))
  interactions$z_statistic[i]<-unlist(exp(tidy(test)[5,3]))
  interactions$pvalue[i]<-unlist(tidy(test)[5,5])
}

View(interactions)

# turn warnings back on
options(warn = warn_on)