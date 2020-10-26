# Download required R libraries
require(tidyverse)
require(lme4)
require(lmerTest)
require(glmmTMB)

# download data
forest_change_data <- read.csv(paste(getwd(),"/Forest_loss_dose_response_analysis/dose_response_data.csv", sep=""))

# Variable `Year` is the last year wherein 5-year forest change was estimated. Year = 2006 is the forest change rates from 2002-2006 in the 5 years before the intervention began in 2007 (considered a "burn-in" year); Year = 2012 is the forest change rates from 2008-2012 in the first 5 years of the intervention; Year = 2017 is forest change ranges from 2013-2017 in the last 5 years of the intervention

## Estimate forest change rates in low-, medium-, and highly-engaged villages (Dusun) over time
low_engaged <- lmer(PercLossAnnual_inside ~ scale(Year) + (1|Dusun), data = forest_change_data[forest_change_data$EngagementBinned=="1low",])
summary(low_engaged)

medium_engaged <- lmer(PercLossAnnual_inside ~ scale(Year) + (1|Dusun), data = forest_change_data[forest_change_data$EngagementBinned=="2medium",])
summary(medium_engaged)

high_engaged <- lmer(PercLossAnnual_inside ~ scale(Year) + (1|Dusun), data = forest_change_data[forest_change_data$EngagementBinned=="3high",])
summary(high_engaged)

## Estimate probability that any forested pixel inside medium- and high-engaged village access areas was lost compared to low-engaged access areas, from the 5-year period before the intervention to the most recent 5-year period, while controlling for average park slope, elevation, distance to nearest road, distance to park boundary, forest loss outside of park, and village population
mod <- glmmTMB(cbind(LossTotal5yr_inside, ForestTotal5yr_inside-LossTotal5yr_inside) ~ 
                log(Population_desa)
              + log(LossTotal5yr_outside+1) 
              + log(Elevation) + log(Slope) + log(DistanceRiver) + log(DistanceRoad) + log(DistanceEdge) 
              + factor(EngagementBinned)*scale(Year)
              + (1|Desa/Dusun),
              family="binomial", 
              data=forest_change_data)
summary(mod)

