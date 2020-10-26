rm(list=ls())

require(tidyverse)
require(ggplot2)
require(microsynth)
require(data.table)
require(here)

#####################
# READ IN DATA FILES:
#####################

allParks <- read.csv(paste(getwd(), "/Synthetic_controls_analysis/Synthetic_controls_allParks.csv", sep=""))
head(allParks)
nrow(allParks)

noMarineParks <- read.csv(paste(getwd(), "/Synthetic_controls_analysis/Synthetic_controls_noMarineParks.csv", sep=""))
head(noMarineParks)
nrow(noMarineParks)

noMarineParksOld <- read.csv(paste(getwd(),"/Synthetic_controls_analysis/Synthetic_controls_noMarineParksOld.csv", sep=""))
head(noMarineParksOld)
nrow(noMarineParksOld)


#############
# Microsynth:
#############

# define variables:

cov.var <- c("Desa_areaKm_In", 
	"Desa_areaKm_Out", 
	"GIS_M_AREA", 
	"GIS_AREA", 
	"STATUS_YR_2", 
	"avg_ParkSlope")

match.outR <- c(
	"Forest_Lost_In", 
	"Forest_Lost_Out", 
	"Total_Forest_In", 
	"Total_Forest_Out", 
	"Fires_In", 
	"Fires_Out", 
	"Population_In", 
	"Population_Out"
	)

match.resR <- c("Forest_Lost_In")


# run synthetic controls:

Microsynth <- microsynth(
	data=allParks, # for donor pool with all parks
	#data=noMarineParks, # for donor pool with no marine parks
	#data=noMarineParksOld, # for donor pool with no marine parks and no parks established after 2001
	idvar="DESA_CODE", timevar="Year", intvar="did", 
	start.pre=2001, end.pre=2007, end.post=2018,
	match.out.min=match.outR, match.covar.min=cov.var,
	result.var=match.resR, 
	omnibus.var=match.resR,
	test="lower",
	perm=500,
	jack=TRUE,
	cut.mse=2,
	use.survey=TRUE,
	use.backup=TRUE
	#,
	#result.file="Microsynth_results" # can name a results file to write to a csv
	)


# explore results
# note that the permuted p-values are based on slightly different placebo groups assembled for the permutation tests, for each run (i.e. the placebo groups for each of the 500 permutations are randomly assembled). Therefore the p-values and CI's obtained can be slightly different for each run.

summary(Microsynth)
Microsynth$Results


# make synthetic controls plots:

plot_microsynth(Microsynth, plot.var="Forest_Lost_In", 
	start.pre=2001, end.pre=2007, end.post=2018,
	legend.spot="topright",
	all=NULL,
	main.tc="Treatment and control",
	main.diff="Difference between groups",
	xlab.tc="Year",
	xlab.diff="Year",
	ylab.tc="Forest Lost",
	ylab.diff="Forest Lost"
	#,
	#file="Microsynth_plots.pdf" # can write plots to a pdf file
	)


