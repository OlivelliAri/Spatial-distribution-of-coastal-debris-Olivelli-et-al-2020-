## Settings and Libraries --------------------------------------------------
#set working directory
setwd("~/Documents/University/MScMarineSciences/CSIRO Internship/Ari_s project")

#load libraries
library(mgcv)
library(plyr)
library(dplyr) 
library(survival)
library(ggplot2)
library(MuMIn)
library(maps)
library(vcd)
library(stargazer)
library(zoo)
library(flexsurv)

## Preliminary Data Analysis -----------------------------------------------

#read in data
#Data <- read.csv('National_debris_distances.csv')
Data <- read.csv('National_debris_distances.csv',
                 colClasses = c(rep("factor",2),
                                "character",
                                rep("numeric",3),
                                rep("character",2),
                                rep("numeric",4),
                                rep("factor",6),
                                rep("numeric",23)))

#Anthropogenic forcing data
ExtraData <- read.csv('Population_data_180119.csv',
                      colClasses = c(rep("factor",2),
                                     "character",
                                     rep("numeric",3),
                                     rep("character",2),
                                     rep("numeric",8),
                                     rep("factor",2),
                                     rep("numeric",2))) 

CountData <- read.csv('National_debris_matrix.csv', 
                      colClasses = c(rep("factor",4),
                                     "character",
                                     rep("numeric",3),
                                     rep("factor",7),
                                     "character",
                                     rep("numeric",2),
                                     rep("character",2),
                                     rep("numeric",2),
                                     rep("character",2),
                                     rep("numeric",6),
                                     rep("factor",6),
                                     rep("numeric", (269-37+1))))

SizeData <- read.csv('National_debris_size_20181204.csv',
                     colClasses = c(rep("factor",5),
                                    "character",
                                    "factor",
                                    rep("numeric",3),
                                    rep("factor",7),
                                    "character",
                                    rep("numeric",2),
                                    rep("character",2),
                                    rep("numeric",2),
                                    rep("character",2),
                                    rep("numeric",6),
                                    rep("factor",7),
                                    rep("numeric",2),
                                    "character",
                                    "numeric",
                                    "character",
                                    "factor"))


#make a unique id to merge the different excel spreadsheets
Data$UniqueID <- paste(Data$BCH_STATE,Data$BCH_SITENAME,Data$BCH_SURVEY_DATE,Data$TRL_TRANSNUM, sep = ".", collapse = NULL)
ExtraData$UniqueID <- paste(ExtraData$BCH_STATE, ExtraData$BCH_SITENA, ExtraData$BCH_SURVEY, ExtraData$TRL_TRANSN, sep=".", collapse = NULL)
CountData$UniqueID <- paste(CountData$BCH_STATE, CountData$BCH_SITENAME, CountData$BCH_SURVEY_DATE,CountData$TRL_TRANSNUM, sep = ".", collapse = NULL)
SizeData$UniqueID <- paste(SizeData$BCH_STATE, SizeData$BCH_SITENAME, SizeData$BCH_SURVEY_DATE, SizeData$TRL_TRANSNUM, sep = ".", collapse = NULL)
index1 <- match(Data$UniqueID,CountData$UniqueID)
index2 <- match(Data$UniqueID,SizeData$UniqueID)
index3 <- match(Data$UniqueID,ExtraData$UniqueID)
Data$BCH_ID <- CountData$BCH_ID[index1]
Data$TRL_ID <- CountData$TRL_ID[index1]
Data$TotalDebris <- CountData$TOTAL[index1]
Data$BCH_WEATHER <- CountData$BCH_WEATHER[index1]
Data$BCH_WINDSPEED <- CountData$BCH_WINDSPEED[index1]
Data$BCH_WINDCOMPASS <- CountData$BCH_WINDCOMPASS[index1]
Data$BCH_WINDSHORE <- CountData$BCH_WINDSHORE[index1]
Data$Pop5km <- ExtraData$POP_5KM[index3]
Data$Pop20km <- ExtraData$POP_20KM[index3]
Data$Pop50km <- ExtraData$POP_50KM[index3]
Data$Road_dist <- ExtraData$Rd_Dist_KM[index3]
Data$Watershd_id <- ExtraData$Wtshd_FID[index3]
Data$Watershd_name <- ExtraData$Wtrshed_Na[index3]
Data$Watershd_size <- ExtraData$Wtrshd_Km2[index3]
Data$Watershd_pop <- ExtraData$Wtrshd_pop[index3]


#shift the data to long format from wide
DataLong <- data.frame((levels(Data$BCH_ID)[rep(Data$BCH_ID, each = 10)]),
                       levels(Data$TRL_ID)[rep(Data$TRL_ID, each = 10)],
                       levels(Data$BCH_STATE)[rep(Data$BCH_STATE,each = 10)],
                       levels(Data$BCH_SITENAME)[rep(Data$BCH_SITENAME, each = 10)],
                       levels(Data$BCH_WEATHER)[rep(Data$BCH_WEATHER, each = 10)],
                       levels(Data$BCH_WINDSPEED)[rep(Data$BCH_WINDSPEED, each = 10)],
                       levels(Data$BCH_WINDCOMPASS)[rep(Data$BCH_WINDCOMPASS, each = 10)],
                       levels(Data$BCH_WINDSHORE)[rep(Data$BCH_WINDSHORE, each = 10)],
                       rep(Data$TRL_STARTLAT, each = 10),
                       rep(Data$TRL_STARTLONG, each = 10),
                       rep(Data$TRL_ENDLAT, each = 10),
                       rep(Data$TRL_ENDLONG, each = 10),
                       levels(Data$TRL_GRADIENT)[rep(Data$TRL_GRADIENT, each = 10)],
                       levels(Data$TRL_SUBSTRATE)[rep(Data$TRL_SUBSTRATE, each = 10)],
                       levels(Data$TRL_SUBCOLOUR)[rep(Data$TRL_SUBCOLOUR, each = 10)],
                       levels(Data$TRL_BACKSHORE)[rep(Data$TRL_BACKSHORE, each = 10)],
                       levels(Data$TRL_SHAPE)[rep(Data$TRL_SHAPE, each = 10)],
                       levels(Data$TRL_ASPECT)[rep(Data$TRL_ASPECT, each = 10)],
                       rep(Data$TRL_LENGTH, each = 10),
                       rep(Data$TRL_DEBRISLINE_M, each = 10),
                       rep(Data$SECT_LENGTH, each = 10),
                       rep(Data$TotalDebris, each = 10),
                       rep(Data$Pop5km, each = 10),
                       rep(Data$Pop20km, each = 10),
                       rep(Data$Pop50km, each = 10), 
                       rep(Data$Road_dist, each = 10),
                       levels(Data$Watershd_id)[rep(Data$Watershd_id, each = 10)],
                       levels(Data$Watershd_name)[rep(Data$Watershd_name, each = 10)],
                       rep(Data$Watershd_size, each = 10),
                       rep(Data$Watershd_pop, each = 10),
                       rep(Data$UniqueID, each = 10),
                       rep(as.factor(c(1,2,3,4,5,6,7,8,9,10))),
                       Data[,c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")][cbind(rep(1:dim(Data)[1],each = 10),rep(1:10,times = dim(Data)[1]))],
                       Data[,c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10")][cbind(rep(1:dim(Data)[1],each = 10),rep(1:10,times = dim(Data)[1]))])


#set names
names(DataLong) <- c("BeachID", "TransID", "State",
                     "BeachName","Weather","WindSpeed",
                     "WindCompass","WindShore","StartLat",
                     "StartLong","EndLat", "EndLong",
                     "Gradient","Substrate","SubColour",
                     "Backshore","Shape","Aspect",
                     "TransectLength", "DebrisLine", "SectionLength",
                     "TotalDebris", "Population_5km", 
                     "Population_20km", "Population_50km", "RoadDistance", 
                     "WatershedID", "WatershedName", "WatershedSize", 
                     "WatershedPop", "ID", "SectionNumber","Distance","Status")

#The size category excel spreadsheet already had 6350 entries (while the others had 635)
SizeCategory <- vector(length = length(index2)*10)
for (i in 1:length(index2)){
  SizeCategory[(i+9*(i-1)):(i*10)] <- SizeData$SIZ_DEBRSIZE[index2[i]:(index2[i]+9)]}
SizeCategory <- SizeCategory - 1 #so that the lowest size category is zero (when debris was not found)
SizeCategory <- factor(SizeCategory, ordered = T)
DataLong$SizeCategory <- SizeCategory
#Make sure that when SizeCatgory = 0, Status = 0
DataLong$Status[DataLong$SizeCategory == "0"] <- 0
DataLong$Status[DataLong$SizeCategory != "0"] <- 1
table(DataLong$Status,DataLong$SizeCategory)

#Order categorical variables (when possible) and make numeric variables
DataLong$WindSpeed <- factor(DataLong$WindSpeed, ordered = T)
DataLong$Gradient <- factor(DataLong$Gradient, ordered = T)
DataLong$NumericGradient <- as.numeric(DataLong$Gradient)
gradient <- matrix(c(0,1,1,2,2,4,4,8,8,Inf), byrow = T, ncol = 2)
DataLong$LowerNumGradient <- gradient[DataLong$NumericGradient,1]
DataLong$WeatherOr <- factor(DataLong$Weather, ordered = T, levels = c("Clear", "Overcast","Drizzle","Rain/Storm"))
DataLong$SectionNumber <- factor(DataLong$SectionNumber, ordered = T)
DataLong$NumericSectionNumber <- as.numeric(DataLong$SectionNumber)
wind <- c(0,5,17.5,37,57,75)
DataLong$NumericWindspeed <- wind[DataLong$WindSpeed]
DataLong$NumericSizeCategory <- as.numeric(DataLong$SizeCategory)
size <- matrix(c(0,1.2,1.2,2.4,2.4,4.9,4.9,9.6,9.6,19.2,19.2,36.3,36.3,Inf),byrow =T, ncol =2)
DataLong$LowerBSize <- size[DataLong$NumericSizeCategory,1]

#Scale anthropogenic forcing covariates (otherwise the model won't converge later)
DataLong$Population_20km <- scale(DataLong$Population_20km)
DataLong$Population_50km <- scale(DataLong$Population_50km)
DataLong$Population_5km <- scale(DataLong$Population_5km)
DataLong$RoadDistance <- scale(DataLong$RoadDistance)

#calculations
DataLong$RelDistance <- DataLong$Distance/DataLong$SectionLength
DataLong$RelDistance[DataLong$RelDistance==0] <- 0.0001
DataLong$RelDebrisLine <- DataLong$DebrisLine/DataLong$TransectLength
DataLong$RelDebrisLine[is.na(DataLong$RelDebrisLine)] <- 1


