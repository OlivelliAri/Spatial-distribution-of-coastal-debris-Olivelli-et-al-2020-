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
source("Dredgeterms.R")
source("Effect Size plots_surv.R")

#load data
load("~/Documents/University/MScMarineSciences/CSIRO Internship/Ari_s project/final_environment.RData")

## Survival analysis for the whole dataset (model excludes SizeCategory). Total model in the paper --------------
#figure out what are the factors playing a role when size is excluded
options(na.action = "na.fail")
M.global_noSize <- survreg(Surv(RelDistance, Status) ~ State + Weather +  Substrate + SubColour + Backshore + LowerNumGradient + Shape + Aspect + NumericSectionNumber + TotalDebris,dist="exponential", data = DataLong) 
M.dredge_noSize <- dredge(M.global_noSize)
plot(M.dredge_noSize, labAsExpr = TRUE)
M.global_noSize.1 <- survreg(Surv(RelDistance, Status) ~  Aspect + Backshore + NumericSectionNumber + State + SubColour + Substrate + TotalDebris + Weather, dist = "exponential", data = DataLong)
summary(M.global_noSize.1)


## Survival analysis for those sections where debris was found (Status = 1). Positives model in the paper -------------
DataLongPositives <- DataLong[DataLong$Status == 1,] #New dataset with positives only

M.global_inclSize <- survreg(Surv(RelDistance, Status) ~ State + Weather + Substrate + SubColour + Backshore + LowerNumGradient  + Shape + Aspect  + TotalDebris + LowerBSize*NumericSectionNumber,dist="exponential", data = DataLongPositives) 
M.dredge_inclSize <- dredge(M.global_inclSize)
plot(M.dredge_inclSize, labAsExpr = TRUE)
M.global_inclSize.1 <- survreg(Surv(RelDistance, Status) ~ LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber, dist = "exponential", data = DataLong)
summary(M.global_inclSize.1)

## Stokes and Wind Components - Functions ----------------------------------------------

#rm(list = ls())

FindDistance <- function(LatsStart, LongsStart,LatsEnd,LongsEnd){
  #this function calculates the great circle distance between two points or vectors of
  #points.  Note that it assumes the earth is spherical, so has an error of ~0.3%
  #formulae are from http://www.movable-type.co.uk/scripts/latlong.html
  #
  #note that things need to be converted to numeric format if they are not already in it.
  #the function can handle a single starting point and a vector of ending points.  In this
  #case it returns a vector of distances
  
  EarthRadius <- 6371 #in km, so answer is in km
  
  #convert lat, long and bearings from degrees to radians using rad = pi*deg/180
  LatsStart <- (22/7)*LatsStart/180
  LongsStart <- (22/7)*LongsStart/180
  LatsEnd <- (22/7)*LatsEnd/180
  LongsEnd <- (22/7)*LongsEnd/180
  
  DeltaLat <- LatsStart - LatsEnd
  DeltaLong <- LongsStart - LongsEnd
  a <- (sin(DeltaLat/2))^2 + (cos(LatsStart)*cos(LatsEnd)*((sin(DeltaLong/2))^2))
  c <- 2*atan2(a^0.5,(1-a)^0.5)
  Distance <- EarthRadius * c
  return(Distance)
}


FindDepartBearing <- function(LatsStart, LongsStart,LatsEnd,LongsEnd,Switch=T){
  #This formula is for the initial bearing (sometimes referred to as forward azimuth) which if
  #followed in a straight line along a great-circle arc will take you from the start point to the 
  #end point.  Taken from http://www.movable-type.co.uk/scripts/latlong.html.  Note that the bearing
  #will change as you are following a great circle, ie taking the shortest route.  Note that it 
  #returns the bearing in 360 degrees, remove the last step to get it from -180 to + 180 which 
  #is deviation off true north. 
  #
  #I added a variable called switch, if its true then the function returns it on 360, if false -180 to 180
  #
  #Note Trig functions take arguments in radians, so latitude, longitude, and bearings in degrees (either decimal or degrees/minutes/seconds) 
  #need to be converted to radians, rad = ??.deg/180. When converting radians back to degrees (deg = 180.rad/??), West is negative if using 
  #signed decimal degrees. For bearings, values in the range -?? to +?? [-180? to +180?] need to be converted to 0 to +2?? [0?-360?]; this 
  #can be done by (brng+2.??)%2.?? [or brng+360)%360] where % is the modulo operator.
  #
  
  
  #convert lat, long and bearings from degrees to radians using rad = pi*deg/180
  LatsStart <- (22/7)*LatsStart/180
  LongsStart <- (22/7)*LongsStart/180
  LatsEnd <- (22/7)*LatsEnd/180
  LongsEnd <- (22/7)*LongsEnd/180
  
  DeltaLat <- LatsStart - LatsEnd
  DeltaLong <- LongsStart - LongsEnd
  
  Bearings <- atan2(sin(DeltaLong) * cos(LatsEnd), cos(LatsStart)*sin(LatsEnd) - sin(LatsStart)*cos(LatsEnd)*cos(DeltaLong))
  
  #Since atan2 returns values in the range -?? ... +?? (that is, -180? ... +180?), to normalise the result to a compass bearing (in the range 0? ... 360?, 
  #with -ve values transformed into the range 180? ... 360?), convert to degrees and then use (??+360) % 360, where % is modulo.
  Bearings <- Bearings*180/(22/7)
  if(Switch) Bearings <- (Bearings + 360) %% 360
  
  return(Bearings)
}

FindDestination <- function(LatsStart, LongsStart, Bearing, Distance){
  #This function finds the end lat and long, given the starting points and a distance and 
  #bearing.  Distance/EarthRadius is the angular distance (in radians) where Distance 
  #is the distance travelled and EarthRadius is obvious.  Taken from http://www.movable-type.co.uk/scripts/latlong.html
  #All input angles are on 360 degrees.
  
  EarthRadius <- 6371
  
  #convert lat, long and bearings from degrees to radians using rad = pi*deg/180
  LatsStart <- (22/7)*LatsStart/180
  LongsStart <- (22/7)*LongsStart/180
  Bearing <- (22/7)*Bearing/180
  
  LatsEnd <- asin(sin(LatsStart)*cos(Distance/EarthRadius) + cos(LatsStart)*sin(Distance/EarthRadius)*cos(Bearing))
  LongsEnd <- LongsStart + atan2(sin(Bearing)*sin(Distance/EarthRadius)*cos(LatsStart),cos(Distance/EarthRadius) - sin(LatsStart)*sin(LatsEnd))
  
  Destinations <- cbind(LatsEnd, LongsEnd)
  
  #convert back from radians
  Destinations <- Destinations*180/(22/7)
  
  return(Destinations)
}

FitBoxCox <- function(l,obs){
  #this function fits the lambda in the box cox transformation based on the 
  #value that maximises the shapiro test 
  BCTrans <- (obs[!is.na(obs)] ^ l - 1)/l
  
  if((!(sum(BCTrans == 0) == length(BCTrans))) &  (length(unique(BCTrans)) != 1) ) {
    shapirovalues <- shapiro.test(BCTrans)
  } else {
    
    shapirovalues <- c(NA, NA)
  }
  
  #output <- shapirovalues[1]	#test statistic
  #output <- shapirovalues[1:2]	#both
  output <- shapirovalues[2]		#p value, used this to work with optim
  
  return(output)
}

#this function does the box cox transform for a given pwer value
BoxCoxTrans <- function(l,obs){
  #added an exemption for l being na, in which case there is no transform
  
  output <- obs
  if (!(is.na(l))){output <- (obs^l - 1)/l}
  
  return(output)
}

FindBearing <- function(U,V,Switch = F){
  #this function takes in a U and V, on the unit circle.  It returns the bearing of the vector from the origin
  #This was designed to be used with North as the origin, so it will give degrees of rotation right of north as positive and 
  #degrees of rotation left of north as negative.  All resulting angles should be in the range -180,180
  #This ignores the curvature of the earth, so its fine, but should not be used with great circle distances. 
  #
  #U is longitude, V is latitude.  Positive U is wind blowing to the east.
  
  #
  #I modified it so
  #that if you needed wind going counterclockwise 360 degrees, you could just change Switch to get it.  This is useful for using
  #stars plotting, which starts at horizontal right and goes counterclockwise.
  
  #to get atan to work, giving positive rotation right, and negative left, rotate 90 degrees to left (ie v goes to u, and keep sign)
  #then multiply by -1 to get correct rotational direction
  Bearing <- atan2(U,V)*180/pi
  
  #this is an earlier version dealing with the annoying thing that positive v is west for oceanographers, east for everyone else
  #if(U >= 0) Bearing <- -1 *asin(V)*(180/pi)
  #if(U < 0 & V >= 0) Bearing <- (asin(V)*(180/pi)) - 180  
  #if(U < 0 & V < 0) Bearing <- 180 + (asin(-1)*(180/pi))
  
  if(Switch == T){
    if(Bearing < 0) Bearing <- 360 + Bearing
  }
  
  return(Bearing)
  
}


## Stokes calculation from physical forcing data ------------------------- 

#first of all, make sure that the dataset have the same structure and the transects are in the same order! 
#To do this, check the lambda and wind indeces with respect to Data (they need to be in order from 1 to 635)
lambda <- read.csv('Wavelength.csv')
lambda_mat <- as.matrix(lambda[,7:8766])
lambda_matT <- t(lambda_mat) #rollapply works on a vertical structure hence the transposed matrix
lambda_mean <- rollapply(lambda_matT, 24, mean, by = 24, align = "left") #calucluates the daily mean (as we have hourly data)
wavenumber <- 1/lambda_mean
#lambda$UniqueID <- paste(lambda$BCH_STATE,lambda$BCH_SITENAME,lambda$BCH_SURVEY_DATE,lambda$TRL_TRANSNUM, sep = ".", collapse = NULL)
#index4 <- match(Data$UniqueID, lambda$UniqueID)
#k <- cbind(lambda$UniqueID, t(wavenumber))
k <- t(wavenumber) #transpose again to go back to the original shape 
k <- k[rep(1:nrow(k), each = 10),]
h <- 20 #water depth
z<- 19.8 #height above seabed in [m]

Hs <- read.csv('Waveheight.csv')
Hs_mat <- as.matrix(Hs[,7:8766])
Hs_matT <- t(Hs_mat)
Hs_mean <- rollapply(Hs_matT, 24, mean, by = 24, align = "left")
#Hs$UniqueID <- paste(Hs$BCH_STATE,Hs$BCH_SITENAME,Hs$BCH_SURVEY_DATE,Hs$TRL_TRANSNUM, sep = ".", collapse = NULL)
#index5 <- match(Data$UniqueID, Hs$UniqueID)
#Hs_mean <- cbind(Hs$UniqueID,t(Hs_mean))
Hs_mean <- t(Hs_mean)
Hs_mean <- Hs_mean[rep(1:nrow(Hs_mean), each = 10),]

#Formulas from Bever et al. (2011)
Amp <- Hs_mean/2
E <- (9.81*Amp^2)/2
c <- sqrt(9.81*tanh(k*h)/k)
Stokes_u = (2*k*E*cosh(2*k*z))/(c*sinh(2*k*h))

Wind_u <- read.csv('Wind_u.csv')
#Wind_u$UniqueID <- paste(Wind_u$BCH_STATE,Wind_u$BCH_SITENAME,Wind_u$BCH_SURVEY_DATE,Wind_u$TRL_TRANSNUM, sep = ".", collapse = NULL)
#index6 <- match(Data$UniqueID,Wind_u$UniqueID)
Wind_u_mat <- as.matrix(Wind_u[,7:8766])
Wind_u_matT <- t(Wind_u_mat)
Wind_u_mean <- rollapply(Wind_u_matT, 24, mean, by = 24, align = "left")
Wind_u_mean <- t(Wind_u_mean)
Wind_u_mean <- Wind_u_mean[rep(1:nrow(Wind_u_mean), each = 10),]

Wind_v <- read.csv('Wind_v.csv')
#Wind_v$UniqueID <- paste(Wind_v$BCH_STATE,Wind_v$BCH_SITENAME,Wind_v$BCH_SURVEY_DATE,Wind_v$TRL_TRANSNUM, sep = ".", collapse = NULL)
#index7 <- match(Data$UniqueID,Wind_v$UniqueID)
Wind_v_mat <- as.matrix(Wind_v[,7:8766])
Wind_v_matT <- t(Wind_v_mat)
Wind_v_mean <- rollapply(Wind_v_matT, 24, mean, by = 24, align = "left")
Wind_v_mean <- t(Wind_v_mean)
Wind_v_mean <- Wind_v_mean[rep(1:nrow(Wind_v_mean), each = 10),]

Wind_velocity <- sqrt(Wind_u_mean^2 + Wind_v_mean^2)
## Stokes and Wind forcing - Calculations -------------------------------

#Working on calculating the relative deviation between wind and coast direction
#Use the bearing function to calculate the deviation from north, note that deviation to right/clockwise is positive, left negative
WindAngle <- FindBearing(Wind_u_mean,Wind_v_mean)
#sometimes U and V are zero, so we get an NA for WindAngle.  These should be replaced with 0, to allow caculations.  The onshore
#component will end up 0, as there is no velocity, so it doesnt really matter.
WindAngle[is.na(WindAngle)] <- 0
#for Stokes drift we only have the shoreward component, so we can't do this

#we want the offshore direction of the coast.  The simplest is to just take the transect gps coordinates
OffshoreBearing <- FindDepartBearing(DataLong$EndLat, DataLong$EndLong, DataLong$StartLat, DataLong$StartLong,Switch = F)
#Get the difference in the angles, which will allow us to calculate the onshore
#componet.  For differences (relative wind angles) greater than 180, we should use the complement, as the angle from north overestimates 
#the difference eg -360 is equivalent to a difference of 0.  So we want the minimum enclosed angle.  We should be able to do this by using
#abs and then just getting the ones greater than 180 and subtracting from 360
RelativeWindAngle <- abs(WindAngle - OffshoreBearing)
### fix this HalfRelativeWindAngle <- matrix(nrow =dim(RelativeWindAngle)[1], ncol=dim(RelativeWindAngle)[2])
HalfRelativeWindAngle <- RelativeWindAngle*NA
HalfRelativeWindAngle[RelativeWindAngle <= 180] <- RelativeWindAngle[RelativeWindAngle <= 180]
HalfRelativeWindAngle[RelativeWindAngle > 180] <- 360 - RelativeWindAngle[RelativeWindAngle > 180]
#again this is not the case for Stokes

#Now add the wind and Stokes data to DataLong
DataLong$MeanHalfWindAngle <- apply(HalfRelativeWindAngle,1,mean)
#not available for Stokes

#might want to consider the velocity, compute the onshore component of the wind
DataLong$MeanOnshoreWind <- apply(Wind_velocity * cos((HalfRelativeWindAngle) * pi/180),1,mean)
DataLong$MeanOnshoreStokes <- apply(Stokes_u,1,mean)

#Now, just use a month, and a week, and a day before the sampling
DataLong$MeanOnshoreWind1Mo <- apply(Wind_velocity[,1:30] * cos((HalfRelativeWindAngle[,1:30]) * pi/180),1,mean)
DataLong$MeanOnshoreWind1Wk <- apply(Wind_velocity[,1:7] * cos((HalfRelativeWindAngle[,1:7]) * pi/180),1,mean)
DataLong$MeanOnshoreWind1Dy <- apply(matrix(Wind_velocity[,1] * cos((HalfRelativeWindAngle[,1]) * pi/180),ncol = 1),1,mean)
DataLong$MeanOnshoreWind10Dy <- apply(Wind_velocity[,1:10] * cos((HalfRelativeWindAngle[,1:10]) * pi/180),1,mean)
DataLong$MeanOnshoreWindHyr <- apply(Wind_velocity[,1:183] * cos((HalfRelativeWindAngle[,1:183]) * pi/180),1,mean)
DataLong$MeanOnshoreStokes1Mo <- apply(Stokes_u[,1:30],1,mean)
DataLong$MeanOnshoreStokes1Wk <- apply(Stokes_u[,1:7],1,mean)
DataLong$MeanOnshoreStokes1Dy <- Stokes_u[,1]
DataLong$MeanOnshoreStokes10Dy <- apply(Stokes_u[,1:10],1,mean)
DataLong$MeanOnshoreStokesHyr <- apply(Stokes_u[,1:183],1,mean)

## Forced survival model for the whole dataset (excluding SizeCategory) ------------
#Physical forcing
numvar1ph <- c("MeanOnshoreWind1Mo","MeanOnshoreWind1Wk","MeanOnshoreWind1Dy","MeanOnshoreWind10Dy", "MeanOnshoreWindHyr",
             "MeanOnshoreStokes1Mo","MeanOnshoreStokes1Wk","MeanOnshoreStokes1Dy", "MeanOnshoreStokes10Dy", "MeanOnshoreStokesHyr",
             "TotalDebris", "NumericSectionNumber")
excl1ph <- Dredgeterms(numvar1ph,DataLong) #Creates a subset of collinear variables that are not called together when dredging the global model
M.global_noSize.Fph <- survreg(Surv(RelDistance, Status) ~ Aspect + Backshore + NumericSectionNumber + State + SubColour + Substrate + TotalDebris + Weather + MeanOnshoreWind1Mo
                               + MeanOnshoreWind1Wk + MeanOnshoreWind1Dy + MeanOnshoreWind10Dy + MeanOnshoreWindHyr + MeanOnshoreStokes1Mo
                               + MeanOnshoreStokes1Wk + MeanOnshoreStokes1Dy + MeanOnshoreStokes10Dy + MeanOnshoreStokesHyr, dist="exponential", data = DataLong)
M.dredge_noSize.Fph  <- dredge(M.global_noSize.Fph, subset = c(excl1ph))
#Best observation + physical forcing model
M.global_noSize.Fph.1 <- survreg(Surv(RelDistance, Status) ~ Aspect + Backshore + MeanOnshoreStokes1Mo + MeanOnshoreWind1Wk 
                               + NumericSectionNumber + State + SubColour + Substrate + TotalDebris + Weather, dist = "exponential", data = DataLongFull)
summary(M.global_noSize.Fph.1)

#Anthropogenic forcing
numvar1so <- c("Population_5km","Population_20km","Population_50km","RoadDistance",
                 "TotalDebris", "NumericSectionNumber")
excl1so <- Dredgeterms(numvar1so,DataLong)
M.global_noSize.Fso <- survreg(Surv(RelDistance, Status) ~ Aspect + Backshore + NumericSectionNumber + State + SubColour + Substrate + TotalDebris + Weather 
                               + Population_5km + Population_20km + Population_50km + Population_50km*RoadDistance + RoadDistance 
                              , dist="exponential", data = DataLong)
M.dredge_noSize.Fso  <- dredge(M.global_noSize.Fso, subset = c(excl1so))
#Best observation + anthropogenic forcing model
M.global_noSize.Fso.1 <- survreg(Surv(RelDistance, Status) ~ Aspect + Backshore + NumericSectionNumber + Population_50km + RoadDistance 
                                 + State + SubColour + Substrate + TotalDebris + Weather + Population_50km*RoadDistance, dist = "exponential", data = DataLong)
summary(M.global_noSize.Fso.1)

#Combine predicting physical and anthropogenic forcing covariates to get the best full model
M.global_noSize.Ftot <- survreg(Surv(RelDistance, Status) ~  Aspect + Backshore + MeanOnshoreStokes1Mo + MeanOnshoreWind1Wk 
                               + NumericSectionNumber + State + SubColour + Substrate + TotalDebris + Weather + Population_50km + RoadDistance 
                               + Population_50km*RoadDistance + NumericSectionNumber*Population_50km,
                               dist = "exponential", data = DataLong)
M.dredge_noSize.Ftot <- dredge(M.global_noSize.Ftot)
M.global_noSize.Ftot.1 <- survreg(Surv(RelDistance, Status) ~  Aspect + Backshore + MeanOnshoreStokes1Mo + MeanOnshoreWind1Wk 
                                  + NumericSectionNumber + NumericSectionNumber*Population_50km + State + SubColour + Substrate + TotalDebris + Weather + RoadDistance + Population_50km,
                                  dist = "exponential", data = DataLong)
summary(M.global_noSize.Ftot.1)

## Forced survival model for those sections where debris was found (Status = 1) ---------
DataLongPositives <- DataLong[DataLong$Status == 1,] #Again, so wind and stokes component get included

#Physical forcing
numvar2ph <- c("MeanOnshoreWind1Mo","MeanOnshoreWind1Wk","MeanOnshoreWind1Dy","MeanOnshoreWind10Dy", "MeanOnshoreWindHyr",
             "MeanOnshoreStokes1Mo","MeanOnshoreStokes1Wk","MeanOnshoreStokes1Dy", "MeanOnshoreStokes10Dy", "MeanOnshoreStokesHyr",
             "TotalDebris", "NumericSectionNumber")
excl2ph <- Dredgeterms(numvar2ph,DataLongPositives)
M.global_inclSize.Fph <- survreg(Surv(RelDistance, Status) ~  LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber + MeanOnshoreWind1Mo
                             + MeanOnshoreWind1Wk + MeanOnshoreWind1Dy + MeanOnshoreWind10Dy + MeanOnshoreWindHyr + MeanOnshoreStokes1Mo
                             + MeanOnshoreStokes1Wk + MeanOnshoreStokes1Dy + MeanOnshoreStokes10Dy + MeanOnshoreStokesHyr, dist="exponential", data = DataLongPositives)
M.dredge_inclSize.Fph  <- dredge(M.global_inclSize.Fph, subset = c(excl2ph))
#Best observation + physical forcing model
M.global_inclSize.F.1ph <- survreg(Surv(RelDistance, Status) ~  LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber, dist="exponential", data = DataLongPositives)

#Anthropogenic forcing
numvar2so <- c("Population_5km","Population_20km","Population_50km","RoadDistance",
              "TotalDebris", "NumericSectionNumber")
excl2so <- Dredgeterms(numvar2so,DataLongPositives)
M.global_inclSize.Fso <- survreg(Surv(RelDistance, Status) ~ LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber + Population_5km + Population_20km + Population_50km 
                             + Population_50km*RoadDistance + RoadDistance, dist="exponential", data = DataLongPositives)
M.dredge_inclSize.Fso  <- dredge(M.global_inclSize.Fso, subset = c(excl2so))
#Best observation + anthropogenic forcing model
M.global_inclSize.F.1so <- survreg(Surv(RelDistance, Status) ~  LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber + Population_50km, dist="exponential", data = DataLongPositives)

#Combine predicting physical and anthropogenic forcing covariates to get the best full model
M.global_inclSize.F.tot <- survreg(Surv(RelDistance, Status) ~  LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber + Population_50km + Population_50km*NumericSectionNumber, dist="exponential", data = DataLongPositives)
M.dredge_inclSize.F.tot <- dredge(M.global_inclSize.F.tot)
M.global_inclSize.F.1tot <- survreg(Surv(RelDistance) ~  LowerBSize + NumericSectionNumber + TotalDebris + LowerBSize*NumericSectionNumber + Population_50km, dist="exponential", data = DataLongPositives)
summary(M.global_inclSize.F.1tot)

## Predictions ----------------------
NewData <- cbind(rep(c(0,1.2,2.4,4.9,9.6,19.2), times = 10), rep(1:10, each = 6), rep(median(DataLongPositives$TotalDebris), times = 60), rep(median(DataLongPositives$Population_50km), times = 60))
NewData <- as.data.frame(NewData)
names(NewData) <- c("LowerBSize","NumericSectionNumber","TotalDebris","Population_50km")
try.pred <- predict(M.global_inclSize.F.1tot, newdata = NewData, se.fit = T, type = c("response"))
NewData$RelDist <- try.pred$fit
NewData$SERelDist <- try.pred$se.fit
NewData$RelDistTrans <- NewData$RelDist + NewData$NumericSectionNumber - 1
NewData$DistMid <- NewData$RelDist - 0.5
