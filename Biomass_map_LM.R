#### Activate the libraries#####

library(raster)
library(rgdal)
library(rgeos)
library(nlme)
library(ggplot2)


#### set the path#####
path1 = "C:/Users/Senait Alazar/Desktop/day 6/ForTheStudents"
path = "C:/Users/Senait Alazar/Desktop/day 6/ForTheStudents/"


## importing the data####
DTM <- raster(paste(path, "DTM_2007_UTM32N.tif", sep=""), band=1) ### Terrain model ###
DSM <- raster(paste(path, "DSM_2007_UTM32N.tif", sep=""), band=1) ### Surface model ###


#### Rapid Eye bands (B1=blue, B2= green, B3=red, B4=red edge, B5=Near infrared) #####
RAP1 <- raster(paste(path, "KA_RAP_2013_5_18_UTM32N_res.tif", sep="" ), band=1)
RAP2 <- raster(paste(path, "KA_RAP_2013_5_18_UTM32N_res.tif", sep="" ), band=2)
RAP3 <- raster(paste(path, "KA_RAP_2013_5_18_UTM32N_res.tif", sep="" ), band=3)
RAP4 <- raster(paste(path, "KA_RAP_2013_5_18_UTM32N_res.tif", sep="" ), band=4)
RAP5 <- raster(paste(path, "KA_RAP_2013_5_18_UTM32N_res.tif", sep="" ), band=5)


### Creation of nDSM layer ####
nDSM <- DSM-DTM

### Create NDVI layer #### we use the red edge band [B4] and NIF
ndvi <- (RAP5-RAP4)/(RAP5+RAP4)

### create composite
RAP <- stack(RAP1, RAP2, RAP3, RAP4, RAP5, ndvi, nDSM)

### import the field inventory data--we use ReadOGR to import vector data
cluster_shp <- readOGR(dsn = path1, layer= "Cluster_with_biomass_KA_2014_UTM")

###extract the training data--creating buffer--create polygon [best option] and use the Hmean value  

GTdata <- extract(RAP, cluster_shp, buffer= 35, fun=mean)
head(GTdata)


### Create a data frame with all variables [response[which is the biomass] + explanatory]
biomass<- as.data.frame(cluster_shp$bm_cl_tha)
biomass

GTdata <- cbind(GTdata, biomass)### rewrite the GTdata by adding one layer that is the biomass
GTdata
head(GTdata)

### renaming...###
names(GTdata) <- c("mean_Blue","mean_Green","mean_Red","mean_RedEdge","mean_NIR","mean_NDVI","mean_nDSM","Biomass")

head(GTdata)
summary(GTdata)


#################### MODELING STEPS" ###################################
########1- Data Exploration ###################

#1a- outlier detection ----checking the data
x11()
op <- par(mfrow=c(4,2))
dotchart(GTdata$mean_Blue) #### Blue band got weird--not good to model with blue
dotchart(GTdata$mean_Green)
dotchart(GTdata$mean_Red)
dotchart(GTdata$mean_RedEdge)
dotchart(GTdata$mean_NIR)
dotchart(GTdata$mean_NDVI)
dotchart(GTdata$mean_nDSM)
dotchart(GTdata$Biomass)
par(op)
dev.off() 

rm(op)
##### no issue is observed, next step #######

#1b- collinearity-- need to drop one of the variables that are strongly related-- use one to predict the other-drop one the affected variables---[pair plot, correlation coefficient and variant inflation factor- used to detect collinearity]

#pair plot-- 
pairs(GTdata)
x11()

# correlation
cor(GTdata)

######## 2- Fitting the model/ Model selection #############

mod_lm1 <- lm(Biomass ~ 
                mean_Blue+
                mean_Green+
                mean_Red+
                mean_RedEdge+
                mean_NIR+
                mean_NDVI+
                mean_nDSM, data=GTdata)
## Model summary
summary(mod_lm1)

## drop the blue band as it got nothing to do with vegetation

mod_lm2 <- lm(Biomass ~ 
                mean_Green+
                mean_Red+
                mean_RedEdge+
                mean_NIR+
                mean_NDVI+
                mean_nDSM, data=GTdata)
#model1 summary
summary(mod_lm2)

AIC(mod_lm1,mod_lm2) ### the one with lower AIC value is the better model..mod_lm1 is better here..


##### step procedure to drop variables in a model-- using a step wise fitting

step(mod_lm1)

####### fit the best model###########

### from the results..lm(formula = Biomass ~ mean_Green + mean_RedEdge + mean_nDSM, data = GTdata)

mod_lm <- lm(Biomass ~
               mean_Green+
               mean_RedEdge+
               mean_nDSM, data=GTdata)

summary(mod_lm)

####### 3- MOdel Validation #############
## plot the model

#plot - check of homogeneity - try to identify a pattern 

op <- par(mfrow=c(2,2))
plot(mod_lm)
dev.off()


### Histogram check for normality: plotting a histogram of residuals
hist(resid(mod_lm), breaks=20, main=" Histogram of Model Residuals", xlab= "Residuals[model: mod_lm]")


### assumption of independence plot residuals against each explanatory variable used in the final model. If you see a pattern then there is a violation

plot(GTdata$mean_Green, resid(mod_lm))### no problem in this aspect
plot(GTdata$mean_RedEdge, resid(mod_lm))### no problem in this aspect
plot(GTdata$mean_nDSM, resid(mod_lm))## outlier because of  very high tree

### all the above validation suggest no serious violation assumptions--so all good!!

#### Goodness of fit-used for presenting your model in scientific study --to check the quality #####
biomass_m <- cluster_shp$bm_cl_tha
biomass_pred <- predict(mod_lm)

plot(biomass_m,biomass_pred, main=" Goodness of Fit: MLR Model", xlab="Measured Biomass [t/ha]", ylab=" Predicted Biomass [t/ha]", xlim= c(75, 350), ylim=c(100, 250)) ## goodness of the fit plot##

abline(1.1)
abline(lm(biomass_pred ~ biomass_m), col=2)## adding abline
abline(lm(biomass_m ~ biomass_pred), col=1)

cor(biomass_m,biomass_pred)## correlation of between measured and predicted

r2 <- cor(biomass_m,biomass_pred)*cor(biomass_m,biomass_pred)# coeff. of determination for the model goodness of fit


###### 4- Plotting the biomass map using our best fit model (mod_lm)

## create a stack with the raster bands/ layers from which the variables used in the final model were extracted
RAP_mod <-stack(RAP2, RAP4, nDSM) 
names(RAP_mod) <- c("mean_Green","mean_RedEdge","mean_nDSM")
writeRaster(predict(object=RAP_mod, model=mod_lm, na.rm=T), filename = paste(path, "b_RAP-L_lm_2.tif", sep=""), format="GTiff", overwrite= TRUE)
RAP_mod


