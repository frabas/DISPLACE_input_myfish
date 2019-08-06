###################################################################################
#                                                                                 #
#                    Coupling spatio-temporal model to DISPLACE                   #
#                             Part 1: data pre-processing                         #
#                                                                                 #
###################################################################################


# Process model results to retain only the most important objects that will be used
# for the NBCP-DISPLACE coupling (script 2)
# This includes: 
# * spatio-tempral correlation parameters (phi, delta, scale)
# * Precision matrix (Q)
# * time periods (to identify the last time-period)
# * estimated abundances for each time period (transformed to natural scale)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) Loading libraries & set working directory 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(TMB)
library(gridConstruct)
library(Matrix)
library(fields)
library(raster)
library(tidyr)


#setwd("H:/FB_MR/Results/Cod/WBS") #Full results to be uploaded
setwd("C:/Users/mruf/OneDrive - Danmarks Tekniske Universitet/Results_WBScod/DISPLACE")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2) Load model results for each size group and set into a list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIXME: the loaded model is saved in "env2" because it corresponds to model 2 in the LGCP script;
# The best chosen model will note be necessarily no.2, but could be no.4 instead and thus will be saved
# in the form of "env4". ADAPT SCRIPT TO AUTOMATICALY RECOGNIZE WHICH ENV IS LOADED.


# SizeGroup 1 (0-5cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S1.RData") 
SG0 <- env2
time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG0 <- list() # Put everyhing into a list
lSG0$abundance <- as.data.frame(as.list(SG0$sdr, "Estimate")[1]); lSG0$abundance <- apply(lSG0$abundance,2,exp) #Put on natural scale
lSG0$gr <- gr
lSG0$time_period <-SG0$data$time

lSG0$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG0$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG0$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG0$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0")))



# SizeGroup 2 (5-10cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S2.RData") 
SG1 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG1 <- list() # Put everyhing into a list
lSG1$abundance <- as.data.frame(as.list(SG1$sdr, "Estimate")[1]); lSG1$abundance <- apply(lSG1$abundance,2,exp) #Put on natural scale
lSG1$gr <- gr
lSG1$time_period <-SG1$data$time

lSG1$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG1$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG1$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG1$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1")))


# SizeGroup 3 (10-15cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S3.RData") 
SG2 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG2 <- list() # Put everyhing into a list
lSG2$abundance <- as.data.frame(as.list(SG2$sdr, "Estimate")[1]); lSG2$abundance <- apply(lSG2$abundance,2,exp) #Put on natural scale
lSG2$gr <- gr
lSG2$time_period <-SG2$data$time

lSG2$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG2$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG2$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG2$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2")))


# SizeGroup 4 (15-20cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S4.RData") 
SG3 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG3 <- list() # Put everyhing into a list
lSG3$abundance <- as.data.frame(as.list(SG3$sdr, "Estimate")[1]); lSG3$abundance <- apply(lSG3$abundance,2,exp) #Put on natural scale
lSG3$gr <- gr
lSG3$time_period <-SG3$data$time

lSG3$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG3$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG3$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG3$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3")))



# SizeGroup 5 (20-25cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S5.RData") 
SG4 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG4 <- list() # Put everyhing into a list
lSG4$abundance <- as.data.frame(as.list(SG4$sdr, "Estimate")[1]); lSG4$abundance <- apply(lSG4$abundance,2,exp) #Put on natural scale
lSG4$gr <- gr
lSG4$time_period <-SG4$data$time

lSG4$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG4$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG4$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG4$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4")))



# SizeGroup 6 (25-30cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S6.RData") 
SG5 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG5 <- list() # Put everyhing into a list
lSG5$abundance <- as.data.frame(as.list(SG5$sdr, "Estimate")[1]); lSG5$abundance <- apply(lSG5$abundance,2,exp) #Put on natural scale
lSG5$gr <- gr
lSG5$time_period <-SG5$data$time

lSG5$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG5$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG5$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG5$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5")))



# SizeGroup 7 (30-35cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S7.RData") 
SG6 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG6 <- list() # Put everyhing into a list
lSG6$abundance <- as.data.frame(as.list(SG6$sdr, "Estimate")[1]); lSG6$abundance <- apply(lSG6$abundance,2,exp) #Put on natural scale
lSG6$gr <- gr
lSG6$time_period <-SG6$data$time

lSG6$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG6$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG6$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG6$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6")))



# SizeGroup 8 (35-40cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S8.RData") 
SG7 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG7 <- list() # Put everyhing into a list
lSG7$abundance <- as.data.frame(as.list(SG7$sdr, "Estimate")[1]); lSG7$abundance <- apply(lSG7$abundance,2,exp) #Put on natural scale
lSG7$gr <- gr
lSG7$time_period <-SG7$data$time

lSG7$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG7$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG7$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG7$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7")))


# SizeGroup 9 (40-45cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S9.RData") 
SG8 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG8 <- list() # Put everyhing into a list
lSG8$abundance <- as.data.frame(as.list(SG8$sdr, "Estimate")[1]); lSG8$abundance <- apply(lSG8$abundance,2,exp) #Put on natural scale
lSG8$gr <- gr
lSG8$time_period <-SG8$data$time

lSG8$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG8$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG8$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG8$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8")))



# SizeGroup 10 (45-50cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S10.RData") 
SG9 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG9 <- list() # Put everyhing into a list
lSG9$abundance <- as.data.frame(as.list(SG9$sdr, "Estimate")[1]); lSG9$abundance <- apply(lSG9$abundance,2,exp) #Put on natural scale
lSG9$gr <- gr
lSG9$time_period <-SG9$data$time

lSG9$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG9$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG9$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG9$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8","lSG9")))


# SizeGroup 11 (50-55cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S11.RData") 
SG10 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG10 <- list() # Put everyhing into a list
lSG10$abundance <- as.data.frame(as.list(SG10$sdr, "Estimate")[1]); lSG10$abundance <- apply(lSG10$abundance,2,exp) #Put on natural scale
lSG10$gr <- gr
lSG10$time_period <-SG10$data$time

lSG10$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG10$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG10$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG10$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8","lSG9","lSG10")))


# SizeGroup 12 (55-60cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S12.RData") 
SG11 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG11 <- list() # Put everyhing into a list
lSG11$abundance <- as.data.frame(as.list(SG11$sdr, "Estimate")[1]); lSG11$abundance <- apply(lSG11$abundance,2,exp) #Put on natural scale
lSG11$gr <- gr
lSG11$time_period <-SG11$data$time

lSG11$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG11$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG11$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG11$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8","lSG9","lSG10","lSG11")))


# SizeGroup 13 (60-65cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S13.RData") 
SG12 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG12 <- list() # Put everyhing into a list
lSG12$abundance <- as.data.frame(as.list(SG12$sdr, "Estimate")[1]); lSG12$abundance <- apply(lSG12$abundance,2,exp) #Put on natural scale
lSG12$gr <- gr
lSG12$time_period <-SG12$data$time

lSG12$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG12$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG12$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG12$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8","lSG9","lSG10","lSG11","lSG12")))


# SizeGroup 14 (65+cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load("cod-both-m2_S14.RData") 
SG13 <- env2

time_corr <- lpb["time_corr"] #untransformed time-correlation parameter
delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)

lSG13 <- list() # Put everyhing into a list
lSG13$abundance <- as.data.frame(as.list(SG13$sdr, "Estimate")[1]); lSG13$abundance <- apply(lSG13$abundance,2,exp) #Put on natural scale
lSG13$gr <- gr
lSG13$time_period <-SG13$data$time

lSG13$phi <- time_corr / sqrt(1.0 + time_corr*time_corr) # phi (time correlation param.)
lSG13$delta <- exp(lpb["logdelta"]) # delta (spatial correlation param.)
lSG13$scale <- exp(lpb["logscale"]) # scale (spatial correlation param)
lSG13$Q <- obj$env$data$Q0 + delta * obj$env$data$I # Precision matrix

rm(list=setdiff(ls(), c("lSG0","lSG1","lSG2","lSG3","lSG4","lSG5","lSG6","lSG7","lSG8","lSG9","lSG10","lSG11","lSG12","lSG13")))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3) Make a list of all size-group lists
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fullres <- list(SG0=lSG0,SG1=lSG1,SG2=lSG2,SG3=lSG3,SG4=lSG4,SG5=lSG5,SG6=lSG6,
                 SG7=lSG7,SG8=lSG8,SG9=lSG9,SG10=lSG10,SG11=lSG11,SG12=lSG12,SG13=lSG13)


# Set colnames based on timesteps
for(i in seq_along(fullres)){
  colnames(fullres[[i]]$abundance) <- levels(fullres[[i]]$time_period) 
}


# Save full list (to be used in script 2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save(fullres,file="H:/FB_MR/Coupling/WBScod.RData") 
