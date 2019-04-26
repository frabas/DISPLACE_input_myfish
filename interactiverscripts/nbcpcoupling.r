###################################################################################
#                                                                                 #
#                    Coupling spatio-temporal model to DISPLACE                   #  
#                                                                                 #
###################################################################################

# FIXME: Adapt script to automatically recognize the environment in which the model results are loaded!!!



# SECTION 1: Load model and model results 
# SECTION 2: Keep original abundance fields from model 
# SECTION 3: Predict abundance field one time-step ahead (core section for DISPLACE coupling)
# SECTION 4: Exract abundances for DISPLACE grid
# SECTION 5: Output for DISPLACE



 # CALLED FROM DISPLACE when dyn_pop_sce.option(Options::nbcpCoupling)
 args <- commandArgs(trailingOnly = TRUE)

 ##args is now a list of character vectors
 ## First check to see if arguments are passed.
 ## Then cycle through each element of the list and evaluate the expressions.
 if(length(args)==0){
    print("No arguments supplied to nbcpCoupling.r: Take pop 0 and tstep 0 and sce nbcpcoupling and simu simu1")
    pop     <- 0
    tstep   <- 745
    sce     <- "nbcpcoupling"
    sim     <- "simu1"
    igraph  <- 56
 }else{
   pop      <- args[1]
   tstep    <- args[2]
   sce      <- args[3]
   sim      <- args[4]
   igraph   <- args[5]
}

 
 
  application           <- "myfish"
  
  FRANCOIS <- TRUE ; MARIECHRISTINE <- FALSE
  
  if(.Platform$OS.type == "unix") {
    if(FRANCOIS)       path <- file.path("~","ibm_vessels", paste0("DISPLACE_input_", application), "interactiverscripts")
    if(MARIECHRISTINE) path <- file.path("H:","FB_MR","Coupling","old_coupling")  
  }
  if(.Platform$OS.type == "windows") {
    if(FRANCOIS)       path <- file.path("C:","Users","fbas","Documents","GitHub", paste0("DISPLACE_input_", application), "interactiverscripts")
    if(MARIECHRISTINE) path <- file.path("H:","FB_MR","Coupling","old_coupling")  
}

cat(paste("looking at the R script in ", path, "\n"))


if(FALSE){


#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 1                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.1) Loading libraries & functions, set working directory and load the results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(TMB)
library(gridConstruct)
library(Matrix)
library(fields)
library(raster)
library(tidyr)


#~~~~~~~~~~~~~~~~~
# 1.2) Load model 
#~~~~~~~~~~~~~~~~~
#compile(file.path(path, "model.cpp")) # Compile model only when running script for the first time
dyn.load(dynlib(file.path(path, "model"))) # Load the same C++ model used to generate the loaded results!


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.3) Load model results for each size group 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FIXME: the loaded model is saved in "env2" because it corresponds to model 2 in the LGCP script;
# The best chosen model will note be necessarily no.2, but could be no.4 instead and thus will be saved
# in the form of "env4". ADAPT SCRIPT TO AUTOMATICALY RECOGNIZE WHICH ENV IS LOADED.


# SizeGroup 0 (0-5cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "_toy_resTest_SG0.RData"))) 
SG0 <- env2
rm(list=setdiff(ls(), c("SG","gr")))


# SizeGroup 1 (5-10cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG1.RData"))) 
SG1 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","gr")))


# SizeGroup 2 (10-15cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, "toy_resTest_SG2.RData"))) 
SG2 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","gr")))


# SizeGroup 3 (15-20cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG3.RData"))) 
SG3 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","gr")))


# SizeGroup 4 (20-25cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG4.RData"))) 
SG4 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","gr")))


# SizeGroup 5 (25-30cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG5.RData"))) 
SG5 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","gr")))


# SizeGroup 6 (30-35cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG6.RData"))) 
SG6 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","gr")))


# SizeGroup 7 (35-40cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG7.RData"))) 
SG7 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","gr")))


# SizeGroup 8 (40-45cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG8.RData"))) 
SG8 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","gr")))


# SizeGroup 9 (45-50cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG9.RData"))) 
SG9 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","SG9","gr")))


# SizeGroup 10 (50-55cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG10.RData"))) 
SG10 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","SG9","SG10","gr")))


# SizeGroup 11 (55-60cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG11.RData"))) 
SG11 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","SG9","SG10","SG11","gr")))


# SizeGroup 12 (60-65cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG12.RData"))) 
SG12 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","SG9","SG10","SG11","SG12","gr")))


# SizeGroup 13 (65+cm)
#~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path, paste0(pop, "toy_resTest_SG13.RData"))) 
SG13 <- env2
rm(list=setdiff(ls(), c("SG0","SG1","SG2","SG3","SG4","SG5","SG5","SG7","SG8","SG9","SG10","SG11","SG12","SG13","gr")))


snames <- ls(pattern="SG")
Sizes <-  lapply(ls(pattern="SG"), function(x) get(x)) #Include all objects into a list
names(Sizes) <- c(snames)



######################                END OF SECTION 1            #######################






#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 2                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1) Retrieving objects of the model (to be used along the script)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obj <- list()
for(i in seq_along(Sizes)){
  obj[[i]] <- Sizes[[i]]$obj # Retrieve obj from model
  obj[[i]]$fn(Sizes[[i]]$fit$par) # Evaluate again the obj function to retrieve information that were lost during the result saving process 
}

lpb <- list()
for(i in seq_along(obj)){
  lpb[[i]] <- obj[[i]]$env$last.par.best #Get spatio-temporal parameters
}

r <- list()
h <- list()
for(i in seq_along(obj)){
  for(p in seq_along(lpb)){
    r[[i]] <- obj[[i]]$env$random # Get latent (random) variables
    ##obj$retape()
    h[[i]] <- obj[[i]]$env$spHess(lpb[[p]], random=TRUE) #Hessian matrix based on the spatio-temporal models;This is where we need the loaded model.cpp file (done in section 1.2)
    ##image(h)
  }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.2) Retrieve predicted abundance fields from the models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
abundanceSG0 <- as.data.frame(as.list(SG0$sdr,"Estimate")[1]); abundanceSG0 <- apply(abundanceSG0,2,exp)#Transform to natural scale
abundanceSG1 <- as.data.frame(as.list(SG1$sdr,"Estimate")[1]); abundanceSG1 <- apply(abundanceSG1,2,exp)#Transform to natural scale
abundanceSG2 <- as.data.frame(as.list(SG2$sdr,"Estimate")[1]); abundanceSG2 <- apply(abundanceSG2,2,exp)#Transform to natural scale
abundanceSG3 <- as.data.frame(as.list(SG3$sdr,"Estimate")[1]); abundanceSG3 <- apply(abundanceSG3,2,exp)#Transform to natural scale
abundanceSG4 <- as.data.frame(as.list(SG4$sdr,"Estimate")[1]); abundanceSG4 <- apply(abundanceSG4,2,exp)#Transform to natural scale
abundanceSG5 <- as.data.frame(as.list(SG5$sdr,"Estimate")[1]); abundanceSG5 <- apply(abundanceSG5,2,exp)#Transform to natural scale
abundanceSG6 <- as.data.frame(as.list(SG6$sdr,"Estimate")[1]); abundanceSG6 <- apply(abundanceSG6,2,exp)#Transform to natural scale
abundanceSG7 <- as.data.frame(as.list(SG7$sdr,"Estimate")[1]); abundanceSG7 <- apply(abundanceSG7,2,exp)#Transform to natural scale
abundanceSG8 <- as.data.frame(as.list(SG8$sdr,"Estimate")[1]); abundanceSG8 <- apply(abundanceSG8,2,exp)#Transform to natural scale
abundanceSG9 <- as.data.frame(as.list(SG9$sdr,"Estimate")[1]); abundanceSG9 <- apply(abundanceSG9,2,exp)#Transform to natural scale
abundanceSG10 <- as.data.frame(as.list(SG10$sdr,"Estimate")[1]); abundanceSG10 <- apply(abundanceSG10,2,exp)#Transform to natural scale
abundanceSG11 <- as.data.frame(as.list(SG11$sdr,"Estimate")[1]); abundanceSG11 <- apply(abundanceSG11,2,exp)#Transform to natural scale
abundanceSG12 <- as.data.frame(as.list(SG12$sdr,"Estimate")[1]); abundanceSG12 <- apply(abundanceSG12,2,exp)#Transform to natural scale
abundanceSG13 <- as.data.frame(as.list(SG13$sdr,"Estimate")[1]); abundanceSG13 <- apply(abundanceSG13,2,exp)#Transform to natural scale


# Set all abundances dfs into a list
Allsinit <- list(SG0=abundanceSG0,SG1=abundanceSG1,SG2=abundanceSG2,SG3=abundanceSG3,SG4=abundanceSG4,SG5=abundanceSG5,SG6=abundanceSG6,
                 SG7=abundanceSG7,SG8=abundanceSG8,SG9=abundanceSG9,SG10=abundanceSG10,SG11=abundanceSG11,SG12=abundanceSG12,SG13=abundanceSG13)
  

# Set colnames based on timesteps
for(i in seq_along(Allsinit )){
    colnames(Allsinit [[i]]) <- substr(levels(SG0$data$time), start=6,stop=6) #Choose randomly one model output; works for any output (SG0 or SG5,etc.)
  }


#image(gr, Allsinit[[2]][,1]) #To see the predicted abundance fields




######################                END OF SECTION 2            #######################





#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 3                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>

# Here we predict the abundance field one time-step ahead (t+1) for all 4 quarters of the last year.
# This is based on the precited abundance fields that were retrieved in SECTION 2.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.1) Retrieve spatio-temporal parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
time_corr <- list() 
phi <- list() 
delta <- list() 
scale <- list()
Q <- list() # Precision matrix (inverse of covariance matrix)


# Transform parameters into their natural scale
for (i in seq_along(lpb)){
  time_corr[[i]] <- lpb[[i]]["time_corr"] #untransformed time-correlation parameter
  phi[[i]] <- time_corr[[i]] / sqrt(1.0 + time_corr[[i]]*time_corr[[i]]) # phi (time correlation param.)
  delta[[i]] <- exp(lpb[[i]]["logdelta"]) # delta (spatial correlation param.)
  scale[[i]] <- exp(lpb[[i]]["logscale"]) # scale (spatial correlation param)
  Q[[i]] <- obj[[i]]$env$data$Q0 + delta[[i]] * obj[[i]]$env$data$I # Precision matrix
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.2) Predict one-time step ahead
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NROW <- dim(Allsinit[[1]])[1] #We can take any element of the list, as they have all the same dimension
NCOL <- dim(Allsinit[[1]])[2] #We can take any element of the list, as they have all the same dimension
DIM <- length(Allsinit)

tmp <- array(1:(NCOL*NROW), c(NROW,NCOL,DIM)) 

# Simulate in t+1 - Note that in the rmvnorm_prec function one can choose as many simulation as desired (n.sims argument)
for(i in seq_along(Allsinit)){
  for(j in 1:NCOL){
    Abundance <- Allsinit[[i]][,j]
    tmp[,j,i] <- rmvnorm_prec(mu = Abundance, prec = scale[[i]]*sqrt(1-phi[[i]]^2)*Q[[i]], n.sims = 1) #Simulate and store results in the array
    Allspred <- lapply(seq(dim(tmp)[3]), function(x) tmp[ , , x]) # Convert array back to list (easier to maniupulate)
  }
}


names(Allspred) <- names(Sizes)
# image(gr, Allspred[[1]][,1]) #plotting example



######################                END OF SECTION 3            #######################





#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 4                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>



# First we need to create individual raster files, as the abundance output given in 
# SECTION 3 isn't a raster file. After that we use the lon/lat of the DISPLACE graphs 
# and extract the abundances according to each graph node

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.1) Create an empty rasterfile 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
e <- extent(as.matrix(gr))
r <- raster(e,ncol=80,nrow=80) #Adapt raster size according to grid size; These numbers are good for 5x5km grid (as model output)


# Include grid positions to dataframes
coords <- data.frame(lon=gr$lon,lat=gr$lat) #Retrieve grid positions from model output
Allspred <- lapply(Allspred, cbind, coords) #Bind grid positions


# Convert df from wide to long format
Allspredw <- list(NULL)
for(i in seq_along(Allspred)){
  Allspredw[[i]] <- gather(Allspred[[i]], key=Quarter,value=Abundance,-c(lon,lat))
}
names(Allspredw) <- names(Allspred) #name new list accordingly


# Include size-specific column to each df
for( i in seq_along(Allspredw)){
  Allspredw[[i]]$Size <- rep(names(Allspredw)[i],nrow(Allspredw[[i]]))
}

Allspredw <- do.call("rbind", Allspredw) #Unlist
rownames(Allspredw) <- NULL


# Identify properly all quarters
Allspredw$Quarter <- ifelse(Allspredw$Quarter=="1","Q1",
                            ifelse(Allspredw$Quarter=="2","Q2",
                                   ifelse(Allspredw$Quarter=="3","Q3","Q4")))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.2) Create raster based on abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dflist <- split(Allspredw,list(Allspredw$Quarter,Allspredw$Size))

abufields <- list() 
for(i in 1:length(dflist)){
  abufields[[i]] <- rasterize(dflist[[i]][,c("lon","lat")], r, dflist[[i]][,"Abundance"], fun=mean)
  abufields[[i]] <- disaggregate(abufields[[i]],2,method="bilinear")
}
names(abufields) <- names(dflist) #name new list accordingly


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4.3) Extract abundances based on DISPLACE graph coords
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loc <- read.table(file=file.path(path, paste0("coord", igraph,".dat"),sep="")) # DISPLACE coords
loc <- as.matrix(as.vector(loc))
loc <- matrix(loc, ncol=3)
loc <- cbind(loc, 1:nrow(loc))
colnames(loc) <- c('lon', 'lat', 'harb', 'pt_graph')

loc <- as.data.frame(loc)
loc <- loc[loc$harb==0,] # remove points on harbours (all harb>0)
loc <- loc[,c("lon","lat")] # keep only lon/lat columns


dfa <- list()
for(i in seq_along(abufields)){
  dfa[[i]]<- raster::extract(abufields[[i]],loc)
  dfa[[i]]<- data.frame(abundance=dfa[[i]], lon=loc$lon, lat=loc$lat)
}

names(dfa) <- names(abufields)



######################                END OF SECTION 4            #######################




#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 5                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>

# Organize predicted abundances into suitable format for DISPLACE
# The output will be a dataframe cotaining the following columns:
# DISPLACE lon $ lat, Quarter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.1) Include Quarter and Size columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(dfa)){
  dfa[[i]]$Quarter <- factor(substr(names(dfa)[i], start = 1, stop = 2)) # Doubble check if it's correct!!
  dfa[[i]]$Size <- factor(substr(names(dfa)[i], start=4, stop=5)) # Doubble check if it's correct!!
}

dfa <- do.call("rbind",dfa) 
rownames(dfa) <- NULL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.2) Transform to wide format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
displace_dat <- spread(dfa, Size, abundance)


#~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.3) Normalize abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~
NormAbu <- function (x) {
  x <- x/sum(x)
  x
}

# DOUBBLE CHECK THIS!!!
displace_dat[,c(4:ncol(displace_dat))] <- apply(displace_dat[,c(4:ncol(displace_dat))],2,NormAbu)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5.4) Save predcited resuls
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# CHECK!!!
# There are many DISPLACE coords that were outside the range of the predicted abundances by the model which, in turn,
# resulted to NA in the "displace_dat" data.frame. Check wheter these rows should be removed or set to 0 prior to result saving.

saveRDS(displace_dat,"AbuDISPLACE_WBSCod.rds")

} # end FALSE

cat("coupling to nbcp...done\n")

