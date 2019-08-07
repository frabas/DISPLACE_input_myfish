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
    if(FRANCOIS){
           path <- file.path("~","ibm_vessels", paste0("DISPLACE_input_", application))
           }
    if(MARIECHRISTINE){
           path <- file.path("H:","FB_MR","Coupling")  
           }
  }
  if(.Platform$OS.type == "windows") {
    if(FRANCOIS){
           path <- file.path("C:","Users","fbas","Documents","GitHub", paste0("DISPLACE_input_", application))
           }
    if(MARIECHRISTINE){
          path <- file.path("H:","FB_MR","Coupling")  
          }
}

cat(paste("looking at the R script in ", path, "\n"))


if(TRUE){


###################################################################################
#                                                                                 #
#                    Coupling spatio-temporal model to DISPLACE                   #
#                             Part 2: Coupling results                            #
#                                                                                 #
###################################################################################


#library(TMB)
#library(gridConstruct)
suppressWarnings(suppressMessages(library(Matrix)))
suppressWarnings(suppressMessages(library(fields)))
suppressWarnings(suppressMessages(library(raster)))
suppressWarnings(suppressMessages(library(tidyr)))


# SECTION 1: Predict abundance field one time-step ahead (core section for DISPLACE coupling)
# SECTION 2: Extract abundances for DISPLACE grid
# SECTION 3: Output for DISPLACE



# Multivariate normal distribution simulation function (based on precision rather variance)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# @ mu = mean of the abundance field
# @ prec = precision of the spatio-temporal covariance matrix
# @ n.sims = number of simulations
rmvnorm_prec <- function(mu, prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- Cholesky(prec, super=TRUE)
  z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  mu + z
}





#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 1                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>


# Load results
#~~~~~~~~~~~~~~~~
if(pop ==10){
  load(file.path(path, "interactiverscripts", "WBScod.RData")) #fullres is a list of lists, where each individual list stores the results of a particular size-group
  cat("load data input to nbcp...done\n")
  } else{
  cat(paste("no nbcp input file available for this stock...\n"))
  stop()
  }
# names(fullres) #names of the lists(related to the DISPLACE size-groups)
# names(fullres[[1]]) #names of the objects present in each list
# image(fullres[[1]]$gr,fullres[[1]]$abundance[,1]) # To see an example of the estimated abundance fields




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.1) Retrieve spatio-temporal parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get abundances for the last time-step (2016, with all 4 quarters)
abundances <- vector("list", 14)
Q1 <- ncol(fullres[[1]]$abundance)-4+1
Q4 <- ncol(fullres[[1]]$abundance)
for(i in seq_along(fullres)){
  abundances[[i]] <- fullres[[i]]$abundance[,Q1:Q4]
}
names(abundances) <- names(fullres)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.2) Predict one-time step ahead
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To predict the abundance fields in t+1, we need to "draw" a new abundance field that is based
# on the spatio-temporal correlation parameters + precision matrix that were estimated for each stock-sizeGroup
# (retrieved in the coupling script part 1 (NBCP_DISPLACE_coupling_part01.R)).

# As such, we need to apply the rmvnorm_prec function, which basically simulates
# a multivariate Gaussian random distribution based on the input spatio-temporal corr. parameters and precision matrix.
# The function requires that the user specifies the number of simulations; usually,the more precise the abundance
# estimates, the more precise will be the forward prediction (thus, only a few simulations, or even only 1, would be required).
# However, when abundance estimates have high uncertainties (such as those of very small size-groups), the forward-predicted
# abundance fields tend to be very different from one simulation to another, and therefore a higher number of simulation
# is required to stabilize the predicted abundance field. Provided that it is a "half-simulation" framework,
# we can compute forward-predictions for all four quarters of the last time-period. 


Allspred <- vector("list", 14)


NROW <- dim(fullres[[1]]$abundance)[1]  # nb of locations
NCOL <- ncol(abundances[[1]])           # nb of quarters
DIM <- length(fullres)                  # nb of szgroups
tmp <- array(1:(NCOL*NROW), c(NROW,NCOL,DIM)) 

NSIM <- 100 #Choose the number of simulations and take the mean across the simulated fields



for(i in seq_along(abundances)){  # i is szgroup
  for(j in 1:NCOL){  # j is Quarter
    Abundance <- abundances[[i]][,j]
    #tmp[,j,i] <- rmvnorm_prec(mu = Abundance, prec = fullres[[i]]$scale*sqrt(1-fullres[[i]]$phi^2)*fullres[[i]]$Q, n.sims = 1) #Simulate and store results in the array
    tmpres <- as.data.frame(rmvnorm_prec(mu = Abundance, prec = fullres[[i]]$scale*sqrt(1-fullres[[i]]$phi^2)*fullres[[i]]$Q, n.sims = NSIM)) #Simulate and store results in the array
    tmpres <- rowMeans(tmpres[,1:NSIM])
    tmp[,j,i]  <- tmpres
    Allspred <- lapply(seq(dim(tmp)[3]), function(x) tmp[ , , x]) # Convert array back to list (easier to manipulate)
  }
}
names(Allspred) <- names(fullres)
# image(fullres[[1]]$gr, Allspred[[1]][,1],col=tim.colors(99)) #plotting example

cat("simulate a nbcp pop ",pop," abundance field for t+1...done\n")



######################                END OF SECTION 1            #######################





#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 2                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>



# First we need to create individual raster files, as the abundance output given in 
# SECTION 1 isn't a raster file. After that we use the lon/lat of the DISPLACE graphs 
# and extract the abundances according to each graph node


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1) Create an empty rasterfile 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gr <- fullres[[1]]$gr # Get the grid from an arbitrary size-class (gr is the same throughout all size-classes)
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
# 2.2) Create raster based on abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dflist <- split(Allspredw,list(Allspredw$Quarter,Allspredw$Size))

abufields <- list() 
for(i in 1:length(dflist)){
  abufields[[i]] <- rasterize(dflist[[i]][,c("lon","lat")], r, dflist[[i]][,"Abundance"], fun=mean)
  abufields[[i]] <- disaggregate(abufields[[i]],2,method="bilinear")
}
names(abufields) <- names(dflist) #name new list accordingly


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.3) Extract abundances based on DISPLACE graph coords
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
loc <- read.table(file=file.path(path, "graphsspe", paste("coord",igraph,".dat",sep=""))) # DISPLACE coords
loc <- as.matrix(as.vector(loc))
loc <- matrix(loc, ncol=3)
loc <- cbind(loc, 1:nrow(loc))
colnames(loc) <- c('lon', 'lat', 'harb', 'pt_graph')

loc <- as.data.frame(loc)
loc <- loc[loc$harb==0,] # remove points on harbours (all harb>0)
loc2 <- loc[,c("lon","lat")] # keep only lon/lat columns



if(pop==10){
# clip to the real stock area definition
clip <- raster(file.path(path, "interactiverscripts", "sd222324.tif"))
for(i in seq_along(abufields)){
   abufields[[i]] <- crop(abufields[[i]], clip)
   }
}






dfa <- list()
for(i in seq_along(abufields)){
  dfa[[i]]<- raster::extract(abufields[[i]],loc2)
  dfa[[i]]<- data.frame(abundance=dfa[[i]], lon=loc2$lon, lat=loc2$lat,pt_graph=loc$pt_graph)
}

names(dfa) <- names(abufields)







######################                END OF SECTION 2            #######################




#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>
########################                SECTION 3                ########################
#><><><><><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><><><>><>

# Organize predicted abundances into suitable format for DISPLACE
# The output will be a dataframe cotaining the following columns:
# DISPLACE lon $ lat, Quarter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.1) Include Quarter and Size columns
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in seq_along(dfa)){
  dfa[[i]]$Quarter <- factor(substr(names(dfa)[i], start = 1, stop = 2)) # Double check if it's correct!!
  dfa[[i]]$Size <- factor(substr(names(dfa)[i], start=4, stop=7)) # Double check if it's correct!!
}

dfa <- do.call("rbind",dfa) 
rownames(dfa) <- NULL

dfa$Size <- factor(dfa$Size, levels=paste0("SG", 0:13))  # reorder


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.2) Transform to wide format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
displace_dat <- spread(dfa, Size, abundance)



#~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.3) Useful to later Normalize abundances
#~~~~~~~~~~~~~~~~~~~~~~~~~~

# Normalize function
#~~~~~~~~~~~~~~~~~~~~
#NormAbu <- function(x) {
#  return ((x - min(x)) / (max(x) - min(x)))
#}
NormAbu <- function(x){
  return(exp(x)/sum(exp(x)))
}


idxna <- which(is.na(displace_dat$SG0)) #Can choose any arbitrary size-group column as the result will be the same
displace_dat <- displace_dat[-idxna,]





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.4) Save predicted results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 #saveRDS(displace_dat,"AbuDISPLACE_WBSCod.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3.4) Convert to DISPLACE POPULATIONS\static_avai format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # please respect the file naming, also including "_updated" in the name.
 
 the_selected_szgroups <- c(0, 2, 3, 5, 7)
 
 annual_tstep <- as.numeric(tstep) %% 8761
 quarter <- "Q1"
 if(annual_tstep>=2160 && annual_tstep<4344) quarter <- "Q2"
 if(annual_tstep>=4344 && annual_tstep<6552) quarter <- "Q3"
 if(annual_tstep>=6552 && annual_tstep<8761) quarter <- "Q4"
 
 if(quarter=="Q1" || quarter=="Q2"){
   if(quarter=="Q1"){
       dd      <- displace_dat[displace_dat$Quarter=="Q1",]
       dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
       print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
       ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
       library(doBy)
       ddd      <- orderBy(~pt_graph, ddd)
       dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
       dat_full <- ddd[, c("pt_graph", "avai")]
   } else{
       dd      <- displace_dat[displace_dat$Quarter=="Q2",]
       dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
       print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
       ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
       library(doBy)
       ddd      <- orderBy(~pt_graph, ddd)
       dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
       dat_full <- ddd[, c("pt_graph", "avai")]
     } 
   write.table(dat, file=file.path(path, paste0("popsspe_", application), "static_avai",
                          paste(pop, "spe_avai_szgroup_nodes_semester1_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
   write.table(dat_full, file=file.path(path, paste0("popsspe_", application), "static_avai",
                         paste(pop, "spe_full_avai_szgroup_nodes_semester1_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
 } else{
   if(quarter=="Q3"){
       dd      <- displace_dat[displace_dat$Quarter=="Q3",]
       dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
       print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
       ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
       library(doBy)
       ddd      <- orderBy(~pt_graph, ddd)
       dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
       dat_full <- ddd[, c("pt_graph", "avai")]
  } else{
       dd      <- displace_dat[displace_dat$Quarter=="Q4",]
       dd[,c(5:ncol(dd))] <- apply(dd[,c(5:ncol(dd))],2, NormAbu)
       print(apply(dd[,c(5:ncol(dd))], 2, sum))   # should return 1
       ddd     <- gather(dd, key=szgroup,value=avai,5:ncol(dd)) # reshape to long format
       library(doBy)
       ddd      <- orderBy(~pt_graph, ddd)
       dat      <- ddd[ddd$szgroup %in% paste0('SG',the_selected_szgroups), c("pt_graph", "avai")]
       dat_full <- ddd[, c("pt_graph", "avai")]
  } 
   write.table(dat, file=file.path(path, paste0("popsspe_", application), "static_avai",
                         paste(pop, "spe_avai_szgroup_nodes_semester2_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
   write.table(dat_full, file=file.path(path, paste0("popsspe_", application), "static_avai",
                        paste(pop, "spe_full_avai_szgroup_nodes_semester2_updated.dat", sep="")), quote=FALSE, col.names=TRUE, row.names=FALSE)
 }
 
 
 




} # end FALSE

cat("coupling to nbcp...done\n")

