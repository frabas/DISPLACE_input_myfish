
 # CALLED FROM DISPLACE when dyn_pop_sce.option(Options::avai_updater_on)

 args <- commandArgs(trailingOnly = TRUE)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    print("No arguments supplied to input2AvaiUpdater. Take pop 0")
 pop <-0
 }else{
   pop <- args[1]
}



 general <- list()

 general$application           <- "myfish"

 if(.Platform$OS.type == "unix") {
  general$main_path         <- file.path("~","ibm_vessels",paste("DISPLACE_input_", general$application, sep=''))
  }
 if(.Platform$OS.type == "windows") {
    general$main_path         <- file.path("C:","Users","fbas","Documents","GitHub", paste("DISPLACE_input_", general$application, sep=''))
  }

 obj <- data.frame(Survey=NA, Year=c(2014,2014), ShootLon=c(18.79,11.0450), ShootLat=c(55.8667,54.0875),
                      Stock="COD.2224", StockId=10, nb_indiv.0=100, nb_indiv.1=10, nb_indiv.2=1,
                      nb_indiv.3=1, nb_indiv.4=1, nb_indiv.5=1, nb_indiv.6=1, nb_indiv.7=1, nb_indiv.8=1, nb_indiv.9=1,
                      nb_indiv.10=1, nb_indiv.11=1, nb_indiv.12=1, nb_indiv.13=1)  # caution: THE FORMAT IS NOT FLEXIBLE AT ALL

 
 # ...or a more realistic input file:
 if(TRUE) obj <- read.table(file.path(general$main_path, paste("popsspe_", general$application, sep=''), 
                                  "static_avai", "input_file_for_displace_merger.csv"), sep=";", header=TRUE)
 
 obj <- obj[,c("Survey", "Year", "ShootLon", "ShootLat", "Stock", "StockId", paste("nb_indiv.", 0:13, sep=""))]   # caution: THE FORMAT IS NOT FLEXIBLE AT ALL
 
 # make lighter data for the testing purpose....
 obj <-  obj[obj$StockId %in% pop,]
 
 if(nrow(obj)==0) stop(paste("This StockId", pop, "is not in the displace_input_for_data_merger.csv file!!"))

 write.table(obj, file= file.path( general$main_path, paste("popsspe_", general$application, sep=''), "static_avai",
             "displace_input_for_data_merger.dat"), sep=" ", row.names=FALSE, quote=FALSE)   ## CAUTION file is .dat and sep is white space




