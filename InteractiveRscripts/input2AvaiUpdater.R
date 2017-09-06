
 # CALLED FROM DISPLACE when dyn_pop_sce.option(Options::avai_updater_on)

 general <- list()

 general$application           <- "myfish"

 if(.Platform$OS.type == "unix") {
  general$main_path         <- file.path("~","ibm_vessels",paste("DISPLACE_input_", general$application, sep=''))
  }
 if(.Platform$OS.type == "windows") {
    general$main_path         <- file.path("C:","Users","fbas","Documents","GitHub", paste("DISPLACE_input_", general$application, sep=''))
  }

 fake <- data.frame(Survey=NA, Year=c(2014,2014), ShootLon=c(18.79,11.0450), ShootLat=c(55.8667,54.0875),
                      Stock="COD.2224", StockId=10, nb_indiv.0=100, nb_indiv.1=10, nb_indiv.2=1,
                      nb_indiv.3=1, nb_indiv.4=1, nb_indiv.5=1, nb_indiv.6=1, nb_indiv.7=1, nb_indiv.8=1, nb_indiv.9=1,
                      nb_indiv.10=1, nb_indiv.11=1, nb_indiv.12=1, nb_indiv.13=1)  # caution: THE FORMAT IS NOT FLEXIBLE AT ALL

 write.table(fake, file= file.path( general$main_path, paste("popsspe_", general$application, sep=''), "static_avai",
             "displace_input_for_data_merger.dat"), sep=" ", row.names=FALSE, quote=FALSE)   ## CAUTION file is .dat and sep is white space


