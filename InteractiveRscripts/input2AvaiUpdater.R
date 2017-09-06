
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
                      Stock="COD.2224", StockId=10, nb_indiv.0=100, nb_indiv.1=10, nb_indiv2=1)

 write.table(fake, file= file.path( general$main_path,
             "input_file_for_displace_merger.csv"), sep=";", row.names=FALSE, quote=FALSE)


