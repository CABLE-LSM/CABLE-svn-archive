#####################################################################
infile <- 'soil_snow.nc'
outfile <- 'newsoil_snow.nc'
source('proc_nc.R')

#####################################################################
###           Load R libraries to use here                        ###      
#####################################################################
library(ncdf) # load lib to interpret netcdf format

#####################################################################
   fin <- paste(infile,sep='') 
   ###############################################
   ###  reads all data and comp.               ###
   ###############################################
   rdata <- read_data( fin )
   ndata <- comp_data( rdata )
   ndata <- write_data( rdata, ndata ) 
  
 

