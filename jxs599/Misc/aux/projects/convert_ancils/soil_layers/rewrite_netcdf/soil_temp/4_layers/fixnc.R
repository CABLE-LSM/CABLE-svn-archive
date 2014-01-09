#####################################################################
infile <- 'qrclim.nc'
outfile <- 'newqrclim.nc'
source('proc_nc.R')

#####################################################################
###           Load R libraries to use here                        ###      
#####################################################################
library(ncdf) # load lib to interpret netcdf format

#####################################################################
   fin <- paste(infile,sep='') 
   ###############################################
   ###  reads all data and comp. daily formats ###
   ###############################################
   rdata <- read_data( fin )
   ndata <- comp_data( rdata )
   ndata <- write_data( rdata, ndata ) 
  
 

