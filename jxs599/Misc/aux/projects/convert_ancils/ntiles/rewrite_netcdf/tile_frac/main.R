#####################################################################
infile <- 'snow_depth.nc'
outfile <- 'new_depth.nc'
source('rdwr_nc.R')
source('comp.R')

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
  
 

