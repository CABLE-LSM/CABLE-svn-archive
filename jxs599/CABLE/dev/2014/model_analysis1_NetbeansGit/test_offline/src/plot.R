#####################################################################
### code starts here. main loop over the sites listed in main.nml ###      
### reads data from cable and obs. also converts to daily column  ###
### format for use in later functions. then for each site loops   ###
### over plot type, and comp. reduced chi squared vals for cases  ###
### within each plot type (diurnal and season only). can compare  ###
### older version of cable if given (diurnal only). finally       ###
### collates chi data for overall performance review at a glance  ###            
#####################################################################

#####################################################################
###           Load aux. functions and namelist:                   ###      
#####################################################################
source('../main.nml')
source('read_funcs.R')
source('support_funcs.R')
source('diurnal.R')
source('seasonal.R')
source('timeseries.R')
source('window.R')

#####################################################################
###           Load R libraries to use here                        ###      
#####################################################################
library(ncdf) # load lib to interpret netcdf format
library(Hmisc) # load lib to plot error bars 

#####################################################################
###  declare Goodness Of Fit arrays comp. per site and init. to 0 ###
###  1st dim. is over sites (i below in call), 2nd  dim. is over  ###
### flux types, 3rd dim. is over no. of plots per flux type * no. ### 
###  of models(2 here)                                            ###
#####################################################################
di_plotsperflux =8 #num. of seasons *2 models
digof = array( 0, c( length(sitenames), length(flux), di_plotsperflux ) ) 
seas_plotsperflux =2 #num. of models
seasgof = array( 0, c( length(sitenames), length(flux), seas_plotsperflux ) ) #1=annual * 2  models


#####################################################################
### main loop over the sites listed in main.nml                   ###      
#####################################################################
for(i in 1:length(sitenames) ) { # plots for each site listed in  main.nml
   setwd( paste('../out/',sitenames[i],sep='') )
   #fobs <- paste('../../',obsfiles[i],sep='') 
   fobs <- paste('',obsfiles[i],sep='') 
   ###############################################
   ###  reads all data and comp. daily format  ###
   ###############################################
   rdata <- read_data( fobs, cablefile, olderversion, flux )
   cday <- cdaily_data( rdata ) #puts data into columns=day
   oday <- odaily_data( rdata ) #puts data into columns=day
   if( oldFLAG)   cdayB <- cdaily_dataB( rdata ) #puts data into columns=day

   ###############################################
   ###  loops over plot types spec. in main.nml###
   ###############################################
   for( ii in 1:length(plottypes) ) {
      #####################################     
      ###    mean diurnal fluxes       ####     
      #####################################     
      if( plottypes[ii]=='diurnal')
         digof[i,,] <- diurnal( sitenames[i], rdata, cday, oday,cdayB, outtype, flux )
      #####################################     
      ###    mean monthly fluxes       ####     
      #####################################     
      if( plottypes[ii]=='seasonal')
         seasgof[i,,] <- seasonal( sitenames[i], rdata, cday, oday, cdayB, outtype, flux )
      #####################################     
      ###     fluxes windowed on        ###
      ###      different scales         ###     
      #####################################     
      if( plottypes[ii]=='window')
         window( sitenames[i], rdata, outtype, length(flux) ) 
      #####################################     
      ###  time series analysis see     ###
      ### timeseries.R re: args list    ###
      #####################################     
      if( plottypes[ii]=='timeseries')
         timeseries(sitenames[i],cablefile,outtype,-1,-1) 
   }
   rm ( rdata, cday, oday )
   setwd('../')
}

