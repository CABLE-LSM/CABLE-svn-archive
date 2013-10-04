#########################################################################
###  in this file of functions which are generic to most of the plot  ###
###  type - and hence done only once where possible - we read the     ###
###  data from cable outputs (new and old) and convert them into      ###
###  workable formats for the rest of the code. the major difference  ###
###  b/n this code and earlier version requirements is that now we    ###
###  account for patches of different vegetation etc. in cable data   ###                       
#########################################################################

#################################################################
### this is main read function called from plot.R to org.     ###
### read of cable (conditionally old version as well -        ###
### - oldFLAG) and obs data from netcdf binaries.             ###
#################################################################
read_data <- function( obsfile, cablefile, oldcable, flux ) {
   #########################################################
   ###           open netcdf files                       ###
   #########################################################
   cnc <- open.ncdf( cablefile, readunlim=FALSE ) 
   onc <- open.ncdf( obsfile, readunlim=FALSE ) 
   cncB <- open.ncdf( oldcable, readunlim=FALSE ) 
   #########################################################
   ### put cable data info into list deref. by $ syntax  ###    
   #########################################################
   fLcdata <- read_cabledata( cnc, flux ) 
   #########################################################
   ### read old cable data given FLAG=TRUE in main.nml   ###    
   #########################################################
   if( oldFLAG) fLcdataB <- read_cabledata( cncB, flux ) 
   else fLcdataB =0
   #########################################################
   ### determine num. of timesteps in a day  from data   ###    
   #########################################################
   fdt <- set_timing( cnc )
   #########################################################
   ### read obs data, send data length for loop therein  ###
   #########################################################
   fLodata <- read_obsdata( onc, flux, fLcdata[[9]] )
   #########################################################
   ###          close netcdf files                       ###
   #########################################################
	close.ncdf(cnc); 	close.ncdf(onc); close.ncdf(cncB)
   #########################################################
   ###               return "list" to plot.R             ###
   #########################################################
   list( Lcdata=fLcdata, ndt=fdt, Lodata=fLodata, LcdataB=fLcdataB ) 
}

#################################################################
### read cable data for given file - can be older version as  ###
### well. reads data for each flux var. (NEE, Qle, Qh, W/m2)  ###
### from netcdf file into 1D vector, which is then converted  ###
### to reflect patches (cXXX), and also total fluxes from the ###
### site over all patches (tXXX).                             ###
#################################################################
read_cabledata<-function( cnc, flux ) {
   #########################################################
   ### Get number of patches and patch fractions         ###
   #########################################################
   cpatchfrac = c()
   cpatchfrac <- get.var.ncdf(cnc,'patchfrac')
   cnpatch=cnc$dim$patch[[2]]
   #########################################################
   ### read cable data for each flux type into 1D vector ###    
   #########################################################
   for( i in 1:length(flux) ) {
      if(i==1) fcNEE = get.var.ncdf( cnc, flux[i] ) 
	   if(i==2) fcQle = get.var.ncdf( cnc, flux[i] ) 
	   if(i==3) fcQh = get.var.ncdf( cnc, flux[i] ) 
	   if(i==4) fcRnet = get.var.ncdf( cnc, flux[i] ) 
   }
#   cpatchfrac[1] = 0.5
#   cpatchfrac[2] = 0.5
   #########################################################
   ### determine num. entries/timesteps in cable data    ###
   #########################################################
   clen = length(fcQle) / cnpatch  
   #########################################################
   ### convert 1D data into 2D array by  patch type      ###    
   #########################################################
   cNEE = matrix(fcNEE,ncol=cnpatch,byrow=TRUE)
   cQle = matrix(fcQle,ncol=cnpatch,byrow=TRUE)
   cQh = matrix(fcQh,ncol=cnpatch,byrow=TRUE)
   cRnet = matrix(fcRnet,ncol=cnpatch,byrow=TRUE)
   #########################################################
   ### sum total flux contributed by whole site          ###    
   #########################################################
   tNEE = array( 0, c(clen) )
   tQle = array( 0, c(clen) )
   tQh = array( 0, c(clen) )
   tRnet = array( 0, c(clen) )
   for(M in 1:cnpatch){
      tNEE = tNEE  + (cNEE[,M] *cpatchfrac[M] )
      tQle = tQle  + (cQle[,M] *cpatchfrac[M] )
      tQh  = tQh   + (cQh[,M] *cpatchfrac[M]  )
      tRnet= tRnet + (cRnet[,M] *cpatchfrac[M]  )
   }
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet, tNEE=tNEE, tQle=tQle, tQh=tQh, tRnet=tRnet,
         len=clen, patchfrac=cpatchfrac, npatch=cnpatch ) 
}

#################################################################
###  determine number of timesteps per day in cable data      ###  
#################################################################
set_timing = function( cnc ) {
   #########################################################
   ### read time var.                                    ###    
   #########################################################
	time12 = get.var.ncdf( cnc, 'time', start=1, count=2 ) 
	timestepsize = time12[2]-time12[1]
   #########################################################
	### num. time steps/day where there is 86400 secs/day ###
   #########################################################
	tstepinday = 86400/timestepsize 
}

#################################################################
### read   obs data for given file for each flux var. (NEE,   ###
### Qle, Qh, W/m2), from netcdf file into 1D vector, which is ###
### LATER converted to daily-scale fluxes. correction is made ###
### missiing or bad data in some files.                       ###
#################################################################
read_obsdata = function( onc, flux, fclen ) {
   if( is.character(onc$var$Rnet$name) )  fcheck = TRUE  
   else  fcheck = FALSE 
   for( i in 1:length(flux) ) {
      if(i==1) oNEE <- get.var.ncdf(onc,flux[i] )
      if(i==2) oQle <- get.var.ncdf(onc,flux[i] )
      if(i==3) oQh <- get.var.ncdf(onc,flux[i] )
      if (fcheck) {
         if(i==4) oRnet <- get.var.ncdf(onc,flux[i] ) 
      }
      else {
         if(i==4) oRnet <- matrix(0, fclen ) 
      }
   }
   for(i in 1:fclen){ 
      if (oNEE[i]=='NA') oNEE[i]=0.0
      if (oNEE[i] < -9000.0) oNEE[i]=0.0
      if (oQle[i] < -9000.0) oQle[i]=0.0
      if (oQh[i] < -9000.0) oQh[i]=0.0
      if (fcheck) {
         if (oRnet[i] < -9000.0) oRnet[i]=0.0
      }
   }
   list( NEE=oNEE, Qle=oQle, Qh=oQh, Rnet=oRnet, check=fcheck, onc=onc) 
}

#################################################################
### reshape "total" flux data into day per column format      ###  
#################################################################
cdaily_data <- function( rdata ) { 
   cNEE=matrix(rdata$Lcdata$tNEE, ncol = rdata$ndt,byrow=TRUE) 
   cQle=matrix(rdata$Lcdata$tQle, ncol = rdata$ndt,byrow=TRUE) 
   cQh=matrix(rdata$Lcdata$tQh, ncol = rdata$ndt,byrow=TRUE) 
   cRnet=matrix(rdata$Lcdata$tRnet, ncol = rdata$ndt,byrow=TRUE) 
   cndays = length(cNEE[,1])
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet,ndays=cndays )
}

#################################################################
### reshape "patch" flux data into day per column format      ###  
#################################################################
cdaily_data_patches <- function( rdata, ndays ) {
   cNEE=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
   cQle=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
   cQh=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
   cRnet=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
   #########################################################
   ### loop over number of patches, then reshape into    ###
   ### day per column format as before                   ### 
   #########################################################
   for( i in 1:rdata$Lcdata$npatch ) { 
      for( j in 1:ndays ) { 
         for( k in 1:rdata$ndt ) { 
            kk =  ( (j-1)*rdata$ndt ) + k 
            cNEE[i,j,k] = rdata$Lcdata$NEE[kk,i]
            cQle[i,j,k] = rdata$Lcdata$Qle[kk,i]
            cQh[i,j,k] = rdata$Lcdata$Qh[kk,i]
            cRnet[i,j,k] = rdata$Lcdata$Rnet[kk,i]
         }
      }
   }
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
}

#################################################################
### reshape "total" flux data into day per column format      ### 
### for older version of cable data - if so desired           ###   
#################################################################
cdaily_dataB <- function( rdata ) { 
   cNEE=matrix(rdata$LcdataB$tNEE, ncol = rdata$ndt,byrow=TRUE) 
   cQle=matrix(rdata$LcdataB$tQle, ncol = rdata$ndt,byrow=TRUE) 
   cQh=matrix(rdata$LcdataB$tQh, ncol = rdata$ndt,byrow=TRUE) 
   cRnet=matrix(rdata$LcdataB$tRnet, ncol = rdata$ndt,byrow=TRUE) 
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
}

#################################################################
### reshape "patch" flux data into day per column format      ###  
### for older version of cable data - if so desired           ###   
#################################################################
cdaily_dataB_patches <- function( rdata, ndays ) {
   cNEE=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
   cQle=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
   cQh=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
   cRnet=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
   #########################################################
   ### loop over number of patches, then reshape into    ###
   ### day per column format as before                   ### 
   #########################################################
   for( i in 1:rdata$LcdataB$npatch ) { 
      for( j in 1:ndays ) { 
         for( k in 1:rdata$ndt ) { 
            kk =  ( (j-1)*rdata$ndt ) + k 
            cNEE[i,j,k] = rdata$LcdataB$NEE[kk,i]
            cQle[i,j,k] = rdata$LcdataB$Qle[kk,i]
            cQh[i,j,k] = rdata$LcdataB$Qh[kk,i]
            cRnet[i,j,k] = rdata$LcdataB$Rnet[kk,i]
         }
      }
   }
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
}

#################################################################
### reshape " obs " flux data into day per column format      ###  
#################################################################
odaily_data <- function( rdata ) { 
   oNEE = matrix( rdata$Lodata$NEE, ncol=rdata$ndt, byrow=TRUE )
   oQle = matrix(  rdata$Lodata$Qle, ncol=rdata$ndt, byrow=TRUE )
   oQh = matrix( rdata$Lodata$Qh, ncol=rdata$ndt, byrow=TRUE )
   oRnet = matrix( rdata$Lodata$Rnet, ncol=rdata$ndt, byrow=TRUE )
   list( NEE=oNEE, Qle=oQle, Qh=oQh, Rnet=oRnet )
}

