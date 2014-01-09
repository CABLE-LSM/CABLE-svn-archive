
read_data <- function( obsfile, cablefile, oldcable, flux ) {
   cnc <- open.ncdf( cablefile, readunlim=FALSE ) # open CABLE file
   onc <- open.ncdf( obsfile, readunlim=FALSE ) # open observed data file
   cncB <- open.ncdf( oldcable, readunlim=FALSE ) # open observed data file
   fLcdata <- read_cabledata( cnc, flux ) 
   if( oldFLAG) fLcdataB <- read_cabledata( cncB, flux ) 
   else fLcdataB =0
   fdt <- set_timing( cnc )
   fLodata <- read_obsdata( onc, flux, fLcdata[[5]] )
	close.ncdf(cnc)
	close.ncdf(onc)
	close.ncdf(cncB)
   list( Lcdata=fLcdata, ndt=fdt, Lodata=fLodata, LcdataB=fLcdataB ) 
}

#### read model data for each flux var. (NEE, Qle, Qh, W/m2) ###
read_cabledata<-function( cnc, flux ) {
   for( i in 1:length(flux) ) {
      if(i==1) cNEE = get.var.ncdf( cnc, flux[i] ) # read CABLE output data
	   if(i==2) cQle = get.var.ncdf( cnc, flux[i] ) # read CABLE output data
	   if(i==3) cQh = get.var.ncdf( cnc, flux[i] ) # read CABLE output data
	   if(i==4) cRnet = get.var.ncdf( cnc, flux[i] ) # read CABLE output data
   }
   clen = length(cQle) # num timesteps == for all flux vars
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet, len=clen ) 
}

### Get timing details  ????  86400  ???? ###
set_timing = function( cnc ) {
	time12 = get.var.ncdf( cnc, 'time', start=1, count=2 ) # read CABLE 'time' variable
	timestepsize = time12[2]-time12[1]
	tstepinday = 86400/timestepsize # number of time steps in a day#86400 secs/day
}

#### read obs data for each flux var. (NEE, Qle, Qh, W/m2) ###
read_obsdata = function( onc, flux, fclen ) {
   if( is.character(onc$var$Rnet$name) )  fcheck = TRUE  
   else  fcheck = FALSE 
   for( i in 1:length(flux) ) {
      if(i==1) oNEE <- get.var.ncdf(onc,flux[i] )   # read observed data
      if(i==2) oQle <- get.var.ncdf(onc,flux[i] )   # read observed data
      if(i==3) oQh <- get.var.ncdf(onc,flux[i] )   # read observed data
      if (fcheck) {
         if(i==4) oRnet <- get.var.ncdf(onc,flux[i] )   # read observed data
      }
      else {
         if(i==4) oRnet <- matrix(0, fclen ) 
      }
   }
   for(i in 1:fclen){ ### this is bc Tumbarumba data has gaps
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

# reshape data into day per column
cdaily_data <- function( rdata ) { 
   cNEE=matrix(rdata$Lcdata$NEE, ncol = rdata$ndt,byrow=TRUE) 
   cQle=matrix(rdata$Lcdata$Qle, ncol = rdata$ndt,byrow=TRUE) 
   cQh=matrix(rdata$Lcdata$Qh, ncol = rdata$ndt,byrow=TRUE) 
   cRnet=matrix(rdata$Lcdata$Rnet, ncol = rdata$ndt,byrow=TRUE) 
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
}
cdaily_dataB <- function( rdata ) { 
   cNEE=matrix(rdata$LcdataB$NEE, ncol = rdata$ndt,byrow=TRUE) 
   cQle=matrix(rdata$LcdataB$Qle, ncol = rdata$ndt,byrow=TRUE) 
   cQh=matrix(rdata$LcdataB$Qh, ncol = rdata$ndt,byrow=TRUE) 
   cRnet=matrix(rdata$LcdataB$Rnet, ncol = rdata$ndt,byrow=TRUE) 
   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
}
odaily_data <- function( rdata ) { 
   oNEE = matrix( rdata$Lodata$NEE, ncol=rdata$ndt, byrow=TRUE )
   oQle = matrix(  rdata$Lodata$Qle, ncol=rdata$ndt, byrow=TRUE )
   oQh = matrix( rdata$Lodata$Qh, ncol=rdata$ndt, byrow=TRUE )
   oRnet = matrix( rdata$Lodata$Rnet, ncol=rdata$ndt, byrow=TRUE )
   list( NEE=oNEE, Qle=oQle, Qh=oQh, Rnet=oRnet )
}

