
#################################################################
read_data <- function( infile ) {
   Linnc <- open.ncdf( infile, readunlim=FALSE ) 
   Lsoiltemp <- get.var.ncdf(Linnc,'soiltemp')
   #print( paste('soiltemp  ',soiltemp ) )
   #print( dim(soiltemp) )
   lon_vals <- get.var.ncdf(Linnc,'longitude')
   lat_vals <- get.var.ncdf(Linnc,'latitude')
   t_vals <- get.var.ncdf(Linnc,'t')
   surface_vals <-    c(0.1, 0.25, 0.65, 2.0, 3.0, 4.0) 
	close.ncdf(Linnc)
#   print(( Linnc$dim) )
   #########################################################
   ###               return "list" to plot.R             ###
   #########################################################
   list( soiltemp=Lsoiltemp, innc=Linnc, lon=lon_vals,lat=lat_vals,t=t_vals,surf=surface_vals ) 
}

comp_data <- function( rdata ) {
   dlen <- length(rdata$soiltemp[1,1,1,] )
   clen <- 6 
   blen <- length(rdata$soiltemp[1,,1,1] )
   alen <- length(rdata$soiltemp[,1,1,1] )
   newsoil <- array( 0, c(alen,blen,clen,dlen) )
   newsoil[,,1,] <- rdata$soiltemp[,,1,]
   newsoil[,,2,] <- rdata$soiltemp[,,2,]
   newsoil[,,3,] <- rdata$soiltemp[,,3,]
   newsoil[,,4,] <- rdata$soiltemp[,,4,]
   newsoil[,,5,] <- rdata$soiltemp[,,4,]
   newsoil[,,6,] <- rdata$soiltemp[,,4,]
#   print(newsoil[,,6,])
   list(  newst=newsoil ) 
}

write_data <- function( rdata, ndata ) {
   x <- rdata$lon
   print( length(rdata$lon) ) 
   dim1 = dim.def.ncdf( 'longitude','degrees_east', x, unlim=FALSE, create_dimvar=TRUE )
   y <- rdata$lat
   dim2 = dim.def.ncdf( 'latitude','degrees_north', y, unlim=FALSE, create_dimvar=TRUE )
   s <- rdata$surf
   dim3 = dim.def.ncdf( 'surface','level', s, unlim=FALSE, create_dimvar=TRUE )
   t <- rdata$t
   dim4 = dim.def.ncdf( 't','days since 0000-01-01 00:00:00', t, unlim=FALSE, create_dimvar=TRUE )
   varz <- var.def.ncdf( 'soiltemp','K', list(dim1,dim2,dim3,dim4),missval=2.0e20, longname='DEEP SOIL TEMP AFTER TIMESTEP', prec='single' )
   #varz <- var.def.ncdf( 'soiltemp', dim=list(dim1,dim2,dim3,dim4),missval=2.0e20 )
#   varz <- var.def.ncdf( dim=list(dim1,dim2,dim3,dim4) )
   Loutnc <- create.ncdf( outfile,varz ) 
   att.put.ncdf( Loutnc, 'longitude', attname='point_spacing', attval= 'even', prec='char'  )
   att.put.ncdf( Loutnc, 'longitude', attname='modulo', attval= ' ',prec='char'  )
   att.put.ncdf( Loutnc, 'latitude', attname='point_spacing', attval= 'even', prec='char'  )
   att.put.ncdf( Loutnc, 'surface', attname='positive', attval= 'up', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='modulo', attval= ' ', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='calendar', attval= '360_day', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='time_origin', attval= '01-JAN-0000:00:00:00', prec='char'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='Source', attval= 'Unified Model Output (Vn 6.0):', prec='char'  )
#   att.put.ncdf( Loutnc, 'soiltemp', attname='name', attval= 'soiltemp',prec='char'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='title', attval= 'DEEP SOIL TEMP AFTER TIMESTEP',prec='char'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='temp_date', attval= '01/01/00',prec='char'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='temp_time', attval= '00:00', prec='char'  )
#   att.put.ncdf( Loutnc, 'soiltemp', attname='long_name', attval= 'DEEP SOIL TEMP AFTER TIMESTEP',prec='char'  )
#   att.put.ncdf( Loutnc, 'soiltemp', attname='units', attval= 'K',prec='char'  )
#   att.put.ncdf( Loutnc, 'soiltemp', attname='missing_value', attval= 2.0e20, prec='single'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='_FillValue', attval= 2.0e20, prec='single'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='valid_min', attval= 193.9265, prec='single'  )
   att.put.ncdf( Loutnc, 'soiltemp', attname='valid_max', attval= 318.2573, prec='single'  )
   att.put.ncdf( Loutnc, 0, attname='history', attval= 'Thu Jul 23 11:18:42 EST 2009 - XCONV V1.91 Development', prec='char'  )
   put.var.ncdf(Loutnc, vals=ndata$newst)
	close.ncdf(Loutnc)
}



#
##################################################################
#### read cable data for given file - can be older version as  ###
#### well. reads data for each flux var. (NEE, Qle, Qh, W/m2)  ###
#### from netcdf file into 1D vector, which is then converted  ###
#### to reflect patches (cXXX), and also total fluxes from the ###
#### site over all patches (tXXX).                             ###
##################################################################
#read_data<-function( cnc, flux ) {
#   #########################################################
#   ### Get number of patches and patch fractions         ###
#   #########################################################
#   cpatchfrac = c()
#   cpatchfrac <- get.var.ncdf(cnc,'patchfrac')
#   cnpatch=cnc$dim$patch[[2]]
#   #########################################################
#   ### read cable data for each flux type into 1D vector ###    
#   #########################################################
#   for( i in 1:length(flux) ) {
#      if(i==1) fcNEE = get.var.ncdf( cnc, flux[i] ) 
#	   if(i==2) fcQle = get.var.ncdf( cnc, flux[i] ) 
#	   if(i==3) fcQh = get.var.ncdf( cnc, flux[i] ) 
#	   if(i==4) fcRnet = get.var.ncdf( cnc, flux[i] ) 
#   }
##   cpatchfrac[1] = 0.5
##   cpatchfrac[2] = 0.5
#   #########################################################
#   ### determine num. entries/timesteps in cable data    ###
#   #########################################################
#   clen = length(fcQle) / cnpatch  
#   #########################################################
#   ### convert 1D data into 2D array by  patch type      ###    
#   #########################################################
#   cNEE = matrix(fcNEE,ncol=cnpatch,byrow=TRUE)
#   cQle = matrix(fcQle,ncol=cnpatch,byrow=TRUE)
#   cQh = matrix(fcQh,ncol=cnpatch,byrow=TRUE)
#   cRnet = matrix(fcRnet,ncol=cnpatch,byrow=TRUE)
#   #########################################################
#   ### sum total flux contributed by whole site          ###    
#   #########################################################
#   tNEE = array( 0, c(clen) )
#   tQle = array( 0, c(clen) )
#   tQh = array( 0, c(clen) )
#   tRnet = array( 0, c(clen) )
#   for(M in 1:cnpatch){
#      tNEE = tNEE  + (cNEE[,M] *cpatchfrac[M] )
#      tQle = tQle  + (cQle[,M] *cpatchfrac[M] )
#      tQh  = tQh   + (cQh[,M] *cpatchfrac[M]  )
#      tRnet= tRnet + (cRnet[,M] *cpatchfrac[M]  )
#   }
#   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet, tNEE=tNEE, tQle=tQle, tQh=tQh, tRnet=tRnet,
#         len=clen, patchfrac=cpatchfrac, npatch=cnpatch ) 
#}
#
##################################################################
#### reshape "total" flux data into day per column format      ###  
##################################################################
#cdaily_data <- function( rdata ) { 
#   cNEE=matrix(rdata$Lcdata$tNEE, ncol = rdata$ndt,byrow=TRUE) 
#   cQle=matrix(rdata$Lcdata$tQle, ncol = rdata$ndt,byrow=TRUE) 
#   cQh=matrix(rdata$Lcdata$tQh, ncol = rdata$ndt,byrow=TRUE) 
#   cRnet=matrix(rdata$Lcdata$tRnet, ncol = rdata$ndt,byrow=TRUE) 
#   cndays = length(cNEE[,1])
#   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet,ndays=cndays )
#}
#
##################################################################
#### reshape "patch" flux data into day per column format      ###  
##################################################################
#cdaily_data_patches <- function( rdata, ndays ) {
#   cNEE=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
#   cQle=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
#   cQh=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
#   cRnet=array( 0, c(rdata$Lcdata$npatch, ndays, rdata$ndt ) )
#   #########################################################
#   ### loop over number of patches, then reshape into    ###
#   ### day per column format as before                   ### 
#   #########################################################
#   for( i in 1:rdata$Lcdata$npatch ) { 
#      for( j in 1:ndays ) { 
#         for( k in 1:rdata$ndt ) { 
#            kk =  ( (j-1)*rdata$ndt ) + k 
#            cNEE[i,j,k] = rdata$Lcdata$NEE[kk,i]
#            cQle[i,j,k] = rdata$Lcdata$Qle[kk,i]
#            cQh[i,j,k] = rdata$Lcdata$Qh[kk,i]
#            cRnet[i,j,k] = rdata$Lcdata$Rnet[kk,i]
#         }
#      }
#   }
#   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
#}
#
##################################################################
#### reshape "total" flux data into day per column format      ### 
#### for older version of cable data - if so desired           ###   
##################################################################
#cdaily_dataB <- function( rdata ) { 
#   cNEE=matrix(rdata$LcdataB$tNEE, ncol = rdata$ndt,byrow=TRUE) 
#   cQle=matrix(rdata$LcdataB$tQle, ncol = rdata$ndt,byrow=TRUE) 
#   cQh=matrix(rdata$LcdataB$tQh, ncol = rdata$ndt,byrow=TRUE) 
#   cRnet=matrix(rdata$LcdataB$tRnet, ncol = rdata$ndt,byrow=TRUE) 
#   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
#}
#
##################################################################
#### reshape "patch" flux data into day per column format      ###  
#### for older version of cable data - if so desired           ###   
##################################################################
#cdaily_dataB_patches <- function( rdata, ndays ) {
#   cNEE=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
#   cQle=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
#   cQh=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
#   cRnet=array( 0, c(rdata$LcdataB$npatch, ndays, rdata$ndt ) )
#   #########################################################
#   ### loop over number of patches, then reshape into    ###
#   ### day per column format as before                   ### 
#   #########################################################
#   for( i in 1:rdata$LcdataB$npatch ) { 
#      for( j in 1:ndays ) { 
#         for( k in 1:rdata$ndt ) { 
#            kk =  ( (j-1)*rdata$ndt ) + k 
#            cNEE[i,j,k] = rdata$LcdataB$NEE[kk,i]
#            cQle[i,j,k] = rdata$LcdataB$Qle[kk,i]
#            cQh[i,j,k] = rdata$LcdataB$Qh[kk,i]
#            cRnet[i,j,k] = rdata$LcdataB$Rnet[kk,i]
#         }
#      }
#   }
#   list( NEE=cNEE, Qle=cQle, Qh=cQh, Rnet=cRnet )
#}
#
##################################################################
#### reshape " obs " flux data into day per column format      ###  
##################################################################
#odaily_data <- function( rdata ) { 
#   oNEE = matrix( rdata$Lodata$NEE, ncol=rdata$ndt, byrow=TRUE )
#   oQle = matrix(  rdata$Lodata$Qle, ncol=rdata$ndt, byrow=TRUE )
#   oQh = matrix( rdata$Lodata$Qh, ncol=rdata$ndt, byrow=TRUE )
#   oRnet = matrix( rdata$Lodata$Rnet, ncol=rdata$ndt, byrow=TRUE )
#   list( NEE=oNEE, Qle=oQle, Qh=oQh, Rnet=oRnet )
#}
#
