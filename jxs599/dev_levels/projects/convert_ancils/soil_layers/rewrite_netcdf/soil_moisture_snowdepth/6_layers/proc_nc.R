
#################################################################
read_data <- function( infile ) {
   #--- open netcdf file      
   Linnc <- open.ncdf( infile, readunlim=FALSE ) 
   #--- these are the main vars in the file to get
   Lsnowdepth <- get.var.ncdf(Linnc,'snowdepth')
   Lsnowdepth_1 <- get.var.ncdf(Linnc,'snowdepth_1')
   Lsm <- get.var.ncdf(Linnc,'sm')
   #print( paste('soiltemp  ',soiltemp ) )
   #print( paste( 'read: snowdepth(dim)   ', dim(Lsnowdepth) ) )
   #--- these are vars in Linnc -> dimensions of main vars 
   Llon <- get.var.ncdf(Linnc,'longitude')
   Llat <- get.var.ncdf(Linnc,'latitude')
   Lt <- get.var.ncdf(Linnc,'t')
   Lsurface <- get.var.ncdf(Linnc,'surface')
   #Lsurface_1 <-    c(0.1, 0.25, 0.65, 2.0 ) #this is in var srface_l in nc file 
#for 6 layers
   Lsurface_1 <-    c(0.1, 0.25, 0.65, 2.0, 3.0, 4.0 ) #this is in var srface_l in nc file 
	close.ncdf(Linnc)
   #   print(( Linnc$dim) )
   #########################################################
   ###               return "list" to plot.R             ###
   #########################################################
   list( snowdepth=Lsnowdepth, snowdepth_1=Lsnowdepth_1, sm = Lsm, innc=Linnc, lon=Llon, lat=Llat, 
         t=Lt, surface=Lsurface, surface_1=Lsurface_1 ) 
}

comp_data <- function( rdata ) {
   dlen <- length(rdata$sm[1,1,1,] )
   #clen <- 4
#for 6 layers
   clen <- 6 
   blen <- length(rdata$sm[1,,1,1] )
   alen <- length(rdata$sm[,1,1,1] )
   newsm <- array( 0, c(alen,blen,clen,dlen) )
   newsm[,,1,] <- rdata$sm[,,1,]
   newsm[,,2,] <- rdata$sm[,,2,]
   newsm[,,3,] <- rdata$sm[,,3,]
   newsm[,,4,] <- rdata$sm[,,4,]
#for 6 layers
   newsm[,,5,] <- rdata$sm[,,4,]
   newsm[,,6,] <- rdata$sm[,,4,]
   #print(newsoil[,,6,])
   list(  newsm=newsm ) 
}

write_data <- function( rdata, ndata ) {
   #est up format of main var(s) incl. in nc file
   #dimi = dim.def.ncdf( name, units, var, unlim=FALSE, create_dimvar=TRUE )
   x <- (rdata$lon)
   dim1 = dim.def.ncdf( 'longitude','degrees_east', x, unlim=FALSE, create_dimvar=TRUE )
   ##print(paste( 'lon ', length(rdata$lon) ) ) 

   y <- rdata$lat
   dim2 = dim.def.ncdf( 'latitude','degrees_north', y, unlim=FALSE, create_dimvar=TRUE )
   ##print( paste( 'lat ', length(y) ) )

   s <- rdata$surface 
   dim3 = dim.def.ncdf( 'surface','level', s, unlim=FALSE, create_dimvar=TRUE )
   ##print( paste( 'surface ', length(s) ) )
   
   t <- rdata$t
   dim4 = dim.def.ncdf( 't','days since 0000-01-01 00:00:00', t, unlim=FALSE, create_dimvar=TRUE )
   ##print( paste( 't  ', length(t) ))

   s1 <- rdata$surface_1
   dim5 = dim.def.ncdf( 'surface_1','level', s1, unlim=FALSE, create_dimvar=TRUE )

   #vlon <- var.def.ncdf( 'longitude','degrees_east', dim1, prec='single' )
   #vlat <- var.def.ncdf( 'latitude','degrees_north', dim2, prec='single' )
   #vt <- var.def.ncdf( 't','days since 0000-01-01 00:00:00', list(dim1), missval=2.0e20,prec='single' )
   ##print( paste( 'surface_1', length(s1) ) )

   #these are the main vars   
   vsnowdepth <- var.def.ncdf( 'snowdepth','kg m-2', list(dim1,dim2,dim3,dim4),missval=2.0e20, 
                        longname='SNOW AMOUNT OVER LAND AFT TSTP KG/M2', prec='single' )

   vsnowdepth_1 <- var.def.ncdf( 'snowdepth_1','kg m-2', list(dim1,dim2,dim3,dim4),missval=2.0e20, 
                        longname='SNOW EDGE AFTER TIMESTEP          **', prec='single' )

   vsm <- var.def.ncdf( 'sm','kg m-2', list(dim1,dim2,dim5,dim4) ,missval=2.0e20, 
                        longname='SOIL MOISTURE CONTENT IN A LAYER', prec='single' )

   Loutnc <- create.ncdf( outfile, list( vsnowdepth, vsnowdepth_1, vsm)  ) 
   
   #attributes not previously incl.
   att.put.ncdf( Loutnc, 'longitude', attname='point_spacing', attval= 'even', prec='char'  )
   att.put.ncdf( Loutnc, 'longitude', attname='modulo', attval= ' ',prec='char'  )

   att.put.ncdf( Loutnc, 'latitude', attname='point_spacing', attval= 'even', prec='char'  )

   att.put.ncdf( Loutnc, 'surface', attname='positive', attval= 'up', prec='char'  )

   att.put.ncdf( Loutnc, 't', attname='modulo', attval= ' ', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='calendar', attval= '360_day', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='time_origin', attval= '01-JAN-0000:00:00:00', prec='char'  )

   att.put.ncdf( Loutnc, 'snowdepth', attname='source', attval= 'Unified Model Output (Vn 6.0):', prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='title', attval= 'SNOW AMOUNT OVER LAND AFT TSTP KG/M2',prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='date', attval= '01/01/00',prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='time', attval= '00:00', prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='missing_value', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='_FillValue', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='valid_min', attval= 0., prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth', attname='valid_max', attval= 50000., prec='single'  )

   att.put.ncdf( Loutnc, 'snowdepth_1', attname='source', attval= 'Unified Model Output (Vn 6.0):', prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='title', attval= 'SNOW EDGE AFTER TIMESTEP          **',prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='date', attval= '01/01/00',prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='time', attval= '00:00', prec='char'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='missing_value', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='_FillValue', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='valid_min', attval= 0.1111111, prec='single'  )
   att.put.ncdf( Loutnc, 'snowdepth_1', attname='valid_max', attval= 0.8888889, prec='single'  )

   att.put.ncdf( Loutnc, 'surface_1', attname='positive', attval= 'up', prec='char'  )

   att.put.ncdf( Loutnc, 'sm', attname='source', attval= 'Unified Model Output (Vn 6.0):', prec='char'  )
   att.put.ncdf( Loutnc, 'sm', attname='title', attval= 'SOIL MOISTURE CONTENT IN A LAYER',prec='char'  )
   att.put.ncdf( Loutnc, 'sm', attname='date', attval= '01/01/00',prec='char'  )
   att.put.ncdf( Loutnc, 'sm', attname='time', attval= '00:00', prec='char'  )
   att.put.ncdf( Loutnc, 'sm', attname='missing_value', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'sm', attname='_FillValue', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'sm', attname='valid_min', attval= 0., prec='single'  )
   att.put.ncdf( Loutnc, 'sm', attname='valid_max', attval= 916.2976, prec='single'  )
##global attributes
   att.put.ncdf( Loutnc, 0, attname='history', attval= 'Mon Sep 28 14:53:11 EST 2009 - XCONV V1.91 Development', prec='char'  )

   put.var.ncdf(Loutnc,'snowdepth', vals=rdata$snowdepth)
   put.var.ncdf(Loutnc, 'snowdepth_1',vals=rdata$snowdepth_1)
   put.var.ncdf(Loutnc, 'sm',vals=ndata$newsm)

	close.ncdf(Loutnc)
}


