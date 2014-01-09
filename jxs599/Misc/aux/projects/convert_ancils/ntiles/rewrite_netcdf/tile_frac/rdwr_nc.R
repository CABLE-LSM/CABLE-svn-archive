#these files could be altered slightly for different inputs, definitely write headers

#################################################################
read_data <- function( infile ) {

   #--- open netcdf file      
   Linnc <- open.ncdf( infile, readunlim=FALSE ) 

   #--- these are vars in Linnc -> dimensions of main vars 
   Llon <- get.var.ncdf(Linnc,'longitude')
   Llat <- get.var.ncdf(Linnc,'latitude')
   Lt <- get.var.ncdf(Linnc,'t')

   #--- these are the main vars in the file to get
   Lpseudo <- get.var.ncdf(Linnc,'pseudo')
   Lfield <- get.var.ncdf(Linnc,'field1391')

	close.ncdf(Linnc)
   #########################################################
   ###               return "list" to plot.R             ###
   #########################################################
   list( innc=Linnc, lon=Llon, lat=Llat, t=Lt, pseudo=Lpseudo, field=Lfield ) 
}

#################################################################
#################################################################

write_data <- function( rdata, ndata ) {
   #est up format of main var(s) incl. in nc file
   #dimi = dim.def.ncdf( name, units, var, unlim=FALSE, create_dimvar=TRUE )
   x <- (rdata$lon)
   dim1 = dim.def.ncdf( 'longitude','degrees_east', x, unlim=FALSE, create_dimvar=TRUE )
   print(paste( 'lon ', length(rdata$lon) ) ) 
   y <- rdata$lat
   dim2 = dim.def.ncdf( 'latitude','degrees_north', y, unlim=FALSE, create_dimvar=TRUE )
   print( paste( 'lat ', length(y) ) )
   s <- ndata$newpseudo
   dim3 = dim.def.ncdf( 'pseudo','level', s, unlim=FALSE, create_dimvar=TRUE )
   print( paste( 'surface ', length(s) ) )
   t <- rdata$t
   dim4 = dim.def.ncdf( 't','days since 0000-01-01 00:00:00', t, unlim=FALSE, create_dimvar=TRUE )
   print( paste( 't  ', length(t) ))
   #these are the main vars   
   vfield <- var.def.ncdf( 'field1391',' ', list(dim1,dim2,dim3,dim4),missval=2.0e20, 
                        longname='Stash Code = 313', prec='single' )

   Loutnc <- create.ncdf( outfile, vfield  ) 
   
   #attributes not previously incl.
   att.put.ncdf( Loutnc, 'longitude', attname='point_spacing', attval= 'even', prec='char'  )
   att.put.ncdf( Loutnc, 'longitude', attname='modulo', attval= ' ',prec='char'  )

   att.put.ncdf( Loutnc, 'latitude', attname='point_spacing', attval= 'even', prec='char'  )

   att.put.ncdf( Loutnc, 'pseudo', attname='positive', attval= 'up', prec='char'  )

   att.put.ncdf( Loutnc, 't', attname='calendar', attval= '360_day', prec='char'  )
   att.put.ncdf( Loutnc, 't', attname='time_origin', attval= '01-JAN-0000:00:00:00', prec='char'  )

   att.put.ncdf( Loutnc, 'field1391', attname='source', attval= 'Unified Model Output (Vn 6.7):', prec='char'  )
   att.put.ncdf( Loutnc, 'field1391', attname='title', attval= 'Stash Code = 313',prec='char'  )
   att.put.ncdf( Loutnc, 'field1391', attname='date', attval= '01/01/00',prec='char'  )
   att.put.ncdf( Loutnc, 'field1391', attname='time', attval= '00:00', prec='char'  )
   att.put.ncdf( Loutnc, 'field1391', attname='long_name', attval= 'Stash Code = 313',prec='char'  )
   att.put.ncdf( Loutnc, 'field1391', attname='missing_value', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'field1391', attname='_FillValue', attval= 2.e20, prec='single'  )
   att.put.ncdf( Loutnc, 'field1391', attname='valid_min', attval= -1.850372e-18, prec='single'  )
   att.put.ncdf( Loutnc, 'field1391', attname='valid_max', attval= 416.6667, prec='single'  )

   ##global attributes
   att.put.ncdf( Loutnc, 0, attname='history', attval= 'Mon Oct  12 09:57:127EST 2009 - XCONV V1.91 Development', prec='char'  )

   put.var.ncdf(Loutnc,'field1391', vals=ndata$newfield )

	close.ncdf(Loutnc)
}


