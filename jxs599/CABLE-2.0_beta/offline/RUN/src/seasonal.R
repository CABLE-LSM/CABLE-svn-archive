##############################################################################################
###      seasonal ROUTINE seen here produces plots of mean fluxes as a function of month   ###
###      (over the period of a year). the mean flux per month is calc. from the mean daily ###
###      flux, which is computed from the flux/timestep data in the "met" files (obs data) ###
###      and the output from cable (model data), and passed to this function in day/column ###
###      format (_days). note however that monthly averages are computed over # of years   ###
###      present in the data and the mean monthly flux displayed represents the mean over  ###
###      these years. for the sake of comparison between model and obs. data we include in ###
###      each plot a reduced chi squared value, which is computed by considering the       ###
###      statistical variance in the monthlyi obs. flux over the number of years included  ###
###      in the data.         - JS april 2009                                              ###
##############################################################################################

##############################################################################################
### first function called from main program to make seasonal plot - defines output type    ###
### of plot and loops over flux type to appear in each plot.                               ###  
##############################################################################################
seasonal <- function( sitename, rdata, cday, oday, cdayB, outtype, flux ){
   gof = array( 0,c( length(flux), 2 ) ) #goodness of fit 2 <= old model aswell  
	outfilename = paste(sitename, 'Seasonal.', outtype, sep='') # name of plot if not to screen
   outfile_func( outfilename, outtype, sitename ) 
	layout(matrix(1:length(flux),2,2))
   ndt <- rdata$ndt 
	for(j in 1:length(flux)) { 
      if(j==1) gof[j,]<-seas_func(j,cday$NEE,oday$NEE,cdayB$NEE,ndt,flux,sitename,rdata$Lodata$onc) 
      if(j==2) gof[j,]<-seas_func(j,cday$Qle,oday$Qle,cdayB$Qle,ndt,flux,sitename,rdata$Lodata$onc)
      if(j==3) gof[j,]<-seas_func(j,cday$Qh,oday$Qh,cdayB$Qh, ndt,flux,sitename,rdata$Lodata$onc)
      if(j==4) gof[j,]<-seas_func(j,cday$Rnet,oday$Rnet,cdayB$Rnet, ndt,flux,sitename,rdata$Lodata$onc )
	}		
	if(outtype!='screen') dev.off() 	#close graphics file if used:
   return( gof )
} 

##############################################################################################
### given a flux type to plot this function calls subrutines to - compute the mean monthly ###
### flux per year and also the mean monthly flux over all years for obs. and model data;   ###  
### - variance in the obs. data is computed, and then a reduced chi squared value. finally ###  
### a plotting function is called.                                                         ### 
##############################################################################################
seas_func <- function(fluxtype,cable_data,obs_data,old_data,ndt,flux,sitename,onc ) {
   nmonths_pa = 12     
	ndays=length(cable_data[,1]) #find n days in data set
	nyears=as.integer(ndays/365) #find n years in data set
	cable=c() #initialise monthly averages
   cable_i = array( 0, c( nyears, nmonths_pa ) ) #dec.array & init.to 0 monthly vals
	cableB =c() #initialise monthly averages
   cable_iB = array( 0, c( nyears, nmonths_pa ) ) #dec.array & init.to 0 monthly vals
	obs=c()   #
	gof=c()   #
   obs_i = array( 0, c( nyears, nmonths_pa ) ) #dec.array & init.to 0 monthly vals
	sigma=c() #
	plus_sigma=c() #initialise monthly averages
	minus_sigma=c() #initialise monthly averages
   ################################################################
   ####       get monthly flux values from cable data       #######
   ################################################################
   carrays <- compute_seasonal_arrays(ndays, nmonths_pa, nyears, cable_data ) 
   cable <- carrays$xmonthly
   cable_i <- carrays$xmonthly_i
   maxcable <- carrays$maxx
   mincable <- carrays$minx
   ################################################################
   ####       get monthly flux values from old_cable data   #######
   ################################################################
   if( oldFLAG ) {
      carraysB <- compute_seasonal_arrays(ndays, nmonths_pa, nyears, old_data ) 
      cableB <- carraysB$xmonthly
      cable_iB <- carraysB$xmonthly_i
      maxcableB <- carraysB$maxx
      mincableB <- carraysB$minx
   }
   ################################################################
   ####       check that obs data recorded Rnet flux        #######
   ################################################################
   check <- func_ochk( onc, flux[fluxtype] ) 

   ################################################################
   ####       get monthly flux values from obs data         #######
   ################################################################
   if (check) {
      oarrays <- compute_seasonal_arrays(ndays, nmonths_pa, nyears, obs_data )
      obs <- oarrays$xmonthly
      obs_i <- oarrays$xmonthly_i
      maxobs <- oarrays$maxx
      minobs <- oarrays$minx
   } 
   else {
       maxobs = 0.0
       minobs = 200.0
   }
   ################################################################
   ####      compute statisitical scatter in obs data       #######
   ################################################################
   fsigma <- compute_scatter( nmonths_pa, nyears, obs, obs_i ) 
   sigma <- fsigma$sigma
   plus_sigma <- fsigma$plus_sigma
   minus_sigma <- fsigma$minus_sigma
   maxscatter = max( plus_sigma )
   minscatter = min( minus_sigma )

   ################################################################
   ####       compute red. chi squared value                #######
   ################################################################
   gof[1] <- compute_chisq( nmonths_pa, check, cable, obs, sigma )
   if( oldFLAG )gof[2] <- compute_chisq( nmonths_pa, check, cableB, obs, sigma )
   else gof[2]=0.0 
  
   ################################################################
   ####       plot monthly flux data                        #######
   ################################################################
   plot_seasonal( fluxtype, nmonths_pa, cable,cableB, obs_i, obs, plus_sigma, minus_sigma, mincable,
      minobs,minscatter, maxcable,maxobs,maxscatter, check,sitename, flux, gof )
   return(gof)
} 

##############################################################################################
### generic function to compute monthly flux data, from computed daily avgs over number    ###
### of years, computes mean flux per month over all years and also mean flux per month per ###
### year so that scatter in sample can be computed.                                        ###
##############################################################################################
 compute_seasonal_arrays <- function( ndays, nmonths_pa, nyears, xdata ) {
   month_start=c() 	#days on which each month begin:
	month_start[1]=1; month_start[2]=32; month_start[3]=60 
	month_start[4]=91; month_start[5]=121;	month_start[6]=152 
	month_start[7]=182; month_start[8]=213; month_start[9]=244 
	month_start[10]=274; month_start[11]=305; month_start[12]=335
	month_start[13]=366 #i.e. beginning of next year
   xmonthly_i = array( 0, c( nyears, nmonths_pa ) ) #dec.array & init.to 0 monthly vals
	xmonthly=c() #initialise monthly averages
   avday = c() 
   avday[1:ndays] = 0
	for(i in 1:ndays) avday[i]=mean(xdata[i,]) # calc daily average flux
	for(l in 1:nmonths_pa){ #calc mothly averages 
		month_length=month_start[l+1]-month_start[l]
		xmonth=0 #initialise
		for(k in 1:nyears){ #sum daily avgs over months, for each year
			fxmonth <- fsumONEd( xmonth, avday, month_start[l], month_start[l+1]-1, k )
			xmonth  <- fxmonth$cummx
         xmonthly_i[k,l] = fxmonth$ittx / month_length
		}	
		xmonthly[l]=xmonth/(month_length*nyears) #divide by ndays summed over
   }
	maxx = max(xmonthly)
   minx = min(xmonthly)
   list( xmonthly=xmonthly, xmonthly_i=xmonthly_i, maxx=maxx, minxx=minx )
}


##############################################################################################
### plotting function of monthly flux data (obs. and model) and displays chi squared value ###
##############################################################################################
plot_seasonal<- function( fluxtype, n, cable,cableB,obs_i, obs, plus_sigma, minus_sigma, mincable,minobs,minscatter, maxcable,maxobs,maxscatter, check, sitename, flux, gof  ) {   
   require( graphics )
   stid=c(1,60,152,244,335) # seasonal divisions in days (over year)
	fnid=c(59,151,243,334,365)# seasonal divisions in days (over year)
   #labels=c('DJF','MAM','JJA','SON') #Dec,Jan,Feb=summer etc.
	units=c('umol/m2/s','W/m2','W/m2','W/m2') # flux units
	xloc=c(1:n) # set location of x-coords
   ftop =  max(maxcable,maxobs) + ( 0.2 * max(maxcable,maxobs) )
	# Plot CABLE output 
	plot(xloc, cable, type="l",xaxt="n",xlab='Month', 
		   ylab=paste('Average',flux[fluxtype],'flux',units[fluxtype]),lwd=2,col='salmon',
		   ylim=c( min(mincable,minobs), ftop ))
	# Plot obs 
	if (check ) {
      lines(xloc,obs,lwd=2)
      if( oldFLAG ) lines(xloc,cableB,lwd=2,col='green')
      axis(1,at=c(2,4,6,8,10,12),labels=c('2','4','6','8','10','12'))
  	   title(paste(sitename,flux[fluxtype])) # add title
   }
#temp. commented
#   temp1 <- expression(paste (chi[red]^{"2"}, '= ' ) )
#   temp1b = format( gof[1], digits=3 ) 
#   temp2 <- as.character( temp1b )
#   archx <- c(8, 9 )
#   ftop =  max(maxcable,maxobs) + ( 0.1 * max(maxcable,maxobs) )
#   archy <- c( ftop, ftop )
#   archa <- c(temp1, temp2 )
#   text(archx,archy, archa, adj=0.5 )  
   if( oldFLAG ) legend(1,ftop,c('CABLE','obs','old'), lty=1 ,col=c('salmon','black','green'),lwd=2,bty="n")
   else legend(1,ftop,c('CABLE','obs'), lty=1 ,col=c('salmon','black'),lwd=2,bty="n")
}

##############################################################################################
##############################################################################################







