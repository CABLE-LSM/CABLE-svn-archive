##############################################################################################
###      diurnal  ROUTINE seen here produces plots of mean fluxes as a function of hour    ###
###      (over the period of a day). the mean flux per timestep is calc. from the          ### 
###      flux/timestep data in the "met" files (obs data) and the output from CABLE (model ###
###      data), which is passed to this function in day/column format (_days). note        ###
###      however that plotted averages are computed over the number of years present in    ###
###      the data and the mean flux displayed represents the mean over these years. for    ###
###      the sake of comparison between model and obs. data we include in each plot a      ###
###      reduced chi squared value, which is computed by considering the statistical       ###
###      variance in the obs. flux/timestep over the number of years included in the data. ###
##############################################################################################

##############################################################################################
### first function called from main program to make diurnal  plot - defines output type    ###
### of plot and loops over flux type to appear in each plot.                               ###  
##############################################################################################
diurnal <- function( sitename, rdata, cday, oday,cdayB, outtype, flux ){
   seasons=4
   gof = array( 0,c( length(flux), 2*seasons ) ) #goodness of fit - *2 <= old_cable aswell
	outfilename = paste(sitename, 'Diurnal.', outtype, sep='') # name of plot if not to screen
   outfile_func( outfilename, outtype, sitename ) 
   layout(matrix(1:16,4,4))
   ndt <- rdata$ndt
	for(j in 1:length(flux)) { 
      if(j==1) gof[j,]<-di_func(j,cday$NEE,oday$NEE,cdayB$NEE,ndt,flux,sitename,rdata$Lodata$onc,oday$Rnet) 
      if(j==2) gof[j,]<-di_func(j,cday$Qle,oday$Qle,cdayB$Qle,ndt,flux,sitename,rdata$Lodata$onc,oday$Rnet) 
      if(j==3) gof[j,]<-di_func(j,cday$Qh,oday$Qh,cdayB$Qh,ndt,flux,sitename,rdata$Lodata$onc,oday$Rnet) 
      if(j==4) gof[j,]<-di_func(j,cday$Rnet,oday$Rnet,cdayB$Rnet,ndt,flux,sitename,rdata$Lodata$onc,oday$Rnet)
	}		
	if(outtype!='screen') dev.off() 	#close graphics file if used:
   return( gof )
} 

##############################################################################################
### given a flux type to plot this function calls subroutines to compute the mean flux per ###
### timestep (over nyears as well) for obs. and model data; compute the variance in the    ### 
### obs. data, and then a reduced chi squared value. finally a plotting function is called.### 
##############################################################################################
di_func <- function(fluxtype,cable_data,obs_data,old_data,ndt,flux,sitename,onc, Rnet) {
   gof=c()  
	ndays=length(cable_data[,1]) #find n days in data set
	nyears=as.integer(ndays/365) #find n years in data set
	cable=c() #dec.
   cable_i = array( 0, c( nyears, ndt ) ) #dec.array & init.to 0 monthly vals
	cableB=c() #dec.
   cable_iB = array( 0, c( nyears, ndt ) ) #dec.array & init.to 0 monthly vals
	obs=c()   #dec.
   obs_i = array( 0, c( nyears, ndt ) ) #dec.array & init.to 0 monthly vals
	sigma=c() #dec.
	plus_sigma=c() #dec.
	minus_sigma=c() #dec.
	daylight=c() #dec.

   for( season in 1:4 ) {
      ################################################################
      ####      mean flux/timestep values from cable data      #######
      ################################################################
      carrays <- compute_diurnal_arrays(season, ndt, nyears, cable_data ) 
      cable <- carrays$xdt
      cable_i <- carrays$xdt_i
      maxcable <- carrays$maxx
      mincable <- carrays$minx
      ################################################################
      ####      mean flux/timestep values from old version data  #####
      ################################################################
      if( oldFLAG) {
         carraysB <- compute_diurnal_arrays(season, ndt, nyears, old_data ) 
         cableB <- carraysB$xdt
         cable_iB <- carraysB$xdt_i
         maxcableB <- carraysB$maxx
         mincableB <- carraysB$minx
      }
      ################################################################
      ####       check that obs data recorded Rnet flux        #######
      ################################################################
      check <- func_ochk( onc, flux[fluxtype] ) 
   
      ################################################################
      ####      mean flux/timestep values from   obs data      #######
      ################################################################
      if (check) {
         oarrays <- compute_diurnal_arrays(season, ndt, nyears, obs_data ) 
         obs <- oarrays$xdt
         obs_i <- oarrays$xdt_i
         maxobs <- oarrays$maxx
         minobs <- oarrays$minx
      } 
      else {
         maxobs = 0.0 ; minobs = 200.0
      } 
      ################################################################
      ####      compute statisitical scatter in obs data       #######
      ################################################################
      fsigma <- compute_scatter( ndt, nyears, obs, obs_i ) 
      sigma <- fsigma$sigma
      plus_sigma <- fsigma$plus_sigma
      minus_sigma <- fsigma$minus_sigma
      maxscatter = max( plus_sigma )
      minscatter = min( minus_sigma )

      ################################################################
      ####       compute red. chi squared value                #######
      ################################################################
      gof[season] <- compute_chisq( ndt, check, cable, obs, sigma )
      if( oldFLAG ) gof[season+4] <- compute_chisq( ndt, check, cableB, obs, sigma )
      else gof[season+4] =0.0 
  
      ################################################################
      ####       plot monthly flux data                        #######
      ################################################################
      plot_diurnal(season, fluxtype, ndt, cable, cableB,obs_i, obs, plus_sigma, minus_sigma, mincable,
         minobs,minscatter, maxcable,maxobs,maxscatter, check,sitename, flux, gof[season] )
   }
   return(gof)
} 

##############################################################################################
### generic function to compute approp. flux data, from computed daily avgs over number    ###
### of years, computes mean flux per month over all years and also mean flux per month per ###
### year so that scatter in sample can be computed.                                        ###
##############################################################################################
compute_diurnal_arrays <- function( s, ndt, nyears, xdata ) {
   stid=c(1,60,152,244,335) # seasonal divisions in a year
   fnid=c(59,151,243,334,365)# seasonal divisions in a year
   xdt_i = array( 0, c( nyears, ndt ) ) #dec.array & init.to 0 monthly vals
   ixdt_i = array( 0, c( nyears, ndt ) ) #dec.array & init.to 0 monthly vals
	xdt=c() #initialise monthly averages
   ss=s        
   srange = fnid[s]-stid[s] #= num. of days in season
   ixdt=c(); ixdt[1:ndt] = 0 #initialise
	for(k in 1:nyears){	
	   for(l in 1:ndt){
         fxdt <- fsumTWOd(ixdt[l],xdata,stid[s],fnid[s],l,k) # sum func 
         ixdt_i[k,l] =  fxdt$ittx
         xdt_i[k,l] =  fxdt$ittx/(fnid[s]-stid[s])
         ixdt[l] <- fxdt$cummx
		}
		if(s==1){ # i.e. DJF
		   for(l in 1:ndt ){
            fxdt <- fsumTWOd(ixdt[l],xdata,stid[s+4],fnid[s+4],l,k) # sum func 
            xdt_i[k,l] =( ixdt_i[k,l] +  fxdt$ittx) / 90 
            ixdt[l] <- fxdt$cummx
			}
		}
	}
	# Then find the average of these fluxes:
	for(l in 1:ndt ){
	   if(s==1){ # i.e. DJF
         xdt[l]=ixdt[l]/ 90 / nyears #divide by ndays summed over
		}else{
			xdt[l]=ixdt[l]/((fnid[s]-stid[s])*nyears)
		}
   }
	maxx = max(xdt)
   minx = min(xdt)
   list( xdt=xdt, xdt_i=xdt_i, maxx=maxx, minxx=minx )
}

##############################################################################################
### plotting function of diurnal flux data (obs. and model) and displays chi squared value ###
##############################################################################################
plot_diurnal<- function( s,fluxtype, n, cable, cableB,obs_i, obs, plus_sigma, minus_sigma, mincable,
               minobs, minscatter, maxcable, maxobs, maxscatter, check, sitename, flux, gof ) {   
   require( graphics )
   stid=c(1,60,152,244,335) # seasonal divisions in days (over year)
	fnid=c(59,151,243,334,365)# seasonal divisions in days (over year)
	labels=c('DJF','MAM','JJA','SON') #Dec,Jan,Feb=summer etc.
	units=c('umol/m2/s','W/m2','W/m2','W/m2') # flux units
	xloc=c(1:n-1) # set location of x-coords
   fmax = max( maxcable, maxobs, maxscatter )
   fmin = min( mincable, minobs, minscatter )
   ftop = fmax + ( 0.3 * abs( fmax-fmin ) )
	# Plot CABLE output result:
	plot( xloc, cable, type="l",xaxt="n",xlab='Hour', 
		   ylab=paste('Average',flux[fluxtype],'flux', units[fluxtype]),
         lwd=2, col='salmon',
		   ylim=c( min(mincable,minobs), ftop ) )
  # Then plot obs result:
	if (check ) {
      lines( xloc, obs, lwd=2 )#orig
      if( oldFLAG ) lines( xloc, cableB, lwd=2, col='green' )
      axis( 1, at= c(0,6*n/24,12*n/24,18*n/24,23*n/24),labels=c('0','6','12','18','23'))
      title( paste( sitename, labels[s], flux[fluxtype] ) ) # add title
   }
#temp. commented
#   temp1 <- expression( paste (chi[red]^{"2"}, '= ' ) )
#   temp1b = format( gof, digits=3 ) 
#   temp2 <- as.character( temp1b )
#   archx <- c( 8, 13 )
#   ftop = fmax + ( 0.2 * abs( fmax-fmin ) )
#   archy <- c( ftop, ftop )
#   archa <- c( temp1, temp2 )
#   text( archx, archy, archa, adj=0.5 )  

   if( fluxtype==1 & s==1 ) {
      if( oldFLAG) legend(1,ftop,c('CABLE','obs','old'), lty=1 ,col=c('salmon','black','green'),lwd=2,bty="n")
      else legend( ((xloc[n-1]-xloc[1])/2)-3.5 , ftop, c('CABLE','obs'), lty=1 ,col=c('salmon','black'),lwd=2,bty="n")
   }
}
##############################################################################################
##############################################################################################





