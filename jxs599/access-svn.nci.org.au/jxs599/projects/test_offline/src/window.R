
window = function( sitename, rdata, outtype, fluxlength ) {
	library(boot) # load bootstrap library
	outfilename=paste(sitename,'Window.',outtype,sep='')
   outfile_func (outfilename, outtype, sitename ) 
	windowstepsize=2 # resolution of graph
	maxwindow = 30 # in days, the largest averaging window
	#onc=open.ncdf(obsfile,readunlim=FALSE) # open observed data file
   Lch <- chk_obs( rdata$Lodata$onc )
   flux= Lch$flux ; units= Lch$units ; layout= Lch$layout
	tstepinday=rdata$ndt
	numcalcs=maxwindow*tstepinday/windowstepsize
   l = 1
   if(Lch$Rnet) {
	   obs_data = rdata$Lodata$Rnet
		cable_data = rdata$Lcdata$Rnet
      window_func( cable_data, obs_data, sitename, windowstepsize, flux, units,
               layout, numcalcs, tstepinday,maxwindow, l )
   }
   for(l in 2:fluxlength){
		if(l==2) cable_data = rdata$Lcdata$Qle
   	if(l==2) obs_data = rdata$Lodata$Qle
		if(l==3) cable_data = rdata$Lcdata$Qh
   	if(l==3) obs_data = rdata$Lodata$Qh
		if(l==4) cable_data = rdata$Lcdata$NEE
		if(l==4) obs_data = rdata$Lodata$NEE
      window_func( cable_data, obs_data, sitename, windowstepsize, flux, units,
            layout, numcalcs, tstepinday,maxwindow, l )
	}
	if(outtype!='screen') dev.off() 	# close graphics file if used:
} # End function avwinflux
##################################################################################

window_func <-function( cable_data, obs_data, sitename, windowstepsize, flux, units,
                  layout, numcalcs, tstepinday, maxwindow, l ) {
	tsteps=length(cable_data) # total number of timesteps
  	mvals = c() # initialise gradient values
	mwidth =c() # initialise 95% confidence interval for regression
	rvals = c() # initialise correlation coefficient values
	rmse = c()  # intialise RMSE
	xloc = c() # initialise x ticks
	for(i in 1:numcalcs){
		windowsize=windowstepsize*i 			# calculate window size:
      numdump=tsteps %% windowsize # "%%" is modulo
		# reduce data set size if necessary for reshaping:
		nwindows=(tsteps-numdump)/windowsize # number of windows
		# Reshape data:
		cdat=matrix(cable_data[1:(tsteps-numdump)],windowsize,nwindows)
		odat=matrix(obs_data[1:(tsteps-numdump)],windowsize,nwindows)
      Lav_data <- av_window(nwindows, cdat, odat ) 
		# Perform least squares regression (b/w cable ad obs):
      rgrs_data <-rgrs_func( Lav_data, nwindows )  
      mwidth[i] = rgrs_data$mwidth
      mvals[i] = rgrs_data$mvals
		rvals[i]=corr(Lav_data$flav)^2 # get correlation cofficient:
		rmse[i]=sqrt(mean(Lav_data$sqerr))	# Calculate RMSE:
		xloc[i]=windowsize/tstepinday
		rm( cdat, odat, Lav_data ) # clear variables which will have different lengths
  }
  # Work out x-axis ticks:
	tix = c(0,as.integer(maxwindow/4),as.integer(2*maxwindow/4),
	      as.integer(3*maxwindow/4),maxwindow)
   #do the three plots 
   RMSE_plot ( xloc, rmse, maxwindow, sitename, flux, l, tix ) 
   rgrs_gradient_plot ( xloc, mvals, mwidth, sitename, flux, l, tix, numcalcs ) 
   rgrs_corr_plot ( xloc, rvals, sitename, flux, l, tix ) 
}

chk_obs <-function( onc ) {
        if (is.character(rdata$Lodata$onc$var$Rnet$name)) {
          flux=c('Rnet','Qle','Qh','NEE')
          units=c('W/m2','W/m2','W/m2','umol/m2/s') # their units
          layout(matrix(1:12,3,4))
          Rnet = TRUE
        } else {
          flux=c('Qle','Qh','NEE')
          units=c('W/m2','W/m2','umol/m2/s') # their units
          layout(matrix(1:9,3,3))
          Rnet = FALSE
        }
        list( flux=flux, units=units, layout=layout, Rnet=Rnet )
}

av_window<-function(nwindows, cdat, odat ) {
	sqerr=c()
	flav=matrix(0,nwindows,2) # init
	for(j in 1:nwindows){ # calculate average flux for each window
	   flav[j,1]=mean(cdat[,j])	# vector of averages
	   flav[j,2]=mean(odat[,j])	# vector of averages
   	sqerr[j]=(flav[j,1]-flav[j,2])^2
	}
   list( flav=flav, sqerr=sqerr )
}

rgrs_func <-function( fLav_data, nwindows ) { 
	rgrs = lsfit(fLav_data$flav[,2],fLav_data$flav[,1])
	mvals = rgrs$coef[[2]] # store gradient values for plot
	# Get 95% confidence interval:
	linval=c() ; lindev=c() ;	cabdev=c() # init
	meancab=mean(fLav_data$flav[,1])
	for(k in 1:nwindows){
		linval[k]=rgrs$coef[[2]]*fLav_data$flav[k,1]+rgrs$coef[[1]]
		lindev[k]=(linval[k]-fLav_data$flav[k,2])^2
		cabdev[k]=(fLav_data$flav[k,1]-meancab)^2
	}
	# 95% confidence interval:
	mwidth=sqrt(sum(lindev)/nwindows)*1.96/sqrt(sum(cabdev))
   list( mvals=mvals, mwidth=mwidth )
}

RMSE_plot <-function( xloc, rmse, maxwindow, sitename, flux, l, tix) {
 	# Draw RMSE plot:
	plot(xloc,rmse,type="l",xaxt="n",xlab='Averaging window size (days)',
	   ylab = 'RMSE',lwd=1,col='black')
	axis(1,at=tix,labels=as.character(tix))
	title(paste(sitename,'- RMSE for av.',flux[l])) # add title
}

rgrs_gradient_plot <- function( xloc, mvals, mwidth, sitename, flux, l, tix,
   numcalcs ) {
	# Draw gradient plot:
	plot(xloc,mvals,type="l",xaxt="n",xlab='Averaging window size (days)',
		ylab='CABLE vs. obs reg. grad.',lwd=1,col='salmon',
		ylim=c(min(mvals-mwidth),max(mvals+mwidth)))
		# Sort out confidence polygon:
	polyX=c(xloc,xloc[numcalcs:1])
	polyY=c(mvals+mwidth,mvals[numcalcs:1]-mwidth[numcalcs:1])
		polygon(polyX,polyY,col='grey')
	lines(xloc,mvals,lwd=1,col='salmon')		
	axis(1,at=tix,labels=as.character(tix))
	title(paste(sitename,'- gradient for av.',flux[l])) # add title
}

rgrs_corr_plot <- function ( xloc, rvals, sitename, flux, l, tix ) {
	# Draw correlation plot:
	plot(xloc,rvals,type="l",xaxt="n",xlab='Averaging window size (days)',
			ylab=expression(R^2),lwd=1,col='black')
	axis(1,at=tix,labels=as.character(tix))
	title(paste(sitename,'- cor. coeff. for av.',flux[l])) # add title
}
