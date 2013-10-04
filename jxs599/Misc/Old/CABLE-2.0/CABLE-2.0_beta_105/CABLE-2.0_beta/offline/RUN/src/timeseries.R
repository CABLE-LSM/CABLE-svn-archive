##################################################################################
# This function simply plots a timeseries of the entire dataset, for a wide
# range of variables.
# png and jpg types: width scaled down from 3564 to 2000,
#                    height from 2520 to 1414
#                    pointsize from 24 to 12 (BP, Oct 2007)
timeseries=function(sitename,cablefile,outtype,timestep_start,timestep_stop) {
	outfilename=paste(sitename,'Series.',outtype,sep='') # name of plot if not to screen
	if(outtype=='pdf'){
		pdf(file=outfilename,paper='a4r',width=11,height=8)
	}else if(outtype=='ps'){
		postscript(file=outfilename,paper='special',width=11,height=8)
	}else if(outtype=='png'){
		png(file=outfilename,width=2000,height=1414,pointsize=12)
	}else if(outtype=='jpg'){
		jpeg(file=outfilename,width=2000,height=1414,pointsize=12)
	}
	# open CABLE file:
	cnc=open.ncdf(cablefile,readunlim=FALSE) # open CABLE file
	# set number of time steps in plot (data set length):
	if(timestep_start==-1){
		timestep_start=1
		timestep_stop = cnc$dim$time[[2]]
	}
	# sort out x-axis markings for plots:
	ntsteps=timestep_stop-timestep_start+1
	xax=c(timestep_start,as.integer(timestep_start+ntsteps/5),
		as.integer(timestep_start+2*ntsteps/5),
		as.integer(timestep_start+3*ntsteps/5),
		as.integer(timestep_start+4*ntsteps/5),timestep_stop)
	xloc=c(timestep_start:timestep_stop)	
	# set layout of plots (10 in 5*2 matrix):
	layout(matrix(1:12,3,4,byrow=TRUE))
	# Read data for first plot:---------------------------------------
	SoilTemp=get.var.ncdf(cnc,'SoilTemp',start=c(1,1,1,timestep_start),
		count=c(1,1,6,ntsteps)) # read soil temp
	# First plot:	
	plot(xloc,SoilTemp[1,],type="l",xaxt="n",xlab='Time steps',
		ylab='Degrees K',lwd=0.5,col='green',
		ylim=c(min(SoilTemp),max(SoilTemp)))
	lines(xloc,SoilTemp[2,],lwd=0.5,col='forestgreen')
	lines(xloc,SoilTemp[3,],lwd=0.5,col='darkslategrey')
	lines(xloc,SoilTemp[4,],lwd=0.5,col='deepskyblue4')
	lines(xloc,SoilTemp[5,],lwd=0.5,col='blue')
	lines(xloc,SoilTemp[6,],lwd=0.5,col='darkblue')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Soil temperatures -',sitename))
	legend('topleft',c('1','2','3','4','5','6'),lty=1,col=c		('green','forestgreen','darkslategrey','deepskyblue4', 		'blue','darkblue'),lwd=0.5)
	# Read data for second plot:--------------------------------------------
	#Tair=get.var.ncdf(cnc,'Tair',start=c(1,1,1,timestep_start),
	Tair=get.var.ncdf(cnc,'Tair',start=c(1,1,timestep_start),
		#count=c(1,1,1,ntsteps)) # read air temp
		count=c(1,1,ntsteps)) # read air temp
	VegT=get.var.ncdf(cnc,'VegT',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read veg temperature
	# Second plot:	
	plot(xloc,SoilTemp[1,],type="l",xaxt="n",xlab='Time steps',
		ylab='Degrees K',lwd=0.5,col='peru',
		ylim=c(min(SoilTemp[1,],SoilTemp[2,],Tair,VegT),
		max(SoilTemp[1,],SoilTemp[2,],Tair,VegT)))
	lines(xloc,SoilTemp[2,],lwd=0.5,col='sienna4')
	lines(xloc,Tair,lwd=0.5,col='blue')
	lines(xloc,VegT,lwd=0.5,col='green')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Temperatures -',sitename))
	legend('topleft',c('soil1','soil2','airT','vegT'),lty=1,col=c		('peru','sienna4','blue','green'),lwd=0.5)
	# Read data for third plot:--------------------------------------------
	HSoil=get.var.ncdf(cnc,'HSoil',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read sensible heat from soil
	ESoil=get.var.ncdf(cnc,'ESoil',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*2.5104e6 # read latent heat from soil
	Qg=get.var.ncdf(cnc,'Qg',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read ground heat flux
	# Third plot:	
	plot(xloc,HSoil,type="l",xaxt="n",xlab='Time steps',
		ylab=expression(W/m^2),lwd=0.5,col='red',
		ylim=c(min(HSoil,ESoil,Qg),max(HSoil,ESoil,Qg)))
	lines(xloc,ESoil,lwd=0.5,col='blue')
	lines(xloc,Qg,lwd=0.5,col='peru')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Soil fluxes -',sitename))
	legend('topleft',c('SSens','SLat','Gflux'),lty=1,col=c		('red','blue','brown'),lwd=0.5)	
	# Read data for fourth plot:--------------------------------------------
	Evap=get.var.ncdf(cnc,'Evap',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*2.5104e6 # read total evapotranspiration
	ECanop=get.var.ncdf(cnc,'ECanop',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*2.5104e6 # read wet canopy evaporation
	TVeg=get.var.ncdf(cnc,'TVeg',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*2.5104e6 # read vegetation transpiration
	HVeg=get.var.ncdf(cnc,'HVeg',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read sensible heat from vegetation
	# Fourth plot:	
	plot(xloc,Evap,type="l",xaxt="n",xlab='Time steps',
		ylab=expression(W/m^2),lwd=0.5,col='black',
		ylim=c(min(Evap,ECanop,TVeg,HVeg),max(Evap,ECanop,TVeg,HVeg)))
	lines(xloc,ECanop,lwd=0.5,col='blue')
	lines(xloc,TVeg,lwd=0.5,col='green')
	lines(xloc,HVeg,lwd=0.5,col='red')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Canopy fluxes -',sitename))
	legend('topleft',c('TotEvapotr','VEvap','VTransp','VSens'),
		lty=1,col=c	('black','blue','green','red'),lwd=0.5)
	# Read data for fifth plot:---------------------------------------
	SoilMoist=get.var.ncdf(cnc,'SoilMoist',start=c(1,1,1,timestep_start),
		count=c(1,1,6,ntsteps)) # read soil moisture
	zse=get.var.ncdf(cnc,'zse',start=c(1,1,1),count=c(1,1,6)) # read soildepth
        SMoist=SoilMoist
# removed the next 5 lines as there is no need to change units now.
## Change units from kg/m^2 to m^3/m^3:
#SMoist=matrix(0,6,ntsteps)
#for(i in 1:6){
#	SMoist[i,]=SoilMoist[i,]/zse[i]/1000
#}
	# load filed capacity, saturation and wilting point values:
	sfc=vector(mode='numeric',length=ntsteps)
	ssat=vector(mode='numeric',length=ntsteps)
	swilt=vector(mode='numeric',length=ntsteps)
	sfc[]=get.var.ncdf(cnc,'sfc',start=c(1,1),count=c(1,1))
	ssat[]=get.var.ncdf(cnc,'ssat',start=c(1,1),count=c(1,1))
	swilt[]=get.var.ncdf(cnc,'swilt',start=c(1,1),count=c(1,1))
	# Fifth plot:	
	plot(xloc,SMoist[1,],type="l",xaxt="n",xlab='Time steps',
		ylab=expression(m^3/m^3),lwd=0.5,col='green',
		ylim=c(min(SMoist,sfc,swilt,ssat),max(SMoist,sfc,swilt,ssat)))
	lines(xloc,SMoist[2,],lwd=0.5,col='forestgreen')
	lines(xloc,SMoist[3,],lwd=0.5,col='darkslategrey')
	lines(xloc,SMoist[4,],lwd=0.5,col='deepskyblue4')
	lines(xloc,SMoist[5,],lwd=0.5,col='blue')
	lines(xloc,SMoist[6,],lwd=0.5,col='darkblue')
	lines(xloc,sfc,lwd=0.5,col='grey')
	lines(xloc,swilt,lwd=0.5,col='grey')
	lines(xloc,ssat,lwd=0.5,col='grey')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Soil moisture -',sitename))
	legend('topleft',c('1','2','3','4','5','6'),lty=1,col=c		('green','forestgreen','darkslategrey','deepskyblue4', 		'blue','darkblue'),lwd=0.5)
	# Read data for sixth plot:--------------------------------------------
	Rainf=get.var.ncdf(cnc,'Rainf',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*3600 # read rainfall
	CanopInt=get.var.ncdf(cnc,'CanopInt',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read canopy water storage
	Qs=get.var.ncdf(cnc,'Qs',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*3600 # read runoff
	Qsb=get.var.ncdf(cnc,'Qsb',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps))*3600 # read deep drainage
	# Sixth plot:	
	plot(xloc,Rainf,type="l",xaxt="n",xlab='Time steps',
		ylab='mm and mm/h',lwd=0.5,col='blue',
		ylim=c(min(Rainf,CanopInt,Qs,Qsb),max(Rainf,CanopInt,Qs,Qsb)))
	lines(xloc,CanopInt,lwd=0.5,col='green')
	lines(xloc,Qs,lwd=0.5,col='cyan')
	lines(xloc,Qsb,lwd=0.5,col='brown4')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Precip/runoff -',sitename))
	legend('topleft',c('Precip','CanStore','Runoff','Drainage'),
		lty=1,col=c	('blue','green','cyan','brown4'),lwd=0.5)
	# Read data for seventh plot:--------------------------------------------
	SWnet=get.var.ncdf(cnc,'SWnet',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read net shortwave
	LWnet=get.var.ncdf(cnc,'LWnet',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read net longwave
	SWdown=get.var.ncdf(cnc,'SWdown',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read runoff
	LWdown=get.var.ncdf(cnc,'LWdown',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read deep drainage
	# Seventh plot:	
	plot(xloc,SWdown,type="l",xaxt="n",xlab='Time steps',
		ylab=expression(W/m^2),lwd=0.5,col='black',
		ylim=c(min(SWdown,SWnet,LWdown,LWnet),
		max(SWdown,SWnet,LWdown,LWnet)))
	lines(xloc,SWnet,lwd=0.5,col='blue')
	lines(xloc,LWdown,lwd=0.5,col='green')
	lines(xloc,LWnet,lwd=0.5,col='forestgreen')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Radiation -',sitename))
	legend('topleft',c('SWdown','SWnet','LWdown','LWnet'),
		lty=1,col=c	('black','blue','green','forestgreen'),lwd=0.5)
	# Read data for eighth plot:--------------------------------------------
	NEE=get.var.ncdf(cnc,'NEE',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read net ecosystem exchange
	AutoResp=get.var.ncdf(cnc,'AutoResp',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read autotrophic respiration
	HeteroResp=get.var.ncdf(cnc,'HeteroResp',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read heterotrophic respiration
	GPP=get.var.ncdf(cnc,'GPP',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read net primary production
	# Eighth plot:	
	plot(xloc,NEE,type="l",xaxt="n",xlab='Time steps',
		ylab=expression(paste(mu,mol/m^2/s)),lwd=0.5,col='blue',
		ylim=c(min(NEE,AutoResp,HeteroResp,GPP),
		max(NEE,AutoResp,HeteroResp,GPP)))
	lines(xloc,AutoResp,lwd=0.5,col='green')
	lines(xloc,HeteroResp,lwd=0.5,col='cyan')
	lines(xloc,GPP,lwd=0.5,col='brown4')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Carbon fluxes -',sitename))
	legend('topleft',c('NEE','AutoR','HeteroR','GPP'),
		lty=1,col=c	('blue','green','cyan','brown4'),lwd=0.5)
	# Read data for ninth plot:------------------------------------------
	Qle=get.var.ncdf(cnc,'Qle',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read latent heat
	Qh=get.var.ncdf(cnc,'Qh',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read sensible heat
	# Ninth plot:	
	plot(xloc,(SWnet+LWnet),type="l",xaxt="n",xlab='Time steps',
		ylab=expression(W/m^2),lwd=0.5,col='black',
		ylim=c(min(SWnet+LWnet,Qle,Qh,Qg),
		max(SWnet+LWnet,Qle,Qh,Qg)))
	lines(xloc,Qle,lwd=0.5,col='blue')
	lines(xloc,Qh,lwd=0.5,col='red')
	lines(xloc,Qg,lwd=0.5,col='brown4')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Combined fluxes -',sitename))
	legend('topleft',c('Rnet','Lat','Sens','Gfl'),
		lty=1,col=c	('black','blue','red','brown4'),lwd=0.5)
	# Read data for tenth plot:------------------------------------------
	Albedo=get.var.ncdf(cnc,'Albedo',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read albedo
	SnowDepth=get.var.ncdf(cnc,'SnowDepth',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read snow depth
	LAI=get.var.ncdf(cnc,'LAI',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read LAI
	# Tenth plot:	
	plot(xloc,Albedo,type="l",xaxt="n",xlab='Time steps',
		ylab='- and m',lwd=0.5,col='red',
		ylim=c(min(Albedo,SnowDepth,LAI/10),max(Albedo,SnowDepth,LAI/10)))
	lines(xloc,LAI/10,lwd=0.5,col='forestgreen')
	lines(xloc,SnowDepth,lwd=0.5,col='blue')
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Albedo and snow depth -',sitename))
	legend('topleft',c('Albedo','LAI/10','Snowdepth'),
		lty=1,col=c	('red','forestgreen','blue'),lwd=0.5)
	# Read data for eleventh plot:------------------------------------------
	Ebal=get.var.ncdf(cnc,'Ebal',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read albedo
	plot(xloc,Ebal,type="l",xaxt="n",xlab='Time steps',
		ylab=expression(W/m^2),lwd=0.5,col='red',
		ylim=c(min(Ebal),max(Ebal)))
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Cumulative energy balance - ',sitename))
	# Read data for twelfth plot:------------------------------------------
	Wbal=get.var.ncdf(cnc,'Wbal',start=c(1,1,timestep_start),
		count=c(1,1,ntsteps)) # read albedo
	plot(xloc,Wbal,type="l",xaxt="n",xlab='Time steps',
		ylab='mm',lwd=0.5,col='blue',
		ylim=c(min(Wbal),max(Wbal)))
	axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
	title(paste('Cumulative water balance - ',sitename))

	# Close CABLE output netcdf file:	
	close.ncdf(cnc)
	# close graphics file if used:
	if(outtype!='screen')	dev.off()
	# Clear all current local environment variables
	rm(list=ls()) 		
} # End function timeseries
