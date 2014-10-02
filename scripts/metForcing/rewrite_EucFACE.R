# start from a clean slate
rm(list=ls())

# load required libraries
library(ncdf)  # using netcdf files

# set the working directory
setwd("c:/DataWorkshop_R/")

# open met file to read in all values
infile <- open.ncdf("./EucFACE/EucFACE_forcing_1948-2008_ALMA.nc",readunlim=FALSE)

lon      <- get.var.ncdf(infile,"nav_lon")
lat      <- get.var.ncdf(infile,"nav_lat")
GMT      <- get.var.ncdf(infile,"time")
Tair     <- get.var.ncdf(infile,"Tair")
wind     <- get.var.ncdf(infile,"wind")
SWdown   <- get.var.ncdf(infile,"SWdown")
LWdown   <- get.var.ncdf(infile,"LWdown")
Rain     <- get.var.ncdf(infile,"Rainf")
Qair     <- get.var.ncdf(infile,"Qair")
Pres     <- get.var.ncdf(infile,"PSurf")

# shorten timeseries to start from 1949-01-01 00:00:00 local time:
# Firstly, time steps in 1948 = 17568  (leap year, i.e. 17520+48)
# GMT and local time differs by 10 hours, i.e. 20 time steps
# therefore, remove first 17568-20=17548 steps so as to start from 1949
# Also, remove last 20 steps so that the timeseries do not run into year 2009
localTime <- GMT[17549:1069468]-17548*1800

# preparation for the new netcdf file
xvalue <- as.vector(1:1,mode="integer")
yvalue <- as.vector(1:1,mode="integer")
zvalue <- as.vector(1:1,mode="integer")
ALLtime <- as.vector(1:1051920,mode="integer")

# define dimension
xDim <- dim.def.ncdf("x","",xvalue,unlim=FALSE,create_dimvar=TRUE)
yDim <- dim.def.ncdf("y","",yvalue,unlim=FALSE,create_dimvar=TRUE)
zDim <- dim.def.ncdf("z","",zvalue,unlim=FALSE,create_dimvar=TRUE)
tstpDim <- dim.def.ncdf("tstp","timesteps since 1949-01-01 00:00:00",ALLtime,
                        unlim=TRUE,create_dimvar=TRUE)

# define variable for new netcdf file
lonDef  <- var.def.ncdf("nav_lon","degrees_east",dim=list(xDim,yDim),
           missval=-999,prec="single",longname="Longitude")
latDef  <- var.def.ncdf("nav_lat","degrees_north",dim=list(xDim,yDim),
           missval=-999,prec="single",longname="Latitude")
timeDef <- var.def.ncdf("time","seconds since 1949-01-01 00:00:00",
           dim=tstpDim,missval=-999,prec="double",longname="Local Time")
TairDef <- var.def.ncdf("Tair","K",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="Near surface air temperature")
windDef <- var.def.ncdf("wind","m/s",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="Near surface wind speed")
SWdownDef <- var.def.ncdf("SWdown","W/m2",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="Downward shortwave radiation")
LWdownDef <- var.def.ncdf("LWdown","W/m2",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="Downward shortwave radiation")
RainfDef <- var.def.ncdf("Rainf","kg/m2/s",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="precipitation")
QairDef <- var.def.ncdf("Qair","kg/kg",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="specific humidity")
PSurfDef <- var.def.ncdf("PSurf","Pa",dim=list(xDim,yDim,zDim,tstpDim),
           missval=-999,prec="single",longname="Surface pressure")

# open output netcdf file
outfile <- create.ncdf("./EucFACE/EucFACE_forcing_1949-2008_ALMA.nc",
           list(lonDef,latDef,timeDef,TairDef,windDef,SWdownDef,LWdownDef,
                RainfDef,QairDef,PSurfDef))

# write global attribute
att.put.ncdf(outfile,0,"Production","4 July 2013")
att.put.ncdf(outfile,0,"Site_name","EucFACE")
att.put.ncdf(outfile,0,"Length_of_dataset","60 years")
att.put.ncdf(outfile,0,"Contact","Bernard.Pak@csiro.au")
att.put.ncdf(outfile,0,"Source","Princeton 3-hourly dataset interpolated to 30 min, then extracting the grid point closest to the EucFACE expt, changing from GMT to local time")

# write variables
put.var.ncdf(outfile,lonDef, lon)
put.var.ncdf(outfile,latDef, lat)
put.var.ncdf(outfile,timeDef,localTime)
put.var.ncdf(outfile,TairDef,    Tair[17549:1069468])
put.var.ncdf(outfile,windDef,    wind[17549:1069468])
put.var.ncdf(outfile,SWdownDef,SWdown[17549:1069468])
put.var.ncdf(outfile,LWdownDef,LWdown[17549:1069468])
put.var.ncdf(outfile,RainfDef,   Rain[17549:1069468])
put.var.ncdf(outfile,QairDef,    Qair[17549:1069468])
put.var.ncdf(outfile,PSurfDef,   Pres[17549:1069468])

close.ncdf(outfile)
