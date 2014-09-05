#!/bin/bash
# jtk561 (Jatin.Kala.JK@gmail.com) - script to convert Albert's Streamflow data from geotiff to netcdf using gdal, units and more meaningful variable names and descriptions are added using nco tools
module load gdal
module load nco
##########################################################
# input path to GSCD data
dir_GSDC="/g/data1/wd9/BenchMarking/GSCD/RAW/GSCD_v1.8"
# output path to write the nc version
dir_out="/g/data1/wd9/BenchMarking/GSCD/NC/GSCD_v1.8"
##########################################################
# code
# first convert to nc
for f in ${dir_GSDC}/*.tif
    do
       echo "Processing: $f"
       # get fine name without path and extension
       f_only=$(basename $f .tif)
       # delete output nc file in case already exist
       rm -f $dir_out/$f_only.nc
       # convert to nc using gdal_translate
       gdal_translate -of netCDF $f $dir_out/$f_only.nc
    done

# now add units and a bit more description etc
for f in ${dir_out}/BFI1*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Baseflow index" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo
       ncatted -a description,$f_only,c,c,'Baseflow index, defined as ratio of long-term baseflow to total Q (Smakhtin 2001). Computed from daily Q data using the recursive digital filter of Van Dijk (2010) with the window size set to five days (cf. Pena-Arancibia et al. 2010). See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/BFI2*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Baseflow index" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo
       ncatted -a description,$f_only,c,c,'Baseflow index, computed from daily Q data following the local-minimum method described in Pettyjohn and Henning (1979) and Sloto and Crouse (1996), with the duration of surface runoff (N) set to five days. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/BFI3*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Baseflow index" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo      
       ncatted -a description,$f_only,c,c,'Baseflow index, computed from daily Q data following the sliding-interval method described in Pettyjohn and Henning (1979) and Sloto and Crouse (1996), with N set to three days. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/BFI4*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Baseflow index" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo 
       ncatted -a description,$f_only,c,c,'Baseflow index, computed from daily Q data following the procedure described in Institute of Hydrology (1980) and Gustard et al. (1992), which takes the minima at five-day non-overlapping intervals and subsequently connects the valleys in this series of minima to generate baseflow. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/K*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Baseflow recession constant" ${f}_tempo
       ncatted -a units,$f_only,c,c,"d-1" ${f}_tempo 
       ncatted -a description,$f_only,c,c,'Baseflow recession constant, defined as the rate of baseflow decay (Vogel and Kroll 1996). Computed from daily Q data following Van Dijk (2010), with the window size set to five days and days with zero flow ignored. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/Q?.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Daily flow percentiles (exceedance probability)" ${f}_tempo
       ncatted -a units,$f_only,c,c,"mm d-1" ${f}_tempo      
       ncatted -a description,$f_only,c,c,'Daily flow percentiles (exceedance probability) computed from daily Q data. The number refers to the percentage of time that the flow is exceeded. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/Q??.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Daily flow percentiles (exceedance probability)" ${f}_tempo
       ncatted -a units,$f_only,c,c,"mm d-1" ${f}_tempo
       ncatted -a description,$f_only,c,c,'Daily flow percentiles (exceedance probability) computed from daily Q data. The number refers to the percentage of time that the flow is exceeded. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/T50*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"The day of the water year marking the timing of the center of mass of streamflow" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo
       ncatted -a description,$f_only,c,c,'The day of the water year marking the timing of the center of mass of Q (Stewart et al. 2005). The water year is defined in this study sensu the United States Geological Survey (USGS) as the 12-month period from October to September in the Northern Hemisphere and April to March in the Southern Hemisphere. T50 can be computed from both daily Q data and monthly Q data linearly interpolated to daily values.. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/RC*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Runoff coefficient" ${f}_tempo
       ncatted -a units,$f_only,c,c,"-" ${f}_tempo 
       ncatted -a description,$f_only,c,c,'Runoff coefficient, the ratio of Q to mean annual precipitation. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

for f in ${dir_out}/QMEAN*.nc
    do
       echo "Processing: $f"
       f_only=$(basename $f .nc)
       rm -f ${f}_tempo
       ncrename -v Band1,$f_only $f ${f}_tempo
       ncatted -a _FillValue,$f_only,o,f,NaN ${f}_tempo
       ncatted -a _FillValue,$f_only,m,f,-1.0e10 ${f}_tempo
       ncatted -a long_name,$f_only,o,c,"Mean annual streamflow per unit area" ${f}_tempo
       ncatted -a units,$f_only,c,c,"mm year-1" ${f}_tempo              
       ncatted -a description,$f_only,c,c,'Mean annual streamflow per unit area. See Beck et al. 2014, J. Hyrdrometeorology (submitted) for more details.' ${f}_tempo
       rm $f
       mv ${f}_tempo $f
    done

