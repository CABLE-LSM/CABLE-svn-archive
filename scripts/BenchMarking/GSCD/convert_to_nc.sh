#!/bin/bash
# jtk561 (Jatin.Kala.JK@gmail.com) - script to convert Albert's Streamflow data from geotiff to netcdf using gdal
module load gdal
##########################################################
# input path to GSCD data
dir_GSDC="/g/data1/wd9/BenchMarking/GSCD/RAW/GSCD_v1.8"
# output path to write the nc version
dir_out="/g/data1/wd9/BenchMarking/GSCD/NC/GSCD_v1.8"
##########################################################
# code
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
