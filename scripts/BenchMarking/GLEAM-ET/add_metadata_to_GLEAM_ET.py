# -*- coding: utf-8 -*-
# Jatin Kala (Jatin.Kala.JK@gmail.com)
# Purpose of this script is to add metadata to the original ET files provided by Diego Miralles, as they lack lat,lon,units, and proper CF-compliant time variables, the utf-8 at the top is required for use of non-ascii characters when adding references to netcdf file global attribute
import netCDF4 as nc
import glob
import numpy as np
import datetime
import calendar 
import os
############### USER-DEFINED FUNCTIONS ################################################
def get_cf_time_from_year(year,isleap_cal):
   """Return a time variable folloling the CF convention of "days from 1/1/1970" for every Julian day of a given year. Inputs:
   year(integer)
   isleap_cal(string) - must be "leap" or "noleap"
   """
   year=int(year)
   ref_seconds = datetime.datetime(1970,1,1)
   mm = np.zeros(12,dtype=np.int)
   if isleap_cal == "leap":
             cc = 0
             for i in mm:
                 mm[cc] = calendar.monthrange(year,cc+1)[1]
                 cc = cc + 1
   elif isleap_cal == "noleap":
             cc = 0
             for i in mm:
                 mm[cc] = calendar.monthrange(1991,cc+1)[1] #hard code year to non-leap
                 cc = cc + 1  
   else:
       print("isleap_cal must be either set to leap or noleap, exiting")
       exit()
   no_days = int(np.sum(mm))
   time = np.zeros(no_days,dtype='f')
   cc = 0
   cc1 = 0
   for i in mm:
       dd = np.arange(1,i+1,1)
       for j in dd:
            t = datetime.datetime(year,cc+1,j,0,0)
            time[cc1] = ((t-ref_seconds).total_seconds())/(24.0*60.0*60.0)
            cc1 = cc1 + 1
       cc = cc + 1
   return(time)
##########################################################################

############### INPUTS ###################################################
# define dir to v2A files
dir_v2A = "/srv/ccrc/data11/z3381484/GLEAM-ET/v2A/"
# define dir to write new v2A files
dir_write_v2A = "/srv/ccrc/data11/z3381484/GLEAM-ET/v2A_CF_Compliant/"
# define dir to v2B files
dir_v2B = "/srv/ccrc/data11/z3381484/GLEAM-ET/v2B/"
# define dir to write new v2B files
dir_write_v2B = "/srv/ccrc/data11/z3381484/GLEAM-ET/v2B_CF_Compliant/"

# define the domain info from Diego's email
lat_top = -90.0
lat_bot = 90.0
lon_left = -180.0
lon_right = 180.0
dom_res = 0.25
# Strings for global file attributes
str_title = '%s%s%s' % ('GLEAM ET data-set from Diego Miralles (Diego.Miralles@bristol.ac.uk), modified by Jatin Kala (J.Kala@unsw.edu.au or Jatin.Kala.JK@gmail.com), on: ',datetime.datetime.now(),' to include lat/lon, units, _FillValue, and CF-compliant time dimension')
str_references = """Miralles, D.G., J.H. Gash, T.R.H. Holmes, R.A.M. De Jeu & A.J. Dolman (2010). Global canopy interception from satellite observations. J. Geophys. Res. D: Atmos., 115, D16122.
Miralles, D. G., R. A. M. De Jeu, J. H. Gash, T. R. H. Holmes & A. J. Dolman (2011), Magnitude and variability of land evaporation and its components at the global scale, Hydrol. Earth Syst. Sci., 15(3), 967-981.
Miralles, D. G., T. R. H. Holmes, R. A. M. De Jeu, J. H. Gash, A. G. C. A. Meesters, and A. J. Dolman (2011), Global land-surface evaporation estimated from satellite-based observations, Hydrol. Earth Syst. Sci., 15(2), 453-469.
Miralles, D. G., et al. (2014), El Niño–La Niña cycle and recent trends in continental evaporation, Nature Clim. Change, 4, 122-126."""
CF_ver = 'CF-1.6'
########################################################################## 

########## CODE ##########################################################
# define lat-lon arrays, as per Diego's email
latitude = np.arange(lat_bot,lat_top,-dom_res,dtype='f')
longitude = np.arange(lon_left,lon_right,dom_res,dtype='f')

# process v2A files
# get files 
files_v2A = glob.glob('%s%s' % (dir_v2A,"E_*.nc"))
# loop through files
for f in files_v2A:
    # open file
    print('%s%s' % ("Processing: ", f))
    file_open_read = nc.Dataset(f,mode='r')
    # get ET out
    ET = file_open_read.variables['E']
    del(file_open_read)
    # work out the time, get year from string
    year = (f[len(f)-7:len(f)-3])
    # define CF_compliant time dim for each day of the year
    time_cf = get_cf_time_from_year(int(year),"leap")
    file_name_write = '%s%s%s%s' % (dir_write_v2A,"GLEAM-ET_",year,".nc")
    file_name_compress = '%s%s%s%s' % (dir_write_v2A,"GLEAM-ET_C_",year,".nc")
    if os.path.exists(file_name_write):
        os.remove(file_name_write)
    fileobj = nc.Dataset(file_name_write,mode='w')
    fileobj.createDimension('latitude', len(latitude))
    fileobj.createDimension('longitude', len(longitude))
    fileobj.createDimension('time',len(time_cf))
    fileobj.createDimension('nv',2)
    lat_var = fileobj.createVariable('latitude', 'f', ('latitude',))
    lat_var.units = 'degrees_north'
    lat_var.long_name = 'latitude'
    lat_var[:] = latitude
    lon_var = fileobj.createVariable('longitude', 'f', ('longitude',))
    lon_var.units = 'degrees_east'
    lon_var.long_name = 'longitude'
    lon_var[:] = longitude
    time_var = fileobj.createVariable('time','f',('time',))
    time_var.units = 'days since 1970-1-1 0:0:0'
    time_var.calendar = 'gregorian'
    time_var.long_name = 'time'
    time_var.bounds = 'time_bnds'
    time_var[:] = time_cf[:]
    bnds_var = fileobj.createVariable('time_bnds', 'f', ('time','nv',))
    bnds_var.units = time_var.units
    bnds_var.long_name = 'Time boundaries'
    bnds = np.zeros((len(time_cf),2))
    bnds[:,0] = time_cf - 0.5
    bnds[:,1] = time_cf + 0.5
    bnds_var[:,:] = bnds[:,:]
    del(bnds)
    del(time_cf)
    del(time_var)
    data_var = fileobj.createVariable('ET', 'd', ('time','latitude','longitude'),fill_value = -9)
    data_var.units = 'mm day-1'
    data_var.long_name = 'Evapotranspiration'
    data_var.standard_name = data_var.long_name
    data_var.cell_methods = 'time: mean'
    data_var[:,:,:] = np.swapaxes(ET[:,:,:],1,2)
    del(ET)
    del(data_var)
    fileobj.history = str_title 
    fileobj.references = str_references
    fileobj.Conventions = CF_ver
    fileobj.close()
    os.system('%s%s%s%s' % ('nccopy -d 9 ', file_name_write, ' ', file_name_compress))
    os.system('%s%s%s%s' % ('mv ', file_name_compress, ' ', file_name_write))

del(files_v2A)

# process v2B files, these have E(total ET) and T(Transpiration) files
# get files 
files_v2B = glob.glob('%s%s' % (dir_v2B,"E_*.nc"))
# loop through files
for f in files_v2B:
    # open file
    print('%s%s' % ("Processing: ", f))
    file_open_read = nc.Dataset(f,mode='r')
    # get ET out
    ET = file_open_read.variables['E']
    del(file_open_read)
    # work out the time, get year from string
    year = (f[len(f)-7:len(f)-3])
    # define CF_compliant time dim for each day of the year
    time_cf = get_cf_time_from_year(int(year),"leap")
    file_name_write = '%s%s%s%s' % (dir_write_v2B,"GLEAM-ET_",year,".nc")
    file_name_compress = '%s%s%s%s' % (dir_write_v2B,"GLEAM-ET_C_",year,".nc")
    if os.path.exists(file_name_write):
        os.remove(file_name_write)
    fileobj = nc.Dataset(file_name_write,mode='w')
    fileobj.createDimension('latitude', len(latitude))
    fileobj.createDimension('longitude', len(longitude))
    fileobj.createDimension('time',len(time_cf))
    fileobj.createDimension('nv',2)
    lat_var = fileobj.createVariable('latitude', 'f', ('latitude',))
    lat_var.units = 'degrees_north'
    lat_var.long_name = 'latitude'
    lat_var[:] = latitude
    lon_var = fileobj.createVariable('longitude', 'f', ('longitude',))
    lon_var.units = 'degrees_east'
    lon_var.long_name = 'longitude'
    lon_var[:] = longitude
    time_var = fileobj.createVariable('time','f',('time',))
    time_var.units = 'days since 1970-1-1 0:0:0'
    time_var.calendar = 'gregorian'
    time_var.long_name = 'time'
    time_var.bounds = 'time_bnds'
    time_var[:] = time_cf[:]
    bnds_var = fileobj.createVariable('time_bnds', 'f', ('time','nv',))
    bnds_var.units = time_var.units
    bnds_var.long_name = 'Time boundaries'
    bnds = np.zeros((len(time_cf),2))
    bnds[:,0] = time_cf - 0.5
    bnds[:,1] = time_cf + 0.5
    bnds_var[:,:] = bnds[:,:]
    del(bnds)
    del(time_cf)
    del(time_var)
    data_var = fileobj.createVariable('ET', 'd', ('time','latitude','longitude'),fill_value = -9)
    data_var.units = 'mm day-1'
    data_var.long_name = 'Evapotranspiration'
    data_var.standard_name = data_var.long_name
    data_var.cell_methods = 'time: mean'
    data_var[:,:,:] = np.swapaxes(ET[:,:,:],1,2)
    del(ET)
    del(data_var)
    fileobj.history = str_title 
    fileobj.references = str_references
    fileobj.Conventions = CF_ver
    fileobj.close()
    os.system('%s%s%s%s' % ('nccopy -d 9 ', file_name_write, ' ', file_name_compress))
    os.system('%s%s%s%s' % ('mv ', file_name_compress, ' ', file_name_write))

# get files 
files_v2B = glob.glob('%s%s' % (dir_v2B,"T_*.nc"))
# loop through files
for f in files_v2B:
    # open file
    print('%s%s' % ("Processing: ", f))
    file_open_read = nc.Dataset(f,mode='r')
    # get T out
    T = file_open_read.variables['T']
    del(file_open_read)
    # work out the time, get year from string
    year = (f[len(f)-7:len(f)-3])
    # define CF_compliant time dim for each day of the year
    time_cf = get_cf_time_from_year(int(year),"leap")
    file_name_write = '%s%s%s%s' % (dir_write_v2B,"GLEAM-T_",year,".nc")
    file_name_compress = '%s%s%s%s' % (dir_write_v2B,"GLEAM-T_C_",year,".nc")
    if os.path.exists(file_name_write):
        os.remove(file_name_write)
    fileobj = nc.Dataset(file_name_write,mode='w')
    fileobj.createDimension('latitude', len(latitude))
    fileobj.createDimension('longitude', len(longitude))
    fileobj.createDimension('time',len(time_cf))
    fileobj.createDimension('nv',2)
    lat_var = fileobj.createVariable('latitude', 'f', ('latitude',))
    lat_var.units = 'degrees_north'
    lat_var.long_name = 'latitude'
    lat_var[:] = latitude
    lon_var = fileobj.createVariable('longitude', 'f', ('longitude',))
    lon_var.units = 'degrees_east'
    lon_var.long_name = 'longitude'
    lon_var[:] = longitude
    time_var = fileobj.createVariable('time','f',('time',))
    time_var.units = 'days since 1970-1-1 0:0:0'
    time_var.calendar = 'gregorian'
    time_var.long_name = 'time'
    time_var.bounds = 'time_bnds'
    time_var[:] = time_cf[:]
    bnds_var = fileobj.createVariable('time_bnds', 'f', ('time','nv',))
    bnds_var.units = time_var.units
    bnds_var.long_name = 'Time boundaries'
    bnds = np.zeros((len(time_cf),2))
    bnds[:,0] = time_cf - 0.5
    bnds[:,1] = time_cf + 0.5
    bnds_var[:,:] = bnds[:,:]
    del(bnds)
    del(time_cf)
    del(time_var)
    data_var = fileobj.createVariable('T', 'd', ('time','latitude','longitude'),fill_value = -9)
    data_var.units = 'mm day-1'
    data_var.long_name = 'Vegetation Transpiration'
    data_var.standard_name = data_var.long_name
    data_var.cell_methods = 'time: mean'
    data_var[:,:,:] = np.swapaxes(T[:,:,:],1,2)
    del(T)
    del(data_var)
    fileobj.history = str_title 
    fileobj.references = str_references
    fileobj.Conventions = CF_ver
    fileobj.close()
    os.system('%s%s%s%s' % ('nccopy -d 9 ', file_name_write, ' ', file_name_compress))
    os.system('%s%s%s%s' % ('mv ', file_name_compress, ' ', file_name_write))

del(files_v2B)
