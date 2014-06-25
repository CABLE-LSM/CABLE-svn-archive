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

def write_lat_lon_time(fileobj,lat,lon,time):
             """Write lat, lon, and time to netcdf file"""
             fileobj = fileobj
             fileobj.createDimension('latitude', len(lat))
             fileobj.createDimension('longitude', len(lon))
             fileobj.createDimension('time',len(time))
             fileobj.createDimension('nv',2)
             lat_var = fileobj.createVariable('latitude', 'f', ('latitude',))
             lat_var.units = 'degrees_north'
             lat_var.long_name = 'latitude'
             lat_var[:] = lat
             lon_var = fileobj.createVariable('longitude', 'f', ('longitude',))
             lon_var.units = 'degrees_east'
             lon_var.long_name = 'longitude'
             lon_var[:] = lon
             time_var = fileobj.createVariable('time','f',('time',))
             time_var.units = 'days since 1970-1-1 0:0:0'
             time_var.calendar = 'gregorian'
             time_var.long_name = 'time'
             time_var.bounds = 'time_bnds'
             time_var[:] = time[:]
             bnds_var = fileobj.createVariable('time_bnds', 'f', ('time','nv',))
             bnds_var.units = time_var.units
             bnds_var.long_name = 'Time boundaries'
             bnds = np.zeros((len(time),2))
             bnds[:,0] = time - 0.5
             bnds[:,1] = time + 0.5
             bnds_var[:,:] = bnds[:,:]
             del(bnds)
             del(time)
             del(time_var)
  
def write_var(fileobj,var_name,var):
             """ write variable to netcdf file"""
             fileobj = fileobj
             data_var = fileobj.createVariable(var_name, 'd', ('time','latitude','longitude'),fill_value = -9)
             data_var.units = 'mm day-1'
             data_var.long_name = (vars_writenames[var_name])
             data_var.standard_name = data_var.long_name
             data_var.cell_methods = 'time: mean'
             data_var[:,:,:] = np.swapaxes(var[:,:,:],1,2)
             del(var)
             del(data_var)
             fileobj.history = str_title
             fileobj.references = str_references
             fileobj.Conventions = CF_ver


##########################################################################

############### INPUTS ###################################################
# define dir to v2A files
dir_v2A = "/srv/ccrc/data34/z3381484/new-GLEAM-ET/GLEAM/v2A/"
# define dir to write new v2A files
dir_write_v2A = "/srv/ccrc/data34/z3381484/new-GLEAM-ET/GLEAM/CF/v2A/"
# define dir to v2B files
dir_v2B = "/srv/ccrc/data34/z3381484/new-GLEAM-ET/GLEAM/v2B/"
# define dir to write new v2B files
dir_write_v2B = "/srv/ccrc/data34/z3381484/new-GLEAM-ET/GLEAM/CF/v2B/"
# vars to process
vars_proc = ['E','E_b','I','T']
vars_writenames = {'E':'Evapotranspiration','E_b':'Transpiration','I':'Vegetation Interception Loss','T':'Bare Soil Evaporation'}

# define the domain as per Doc
lat_top_v2A = -90.0
lat_bot_v2A = 90.0
lat_top_v2B = -59.0
lat_bot_v2B = 59.0
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
######################################################################

# define lat-lon arrays, as per Diego's email
latitude_v2A = np.arange(lat_bot_v2A,lat_top_v2A,-dom_res,dtype='f')
latitude_v2B = np.arange(lat_bot_v2B,lat_top_v2B,-dom_res,dtype='f')
longitude = np.arange(lon_left,lon_right,dom_res,dtype='f')

# loop through vars
for vars_name in vars_proc:
      # get files
      files_v2A = glob.glob('%s%s%s' % (dir_v2A, vars_name,'?????.nc'))
      for f in files_v2A:
             print('%s%s' % ("Processing: ", f))
             # open file
             file_open_read = nc.Dataset(f,mode='r')
             # get var
             var_out = file_open_read.variables[vars_name][:,:,:]
             del(file_open_read)
             # work out the time, get year from string
             year = (f[len(f)-7:len(f)-3])
             # define CF_compliant time dim for each day of the year
             time_cf = get_cf_time_from_year(int(year),"leap")
             # write time,lat,lon, and var to file
             file_name_write = '%s%s%s%s%s%s' % (dir_write_v2A,"GLEAM-",vars_name,"_",year,".nc")
             if os.path.exists(file_name_write):
                   os.remove(file_name_write)
             # create file handle
             fileobj = nc.Dataset(file_name_write,mode='w')
             # write lat, lon, and time
             write_lat_lon_time(fileobj,latitude_v2A,longitude,time_cf)
             del(time_cf)
             # make nan to fill_value = -9
             var_out_no_nan = np.copy(var_out)
             var_out_no_nan[np.isnan(var_out)] = -9
             del(var_out)
             # write var
             write_var(fileobj,vars_name,var_out_no_nan)
             # close file
             fileobj.close()
             del(var_out_no_nan)
             # run compression
             file_name_compress = '%s%s%s%s%s%s%s' % (dir_write_v2A,"GLEAM-",vars_name,"_","C_",year,".nc")
             if os.path.exists(file_name_compress):
                   os.remove(file_name_compress)
             os.system('%s%s%s%s' % ('nccopy -d 9 ', file_name_write, ' ', file_name_compress))
             os.system('%s%s%s%s' % ('mv ', file_name_compress, ' ', file_name_write))
             del(file_name_write,file_name_compress,year) 
      del(files_v2A)
      files_v2B = glob.glob('%s%s%s' % (dir_v2B, vars_name,'?????.nc'))
      for f in files_v2B:
             print('%s%s' % ("Processing: ", f))
             # open file
             file_open_read = nc.Dataset(f,mode='r')
             # get var
             var_out = file_open_read.variables[vars_name][:,:,:]
             del(file_open_read)
             # work out the time, get year from string
             year = (f[len(f)-7:len(f)-3])
             # define CF_compliant time dim for each day of the year
             time_cf = get_cf_time_from_year(int(year),"leap")
             # write time,lat,lon, and var to file
             file_name_write = '%s%s%s%s%s%s' % (dir_write_v2B,"GLEAM-",vars_name,"_",year,".nc")
             if os.path.exists(file_name_write):
                   os.remove(file_name_write)
             # create file handle
             fileobj = nc.Dataset(file_name_write,mode='w')
             # write lat, lon, and time
             write_lat_lon_time(fileobj,latitude_v2B,longitude,time_cf)
             del(time_cf)
             # make nan to fill_value = -9
             var_out_no_nan = np.copy(var_out)
             var_out_no_nan[np.isnan(var_out)] = -9
             del(var_out)
             # write var
             write_var(fileobj,vars_name,var_out_no_nan)
             # close file
             fileobj.close()
             del(var_out_no_nan)
             # run compression
             file_name_compress = '%s%s%s%s%s%s%s' % (dir_write_v2B,"GLEAM-",vars_name,"_","C_",year,".nc")
             if os.path.exists(file_name_compress):
                   os.remove(file_name_compress)
             os.system('%s%s%s%s' % ('nccopy -d 9 ', file_name_write, ' ', file_name_compress))
             os.system('%s%s%s%s' % ('mv ', file_name_compress, ' ', file_name_write))
             del(file_name_write,file_name_compress,year)
      del(files_v2B)
