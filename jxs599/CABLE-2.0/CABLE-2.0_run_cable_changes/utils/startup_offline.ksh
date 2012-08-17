#!/bin/ksh 

#############################################################################
###  as well as a bit of book-keeping this script will compile CABLE      ###
###  from source code  and deposit the executable binary in               ###
###  this directory. CABLE is then run over each of the sites specified   ### 
###  in "sites_main.txt" and output data moved into a created directory   ###
###  labelled by the name of the site, exact location depends somewhat on ###
### 
###  plots of flux data is called.                                        ### 
#############################################################################



#******************************************************************************
#*************    USER SECTION   **********************************************
#******************************************************************************

# sitename: given to subdirectory in out/ per site 
site_name()
{
   set -A sites \
   Tumbarumba
#   bbb \
#   ccc
} 


# sitename: relative path to met. data 
site_data()
{
   set -A sitefile \
   data/sample_met/Bondville1997 
#   bbb \
#   ccc
}



#******************************************************************************
#*************    END USER SECTION   ******************************************
#******************************************************************************







#==============================================================================


check_nml()
{
   if [[ -e cable.nml ]]; then
      nml_exists=TRUE     
      print '\n\tThe existence of a namelist file takes precedence to'
      print'\tgenerating a new one based on configuration of this script.'    
      print '\n\tPress ENTER to proceed or Control-C to quit'
      read dummy
   fi
}


#==============================================================================


write_nml()
{
   print "&cable" > junk
   print "filename%met = ' {$sitefile[$1]} "
   print "filename%met = 'data/sample_met/Bondville1997_csiro.nc'"
}
#   print "filename%out = 'out_cable.nc'
#   print "filename%log = 'log_cable.txt'
#   print "filename%restart_in  = ' ' 
#   print "filename%restart_out = './restart_out.nc'
#   print "filename%type    = 'data/surface_data/gridinfo_CSIRO_1x1.nc'
#   print "filename%veg    = 'data/core/def_veg_params.txt'
#   print "filename%soil    = 'data/core/def_soil_params.txt'












#==============================================================================

run_cable()
{
#   integer kmax=${#sites[*]}
#   integer k=0
   
#   while [[ $k -lt $kmax ]]; do
#      if [[ $nml_exists = 'TRUE' ]]; then
#         print '\n\tUsing existing namelist.'
#      else
#         write_nml $k
#      fi   
      run_run #$k        
#      (( k = k + 1 ))
#   done 
}

run_run()
{
   # remove any trace of previous runs
   tidy
    
   # link to met data etc 
   if [[ ! -e data ]]; then
      ln -s ~/CABLE-AUX/data
   fi
   
   #mkdir out/${sites[$k]}

   # execute CABLE
   #print '\n*** RUNNING CABLE ***\n'
   #./cable 
}


#==============================================================================

         
book_keeping()
{      
   if [[ -d out.9 ]]; then
      print "\n\ntime to organize your out/ directory"
      exit
   fi 
   
   i=8; ((j=i+1)); k=0;
   while [ $k -lt 8 ]
      do 
         if [[ -d out.$i ]]; then
            mv out.$i out.$j
         fi 
         ((j = j - 1)) 
         ((i = i - 1)) 
         ((k = k + 1)) 
      done

   if [[ -d out ]]; then
      mv out out.1
   fi      
   mkdir out


}


#==============================================================================


tidy()
{
   rm -fr src/qsj.j src/*out qs* bu/  
   rm -f *.*out *.out *csv .qu 
}


#==============================================================================








######################################################################
###  set up directory for this run - make output directory and     ###
### possibly helpful book-keeping to allow for immediate execution ### 
### and avoid accidently over-writing data - up to a point !!      ###
######################################################################

# see above functions()

print 'executing script' 


#return
#exit
#check_nml

#site_name
#site_data

book_keeping

########################################################################
##### call CABLE.R in batch mode to avoid going into R first, and    ###
##### then clean up this directory (these files have already been    ###
##### dealt with in R-script )                                       ###
########################################################################

run_cable







