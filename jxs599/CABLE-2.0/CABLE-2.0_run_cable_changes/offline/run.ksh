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


# sitename given to subdirectory in out/ per site (should match cable.nml)
# this will be re-initiated properly for future release
site_name()
{
   integer i=0
   exec < sites.txt
   
   print "\n\tRunning CABLE over the single sites: \n" 
   while read line
   do
   	echo "$line"
      (( i = i + 1 ))
   done 
   integer isites=i/2
   
   i=0
   exec < sites.txt
   while ((i < isites)) 
   do
      read sites[i]	
      read fsites[i]	
      (( i = i + 1 ))
   done 

} 


#==============================================================================


run_cable()
{
   integer kmax=${#sites[*]}
   integer k=0
   
   while [[ $k -lt $kmax ]]; do
      run_run $k        
      (( k = k + 1 ))
   done 
}

run_run()
{
   # remove any trace of previous runs
   tidy

   # link to met data etc 
   if [[ ! -e data ]]; then
      ln -s ~/CABLE-AUX/data
   fi
   
   mkdir out/${sites[$k]}

   # execute CABLE
   print '\n*** RUNNING CABLE ***\n'
   ./cable ${fsites[$k]}
   print '\n*** CABLE RUN FINISHED ***\n'
   
   # crude test for successful run, move files
   if [[ -f out_cable.nc ]]; then

      print '\n*** CABLE RUN (appears) SUCCESSFULL ***\n'
      		
      # CABLE output + restart if applicable
      mv log_cable.txt out_cable.nc restart_out.nc out/${sites[$k]}
      
      # CABLE namelist file used 
      cp cable.nml out/${sites[$k]}
      
      # pools for CASA-CNP
      mv poolcnpOut.csv cnpfluxOut.csv out/${sites[$k]}
   
   else
      print '\n*** ERROR: RUN FAILED ***\n'     
   fi  
   
   # clean up
   rm -f fort.66
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

book_keeping

site_name

########################################################################
##### call CABLE.R in batch mode to avoid going into R first, and    ###
##### then clean up this directory (these files have already been    ###
##### dealt with in R-script )                                       ###
########################################################################

run_cable







