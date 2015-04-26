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



## Set input paths here
cd /srv/ccrc/data45/z3509830/CABLE_runs/Inputs/
DIRS=`find * -maxdepth 1 -type d -name "PLUMBER_sites_slope_0*"`
cd ../CABLE_2.0_new_hydrology_plus_MedlynStom/offline


#Set gc_opt option (Leuning/Medlyn) and GW option (TRUE/FALSE)
gc_opt="Leuning"

gw="FALSE"
gw_flag="no_gw"

drought_fun="standard"


#!!!!!! SET ALL OPTIONS IN create_cable-nml.sh     !!!!!




#########################################################################

#Loop through input dirs
for D in $DIRS
do


echo $D

out_options=`echo $D | cut -d'_' -f 3-8`       #extracts slope, isoil and lai info from path
out='../../Outputs/CABLE_2.0_new_hydrology_plus_MedlynStom/${gc_opt}_${drought_fun}_${gw_flag}_${out_options}'

#model options 
gw_opt=$gw   #use groundwater?
spin_opt="TRUE" #spin up model?
veg_file="g1_files/def_veg_params_medlyn_mean.txt"

site_file=$D".txt"


#check if directory exists, it not create
if [[ -d $out ]]; then
   echo  "output directory exists"
else
    mkdir -p $out
fi



book_keeping()
{      
   if [[ -d $out.9 ]]; then
      print "\n\ntime to organize your out/ directory"
      exit
   fi 
   
   i=8; ((j=i+1)); k=0;
   while [ $k -lt 8 ]
      do 
         if [[ -d $out.$i ]]; then
            mv $out.$i $out.$j
         fi 
         ((j = j - 1)) 
         ((i = i - 1)) 
         ((k = k + 1)) 
      done

   if [[ -d $out ]]; then
      mv $out $out.1
   fi      
   mkdir $out

   HOST_MACH=`uname -n | cut -c 1-4`
   
   if [[ $HOST_MACH = 'vayu' ]]; then
   
      if [[ ! -e cable.nml ]]; then
         ln -s $CABLE_AUX/CABLE-AUX/offline/cable.nml .
      fi
      if [[ ! -e sites.txt ]]; then
         ln -s $CABLE_AUX/CABLE-AUX/offline/sites.txt .
      fi
      
   fi
}


#==============================================================================


# sitename given to subdirectory in out/ per site (should match cable.nml)
# this will be re-initiated properly for future release
site_name()
{
   integer i=0
   exec < ./../../Inputs/Site_files/${site_file}
   
    echo $site_file

   while read line
   do
   	fchar=`echo "$line" | cut -c 1`  
      if [[ $fchar != '#' ]]; then
         echo $line >> fsites.txt     
         (( i = i + 1 ))
      fi   
   done 
   integer isites=i/3
   
   i=0
   exec < fsites.txt
   while ((i < isites)) 
   do
      read sites[i]	
      read fsites[i]	
      read fpoolsites[i]	
      print "\n\t${sites[i]}"
      print "\n\t${fsites[i]}"
      (( i = i + 1 ))
  done 

   rm -f fsites.txt
} 


#==============================================================================


run_cable()
{
   integer kmax=${#sites[*]}
   integer k=0
   
   #  work around to desired trigerring from cable.nml  
   if [[ $kmax = 0 ]]; then
      kmax=1
      sites[0]='default'
   fi
              
   while [[ $k -lt $kmax ]]; do
      run_run $k        
      (( k = k + 1 ))
   done 
}


#==============================================================================


run_run()
{
   # remove any trace of previous runs
   tidy

   mkdir $out/${sites[$1]}


   # execute CABLE
   if [[ ${fsites[$1]} != '' ]]; then
      ./create_cable-nml.sh -e $spin_opt -g $gw_opt -v $veg_file  #create namelist with options above
      ./cable ${fsites[$1]} ${fpoolsites[$1]}
   else
       echo "bypassing namelist arg"
      ./cable       
   fi

   print '\n*** CABLE RUN FINISHED ***\n'
   
   # crude test for successful run, move files
   if [[ -f out_cable.nc ]]; then

      print '\n*** CABLE RUN (appears) SUCCESSFULL ***\n'
      		
      # CABLE output + restart if applicable
      mv log_cable.txt out_cable.nc restart_out.nc $out/${sites[$1]}
      cp cable.nml  $out/${sites[$1]}
      # pools for CASA-CNP
      if [[ -e poolcnpOut.csv ]]; then
         mv poolcnpOut.csv cnpfluxOut.csv out/${sites[$1]}
      fi 
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


done     #Done looping through input DIRS




