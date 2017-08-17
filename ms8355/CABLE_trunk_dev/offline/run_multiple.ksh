#!/bin/ksh 

#############################################################################
###  as well as a bit of book-keeping this script will compile CABLE      ###
###  from source code  and deposit the executable binary in               ###
###  this directory. CABLE will be executed for each fwsoil function      ###
###  possible configuration and +/- SDs depending on the user's choice    ###
###  CABLE is then run over each of the sites specified in "Sites.txt"    ###
###  and output data moved into a created directory labelled by the name  ###
###  of the site, exact location depends on plots of flux data called.    ### 
#############################################################################


## define functions
#==============================================================================


book_keeping()
{      
   if [[ -d $out.20 ]]; then
      print "\n\ntime to organize your out/ directory"
      exit
   fi 
   
   i=19; ((j=i+1)); k=0;
   while [ $k -lt 19 ]
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


# read sites file
# regenerate if number of sites does not match input folder
site_list_file()
{
    num_sites="$(grep -c 'Fluxnet.1.4_met.nc' $site_file)" # counts the sites in list
    num_site_files="$(ls $in | wc -l)" # counts site files in site folders
    if [[ "$num_sites" != "$num_site_files" ]] ; then
	grep '#' $site_file >> Temp_sites.txt
	find $in  -printf "%f\n" >> Temp_sites_2.txt # lists all input site files from folder
	sed '1d' Temp_sites_2.txt >> tmpfile # removes folder name from the list
	mv tmpfile Temp_sites_2.txt
	sed 's/Fluxnet.1.4_met.nc.*//' Temp_sites_2.txt >> Temp_sites_3.txt # creates list with site names only
	awk -v prefix="${in}/" '{print prefix $0}' Temp_sites_2.txt >> tmpfile # adds path to site files' names
	mv tmpfile Temp_sites_2.txt
	len_2="$(wc -l < Temp_sites_2.txt)" # counts lines in the newly created file
	len_3="$(wc -l < Temp_sites_3.txt)"
	if [[ "$len_2" != "$len_3" ]] ; then
	    echo "there is an issue with how you generate your site file"
	    exit
	fi
	if [[ "$len_2" == "$len_3" ]] ; then
	    for i in {1.."$len_2"} ; do
		extract_line="${i}p"
		sed -n $extract_line Temp_sites_3.txt >> Temp_sites_4.txt
		sed -n $extract_line Temp_sites_2.txt >> Temp_sites_4.txt # new text file has all names & paths
		echo -en '\n' >> Temp_sites_4.txt
	    done
	    sed -e '/#Sites/r Temp_sites_4.txt' Temp_sites.txt >> Temp_sites_5.txt # inserts names & paths after #Sites
	    num_sites_5="$(grep -c 'Fluxnet.1.4_met.nc' Temp_sites_5.txt)" # counts sites in the newly created file
	    if [[ "$num_sites" == "$num_sites_5" ]] ; then
		echo "there is an issue with how you generate your site file"
		exit
	    fi
	    if [[ "$num_sites" != "$num_sites_5" ]] ; then # this is satisfactory, the new site list is different from the previous one
		mv Temp_sites_5.txt Sites.txt
		site_file="Sites.txt"
		rm Temp_site*.txt
		echo "you have successfuly generated your new Site.txt file"
	    fi
	fi
    fi
}

#==============================================================================

# sitename given to subdirectory in out/ per site (should match cable.nml)
# this will be re-initiated properly for future release
site_name()
{
   integer i=0
   exec < ./${site_file}
   
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

   mkdir $out


   # execute CABLE
   if [[ ${fsites[$1]} != '' ]]; then
      ./create_cable-nml_casa_spin_and_run.sh -W $SD_loc/${sites[$1]}"_SD.txt" -w "$which_SD" -q "$drought_fun" -j $out/${sites[$1]}"_"   #create namelist with options above
      ./cable-r4334 ${fsites[$1]} ${fpoolsites[$1]}
   else
       echo "bypassing namelist arg"
      ./cable-r4334
   fi

   print '\n*** CABLE RUN FINISHED ***\n' 
   
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


get_SD()
{
   if [[ -d $output ]]; then # check if the simu output exists in this configuration
       for file in $output/*.nc
       do
	   cdo -selname,Fwsoil $file $SD_loc/out.nc # isolate fwsoil
	   cdo -timstd $SD_loc/out.nc $SD_loc/sd.nc # calculate sd over time
	   file=${file%_*}
	   file=${file%_*}
	   file=${file#../../Outputs/*/}
	   cdo -output $SD_loc/sd.nc >> $SD_loc/"$file"_SD.txt # write in text file
	   rm $SD_loc/out.nc $SD_loc/sd.nc
       done
   else # if the simu output doesn't exist, write SD=0
       for file in $in/*.nc
       do
	   file=${file%Fluxnet*}
	   file=${file#../Inputs/Test_sites/}
	   0.00 >> $SD_loc/"$file"_SD.txt
       done
   fi
}


#==============================================================================


write_SD()
{
   #check if directory exists, it not create
   if [[ -d $SD_loc ]]; then
       echo  "output directory exists for the SDs"
       rm $SD_loc/*.txt
   else
       mkdir -p $SD_loc
   fi

   output="../../Outputs/no_SD/standard"
   get_SD

   output="../../Outputs/no_SD/Haverd2013"
   get_SD

   output="../../Outputs/no_SD/LaiandKtaul2000"
   get_SD

   output="../../Outputs/no_SD/non-linearextrapolation"
   get_SD

   # add the same header to all files
   echo "SDs for fwsoil. Respectively: standard, Haverd2013, Lai and Katul 2000, non-linear extrapolation" >> ../../Outputs/headerfile

   for file in $SD_loc/*.txt
   do
       cat ../../Outputs/headerfile $file >> ../../Outputs/tmpfile2
       mv ../../Outputs/tmpfile2 $file
   done

   rm ../../Outputs/headerfile
}


#==============================================================================

## Cd run directory
cd . #/srv/ccrc/data04/z3509830/Fluxnet_data/Data_for_Manon/CABLE_trunk_dev/offline

# Set fwsoil option for creating output path
drought_fun="Haverd2013"
which_SD="minus_SD"

# paths
out="../../Outputs/""$which_SD""/""$drought_fun"
out="$(echo "$out" | sed 's/ //g')" #removes all white spaces from out path (there are white spaces with some of the drought_fun options)
in="../../Inputs/Test_sites"
SD_loc="../../Outputs/SDs"

# site file
site_file="Sites.txt"


######################################################################
###  set up directory for this run - make output directory and     ###
### possibly helpful book-keeping to allow for immediate execution ### 
### and avoid accidently over-writing data - up to a point !!      ###
######################################################################

#check if directory exists, it not create
if [[ -d $out ]]; then
   echo  "output directory exists for the runs"
else
    mkdir -p $out
fi

# see above functions()

book_keeping

site_list_file

site_name

########################################################################
##### call CABLE.R in batch mode to avoid going into R first, and    ###
##### then clean up this directory (these files have already been    ###
##### dealt with in R-script )                                       ###
########################################################################

run_cable

#write_SD
