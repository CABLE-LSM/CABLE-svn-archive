#!/bin/ksh 
set -x
file="ftemp"
i=0
while read line[$i]
do
   # display line or do somthing on $line
	echo "$line" 
   (( i = i +1 ))
done <"$file"

typeset -L ${line[0]}
typeset -R ${line[0]}

   print '\n*** CABLE RUN FINISHED ***\n'
   
   # crude test for successful run, move files
   if [[ -f out_cable.nc ]]; then
      
      mkdir out/${line[0]}

      print '\n*** CABLE RUN (appears) SUCCESSFULL ***\n'
      		
      # CABLE output + restart if applicable
      mv log_cable.txt out_cable.nc restart_out.nc out/${line[0]}
      
      # CABLE namelist file used 
      cp cable.nml out/${line[0]}
      
      # pools for CASA-CNP
      mv poolcnpOut.csv cnpfluxOut.csv out/${line[0]}
   
   else
      print '\n*** ERROR: RUN FAILED ***\n'     
   fi  
   
   # clean up
   rm -f fort.66

