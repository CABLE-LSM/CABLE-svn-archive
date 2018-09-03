#!/bin/bash

# Run LIS
#PBS -m ae
#PBS -P dt6
#PBS -l walltime=170200
#PBS -l mem=7500MB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l other=gdata1


module load intel-mpi/4.1.1.036
module load netcdf/4.2.1.1
module load R



### Runs through three steps of CASA spin-up ###


#############################
###  Set options for run  ###
#############################

### Set file paths ###

site="Konza"   #"Konza", "grassLTER"

run_dir=${PWD##*/}    #`pwd`

data_dir="/srv/ccrc/data45/z3509830/CABLE_runs//Rainfall_assymmetry/"

if [[ $run_dir =~ "CABLE-2.3.4-GW-Medlyn_spin" ]]; then
    INDIR="${site}${run_dir#*GW-Medlyn_spin}"   #creates string of site name and experiment (e.g. Konza100)
else
    INDIR=$site
fi


cable_spinpath="${data_dir}/Outputs/${site}/${INDIR}_spin_newparams"   #where store CABLE outputs from spin-up?
cable_runpath="${data_dir}/Outputs/${site}/${INDIR}_run_newparams"   #where store CABLE outputs from final run (if step4=true)?
casa_spinpath=$cable_outpath    #where store CASA outputs?


#Create symbolic links to input and outputs directories
#Avoids issues with CABLE character length limitations

#create output folders if don't exist
if [[ ! -d $cable_spinpath ]]; then
mkdir -p $cable_spinpath
fi
if [[ ! -d $cable_runpath ]]; then
mkdir -p $cable_runpath
fi

ln -s ${data_dir}/Inputs Inputs
ln -s $cable_spinpath Outputs_spin
ln -s $cable_runpath Outputs_run

cable_outpath="./Outputs_spin"
casa_outpath="./Outputs_spin"
cable_finalpath="./Outputs_run"


met_file="/Met_inputs/${site}/${INDIR}/CABLE_met_input_${INDIR}.nc"

casa_biome="pftlookup_csiro_v16_17tiles_Ticket2_${site}.csv"


#CO2_file="./Annual_CO2_concentration_until_2010.txt"  #CO2 concentration file if want annually varying
#CO2_con=370

### Run settings ###

LAIflag="TRUE"      #use prognostic LAI?
Vcmax_flag="FALSE"  #use prognostic vcmax?
icycle=1            #which icycle (1=C, 2=CN, 3=CNP)



executable="cable-r4182" #name of CABLE executable

#Set spin-up options
#Set step 2 to true to use fast spinup. Otherwise recycles a normal run.
step1=true #run Step 1 (NPP stabilisation)?
step2=false #run Step 2 (fast spin-up)?
step3=true #run Step 3 (further spin-up and final run)?


#create output folders if don't exist
#if [[ ! -d cable_outpath ]]; then
#mkdir -p $cable_outpath
#fi
#if [[ ! -d $casa_outpath ]]; then
#mkdir -p $casa_outpath
#fi




#############################
### Step 1: Stabilise NPP ###
#############################


if [[ "$step1" = true ]]; then


    iter=1  #loop counter


    while [[ $iter -lt 200 ]]   #while iter smaller than 20, keep going (using 20 so doesn't run all day...)
    do



        echo "Spinning up Step #1, iteration: ${iter} ----------------------------------###"


        ## Set restart_in filenames ##

        if [[ $iter -eq 1 ]]; then
            #first iteration, no input restart files and CASA started from zero
            #Give a random name to restart in files, code will fail if fields left empty
            restart_cable_in="none"
            restart_casa_in="none"
            casa_zero="TRUE"
        else
            #else use restart from previous year
            restart_cable_in="${cable_outpath}/cable_restart_Step1_$(($iter-1)).nc"
            restart_casa_in="${casa_outpath}/casa_restart_Step1_$(($iter-1)).nc"
            casa_zero="FALSE"

            #Check that restart files exist, if not abort code (something's not working?)
            if [[ ! -e "$restart_cable_in"  || ! -e "$restart_casa_in" ]]; then
                echo "Cannot find CABLE or CASA restart file, code not working? ABORTING."
            exit
            fi

        fi

        restart_cable_out="cable_restart_Step1_${iter}.nc"
        restart_casa_out="casa_restart_Step1_${iter}.nc"


        #Read annual CO2 concentration
        #CO2_con=`sed -n "${year}p" < $CO2_file`
        #CO2_con=${CO2_con%.*}
        #echo "Annual CO2 concentration: ${CO2_con}"


        #Create namelist
        #Word of warning: if any input fields are empty, namelist won't be written correctly
        ./create_cable-nml_casa_spin_and_run.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $casa_zero -e false -b false -y $iter -i $restart_cable_in -r $restart_cable_out -a $restart_casa_in -f $restart_casa_out -k $casa_biome

        #Run CABLE
        ./$executable

        #save namelist every year
        cp cable.nml $cable_outpath/cable_namelist_step1_iter_${iter}.nml


        #Check if NPP has reaches equilibrium (varies < 1%)
        if [[ $iter -gt 1 ]]; then


            #Create R script for checking global mean NPP
            cat > npp_check.R << EOF
            #Check if packages exist
            if(!require(raster)) install.packages("raster")
            if(!require(ncdf4)) install.packages("ncdf4")

            #Load packages
            library(raster)
            library(ncdf4)

            #Find outputs for current and previous years
            file_prev <- list.files(path="$casa_outpath", pattern=paste('casa_restart_Step1_$(($iter-1)).nc'), full.names=TRUE)    #pattern must match cable output file pattern
            file_current <- list.files(path="$casa_outpath", pattern=paste('casa_restart_Step1_${iter}.nc'), full.names=TRUE)

            if(length(file_prev)==0 | length(file_current)==0)
            {
                print("Looking for output file in $cable_outpath")
                stop("Cannot find output files Rscript Step1 !!!")
            }

            data <- lapply(list(file_prev,file_current), function(x) round(ncvar_get(nc_open(x), "cplant"), digits=8))


            leaf_diff <- abs((data[[2]][1]-data[[1]][1])/data[[2]][1])
            wood_diff <- abs((data[[2]][2]-data[[1]][2])/data[[2]][2])
            root_diff <- abs((data[[2]][3]-data[[1]][3])/data[[2]][3])

	        difference <- sum(leaf_diff, wood_diff, root_diff, na.rm=TRUE)


	        status <- paste("Leaf current: ", data[[2]][1], ", Leaf previous: ", data[[1]][1], ", Leaf difference: ", leaf_diff,
                            "	Wood current: ", data[[2]][2], ", Wood previous: ", data[[1]][2], ", Wood difference: ", wood_diff,
                            "	Root current: ", data[[2]][3], ", Root previous: ", data[[1]][3], ", Root difference: ", root_diff, sep=" ")
    

            print(status)


            #If difference < threshold, create final restart file
            if (difference < 0.01) {
                file.copy(from='${casa_outpath}/casa_restart_Step1_${iter}.nc', to='${casa_outpath}/CASA_restart_Step1_final.nc')
                file.copy(from='${cable_outpath}/cable_restart_Step1_${iter}.nc', to='${cable_outpath}/CABLE_restart_Step1_final.nc')
            } else {
                write.csv(status, file="NPP_difference_iter_${iter}.csv")
            }
    
EOF

            #not pretty indenting but EOF statement has to be in first column...

            #Run R-script
            Rscript npp_check.R

        fi


        #If r-script creates final restart files, exit loop and continue to Step 2
        if [[ -e "${casa_outpath}/CASA_restart_Step1_final.nc" && "${cable_outpath}/CABLE_restart_Step1_final.nc" ]]; then
            echo "Step 1 finished, plant C pools vary < 1%"
            break
        fi


        #If NPP hasn't stabilised, check if R wrote text file (if it did not, problem in code, exit)
        if [[ ! -e "NPP_difference_iter_${iter}.csv" && $iter -gt 1 ]]; then
            echo "Could not find R NPP difference file, problem in code? ABORTING."
            exit
        fi


        iter=$(($iter+1))

    done #while


    echo "END Step1 of spin-up ----------------------------------###"


fi #Step1





############################
### Step 2: fast spin-up ###
############################


#Loops through years once to create daily met, soil moisture and soil carbon outputs for each year for Step 3

if [[ "$step2" = true ]]; then

    ###################################
    ### 1st step: Create met inputs ###
    ###################################

    echo "Spinning up Step #2 met input ----------------------------------###"


    #Set step1 final restart file as starting point
    restart_cable_in="${cable_outpath}/CABLE_restart_Step1_final.nc"
    restart_casa_in="${cable_outpath}/CASA_restart_Step1_final.nc"

    restart_cable_out="CABLE_restart_Step2_met.nc"
    restart_casa_out="CASA_restart_Step2_met.nc"


    #Check if finds restart file, if not abort code (something's not working?)
    if [[ ! -e "$restart_cable_in" || ! -e "$restart_casa_in" ]]; then
        echo "Cannot find restart file in Step 2, code not working? ABORTING."
        exit
    fi


    #Read annual CO2 concentration
    #CO2_con=`sed -n "${year}p" < $CO2_file`
    #CO2_con=${CO2_con%.*}
    #echo "Annual CO2 concentration: ${CO2_con}"


    #Create namelist
    ./create_cable-nml_casa_spin_and_run.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -e TRUE -i $restart_cable_in -t true -r $restart_cable_out -a $restart_casa_in -f $restart_casa_out -k $casa_biome #-n $CO2_con

    #Run CABLE
    ./$executable

    #save namelist
    cp cable.nml $cable_outpath/cable_namelist_step2_first_stage.nml




    echo "END Step2 of spin-up met input ----------------------------------###"


    ##################################
    ### 2nd step: Run fast spin-up ###
    ##################################

    echo "Spinning up Step #2 fast spin-up ----------------------------------###"


    # Spins up CASA until number of iterations equals mloop (set in cable_mpimaster.F90 for an mpi run or cable_driver.F90 for serial run)

    restart_cable_in=$restart_cable_out
    restart_casa_in=$restart_casa_out

    restart_cable_out="CABLE_restart_Step2_final.nc"
    restart_casa_out="CASA_restart_Step2_final.nc"


    #List output met files from Step 2

#    files=`find ${casa_outpath}/casa_dump_cnpspin*.nc`
#    N=`find ${casa_outpath}/casa_dump_cnpspin*.nc | wc -l`


    #have to be non-indented or cable can't read the created file correctly !!!!!!
#cat > ${casa_outpath}/fcnpspin.lst << EOF
#$N
#$files
#EOF


    #Read CO2 concentration
    #CO2_con=`sed -n "${year}p" < $CO2_file`
    #CO2_con=${CO2_con%.*}
    #echo "Annual CO2 concentration: ${CO2_con}"


    ./create_cable-nml_casa_spin_and_run.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -e false -u true -i $restart_cable_in -r $restart_cable_out -a $restart_casa_in -f $restart_casa_out -k $casa_biome #-n $CO2_con

    ./$executable

    #save namelist
    cp cable.nml $cable_outpath/cable_namelist_step2_second_stage.nml


    #Check that wrote a restart file, if not abort code (something's not working?)
    if [[ ! -e "${cable_outpath}/${restart_out}" ]]; then
        echo "Cannot find restart file at the end of Step 2. ABORTING."
        exit
    fi


    echo "END Step2 of fast spin-up ----------------------------------###"


fi





###############################
### Step 3: Full simulation ###
###############################


#Determine which restart file to use depending on whether
#using fast spin-up
if [[ "$step2" = true ]]; then
     restart_cable_init="${cable_outpath}/CABLE_restart_Step2_final.nc"
     restart_casa_init="${casa_outpath}/CASA_restart_Step2_final.nc"
else
     restart_cable_init="${cable_outpath}/CABLE_restart_Step1_final.nc"
     restart_casa_init="${casa_outpath}/CASA_restart_Step1_final.nc"
fi



if [[ "$step3" = true ]]; then

    #Create output directory
    if [[ ! -d cable_finalpath ]]; then
        mkdir -p $cable_finalpath
    fi


    #############################################
    ### 1st step: stabilise soil carbon pools ###
    #############################################

    #loop counter, initialise if step1 false
    if [[ "$step1" = false ]]; then
        iter=1      
    else 
        iter=$(($iter+1))
    fi

    #save first iter
    first_iter=$iter


    while [[ $iter -lt 2000 ]]   #while iter smaller than 1000, keep going (number picked randomly)
    do

        echo "Spinning up Step #3, iteration: ${iter} ----------------------------------###"


        #Set restart_in filename
        if [[ $iter -eq $first_iter ]]; then
            #first iteration, use initial file
            restart_cable_in=$restart_cable_init
            restart_casa_in=$restart_casa_init
        else
            #else use restart from previous year
            restart_cable_in="${cable_outpath}/cable_restart_Step3_$(($iter-1)).nc"
            restart_casa_in="${casa_outpath}/casa_restart_Step3_$(($iter-1)).nc"
        fi

        restart_cable_out="cable_restart_Step3_${iter}.nc"
        restart_casa_out="casa_restart_Step3_${iter}.nc"



        #Check if finds restart file, if not abort code (something's not working?)
        if [[ ! -e "$restart_cable_in" || ! -e "$restart_casa_in" ]]; then
            echo "Cannot find input restart file for Step 3 spin-up, code not working? ABORTING."
            exit
        fi



        #Create namelist
        ./create_cable-nml_casa_spin_and_run.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -e false -b true -i $restart_cable_in -r $restart_cable_out -a $restart_casa_in -f $restart_casa_out -k $casa_biome  #-n $CO2_con
	

        #Run CABLE
        ./$executable

        #save namelist every year
        cp cable.nml $cable_outpath/cable_namelist_step3_iter_${iter}.nml



        #Check if NPP has reaches equilibrium (varies < 1%)
	    #Use first_iter+1 so gets a chance to write a cable_output.nc file in case
	    #one doesn't already exist (used to determine run length)
        if [[ $iter -gt $first_iter ]]; then


            #Create R script for checking global mean NPP
            cat > csoil_check.R << EOF
            #Check if packages exist
            if(!require(raster)) install.packages("raster")
            if(!require(ncdf4)) install.packages("ncdf4")

            #Load packages
            library(raster)
            library(ncdf4)

            #Find outputs for current and previous years
            file_prev <- list.files(path="$casa_outpath", pattern=paste('casa_restart_Step3_$(($iter-1)).nc'), full.names=TRUE)    #pattern must match cable output file pattern
            file_current <- list.files(path="$casa_outpath", pattern=paste('casa_restart_Step3_${iter}.nc'), full.names=TRUE)

            file_time <- list.files(path="$cable_outpath", pattern=paste('cable_output_0.nc'), full.names=TRUE)

            if(length(file_prev)==0 | length(file_current)==0)
            {
                print("Looking for output file in $cable_outpath")
                stop("Cannot find output files Rscript Step3 !!!")
            }

            data <- lapply(list(file_prev,file_current), function(x) ncvar_get(nc_open(x), "csoil"))

            #Setting this manually now, needs FIXING !!!!!
            tsteps_per_day <- 48


            time_info <- ncvar_get(nc_open(file_time), "SoilCarbPassive")
	    yrs <- round(length(time_info)/365/tsteps_per_day, digits=0)


            #If difference <1%, create final restart file

            difference <- abs(data[[2]][3]-data[[1]][3]) / yrs

	        status <- paste("Csoil passive current: ", data[[2]][3], ", Csoil passive previous: ", data[[1]][3], ", Csoil passive difference: ",
			    difference, sep=" ")
 
            print(status)

		
	        #Use 0.5 gC/m2 per 5 years as the threshold (i.e. 0.1 gC/m2 per yr)

            if (difference < 0.1) {
                file.copy(from='${cable_outpath}/cable_restart_Step3_${iter}.nc', to='${cable_outpath}/CABLE_restart_Step3_final.nc')
                file.copy(from='${casa_outpath}/casa_restart_Step3_${iter}.nc', to='${casa_outpath}/CASA_restart_Step3_final.nc')
            } else {
                write.csv(status, file="Csoil_passive_difference_Step3_iter_${iter}.csv")
            }

EOF

            #not pretty indenting but EOF statement has to be on first column...

            #Run R-script
            Rscript csoil_check.R

        fi


        #If r-script creates restart file, exit loop and continue to Step 2
        if [[ -e "${cable_outpath}/CABLE_restart_Step3_final.nc" ]]; then
            echo "Step 3 spin-up finished, Csoil varies < 0.1%"
            break
        fi

        #If NPP hasn't stabilised, check if R wrote text file (if did not, probably problem in code, abort)
        if [[ ! -e "Csoil_passive_difference_Step3_iter_${iter}.csv" && $iter -gt $first_iter ]]; then
            echo "Could not find R Csoil difference file, problem in code? ABORT"
            exit
        fi


        iter=$(($iter+1))

    done #while


    echo "END Step3 of spin-up ----------------------------------###"



    ###########################
    ### 2nd step: Final run ###
    ###########################

    echo "Performing final model run Step #3 ----------------------------------###"

    restart_cable_in="${cable_outpath}/CABLE_restart_Step3_final.nc"
    restart_casa_in="${casa_outpath}/CASA_restart_Step3_final.nc"


    #Create name list (use spin-up TRUE to spin up soil moisture; also make sure to write CASA outputs to file)
    ./create_cable-nml_casa_spin_and_run.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_finalpath -e true -i $restart_cable_in -a $restart_casa_in -b true -z true -k $casa_biome #-n $CO2_con


    ./$executable


    echo "END Step3 final run FINISHED ----------------------------------###"

fi


