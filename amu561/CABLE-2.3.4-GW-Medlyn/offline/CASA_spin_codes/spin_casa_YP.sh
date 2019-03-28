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

site="Tumbarumba"   #"Konza", "grassLTER"

run_dir=${PWD##*/}    #`pwd`

data_dir="./"


if [[ $run_dir =~ "CABLE-2.3.4-GW-Medlyn_spin" ]]; then
    INDIR="${site}${run_dir#*GW-Medlyn_spin}"   #creates string of site name and experiment (e.g. Konza100)
else
    INDIR=$site
fi


cable_outpath="${data_dir}/Outputs_spin"   #where store CABLE outputs from spin-up?
cable_finalpath="${data_dir}/Outputs_run_new"   #where store CABLE outputs from final run (if step4=true)?
casa_outpath=$cable_outpath    #where store CASA outputs?


met_file="TumbaFluxnet.1.4_met.nc"

casa_biome="pftlookup_csiro_v16_17tiles_Ticket2.csv"


#CO2_file="./Annual_CO2_concentration_until_2010.txt"  #CO2 concentration file if want annually varying
#CO2_con=370

### Run settings ###

LAIflag="TRUE"		#use prognostic LAI?
Vcmax_flag="FALSE"	#use prognostic vcmax?
icycle=1		#which icycle (1=C, 2=CN, 3=CNP)

GWflag="FALSE"		#use groundwater scheme?



executable="cable" #name of CABLE executable

#Set spin-up options
#Set step 2 to true to use fast spinup. Otherwise recycles a normal run.
step0=true #Create initial pool file?
step1=true #run Step 1 (NPP stabilisation)?
step2=true #run Step 2 (fast spin-up)?
step3=true #run a full simulation with spun-up model?


#Create output folders if don't exist
if [[ ! -d cable_outpath ]]; then
mkdir -p $cable_outpath
fi
if [[ ! -d $casa_outpath ]]; then
mkdir -p $casa_outpath
fi



################################
### Create initial pool file ###
################################


#Create namelist file
restart_init="Initial_restart_file.nc"
./create_cable-nml_casaspin_Step1-3.sh -m $met_file -p $casa_outpath -j $cable_outpath -g $GWflag -r $restart_init -k $casa_biome #-n $CO2_con


#File to add CASA variables to (used to initialise Step1)
init_restart_file="${casa_outpath}/CABLE_restart_with_CASA_vars_${INDIR}.nc"  #restart file where to find initial pool sizes


if [[ "$step0" = true ]]; then


	#Run CABLE
	./$executable


	cat > create_initial_poolfile.R << EOF

	    if(!require(raster)) install.packages("raster")
	    if(!require(ncdf4)) install.packages("ncdf4")

	    library(raster)
	    library(ncdf4)

	    rm( list=ls(all=TRUE) )

	    #Read in global pool file
	    pools <- read.csv(file="${data_dir}/cnppool1990_steady.csv", header=FALSE)
	    pools <- pools[,1:46]

	    colnames(pools) <-
	    c("ktau","npt","veg_iveg","soil_isoilm","casamet_isorder","casamet_lat","casamet_lon", "casamet_areacell",
	    "casamet_glai",                                                          #NEEDED: LAI, read from surface data file
	    "casabiome_sla",
	    "phen_phase",            "casapool_clabile",                             #NEEDED: Phase, Clabile
	    "casapool_cplant_leaf",  "casapool_cplant_wood", "casapool_cplant_root", #NEEDED: CASA_Cplant
	    "casapool_clitter_metb", "casapool_clitter_str", "casapool_clitter_cwd", #NEEDED: Clitter
	    "casapool_csoil_mic",    "casapool_csoil_slow",  "casapool_csoil_pass",  #NEEDED: CASA_Csoil
	    "casapool_nplant_leaf",  "casapool_nplant_wood", "casapool_nplant_root", #NEEDED: Nplant
	    "casapool_nlitter_metb", "casapool_nlitter_str", "casapool_nlitter_cwd", #NEEDED: Nlitter
	    "casapool_nsoil_mic",    "casapool_nsoil_slow",  "casapool_nsoil_pass",  #NEEDED: Nsoil
	    "casapool_nsoilmin",                                                     #NEEDED: Nsoilmin
	    "casapool_pplant_leaf",  "casapool_pplant_wood", "casapool_pplant_root", #NEEDED: Pplant
	    "casapool_plitter_metb", "casapool_plitter_str", "casapool_plitter_cdw", #NEEDED: Plitter
	    "casapool_psoil_mic",    "casapool_psoil_slow",  "casapool_psoil_pass",  #NEEDED: Psoil
	    "casapool_psoillab",     "casapool_psoilsorb",   "casapool_psoilocc",    #NEEDED: Psoillab, Psoilsorb, Psoilocc
	    "casabal_sumcbal","casabal_sumnbal","casabal_sumpbal")



	    ###--- Convert pool text file into raster ---###

	    #Read in gridinfo file for dimensions (And use January LAI as starting point)
	    gridinfo <- raster("${data_dir}/CABLE-AUX/offline/gridinfo_CSIRO_1x1.nc", varname="LAI")


 	    data_brick <- brick(nrow=nrow(gridinfo), ncol=ncol(gridinfo), xmn=xmin(gridinfo), xmx=xmax(gridinfo), ymn=ymin(gridinfo), ymx=ymax(gridinfo),
 	    nl=(ncol(pools)), crs='+proj=longlat +datum=WGS84')

 	    #Convert pool data to raster brick
 	    lat_lon <- which(colnames(pools)=="casamet_lat" | colnames(pools)=="casamet_lon")

 	    for(k in 1:nlayers(data_brick))
 	    {
	        data_raster <- raster(nrow=nrow(data_brick), ncol=ncol(data_brick))
	        extent(data_raster) <- extent(data_brick)
	        Cells <- cellFromXY(data_raster, pools[,rev(lat_lon)])
	        data_raster[Cells] <- pools[,k]
	        data_brick[[k]] <- data_raster

	    }

	    #Set layer names
	    names(data_brick) <- colnames(pools)


	    ###--- Extract and write data to file ---###

	    # CABLE variables to add to restart_file
	    vars <- c("LAI", "phase", "Clabile", "CASA_Cplant", "Clitter", "CASA_Csoil", "Nplant", "Nlitter", "Nsoil",
	    "Nsoilmin", "Pplant", "Plitter", "Psoil", "Psoillab", "Psoilsorb", "Psoilocc")

	    #Extract site coordinates from met file
	    met_nc <- nc_open("$met_file")
	    lon    <- ncvar_get(met_nc, "longitude")
	    lat    <- ncvar_get(met_nc, "latitude")

	    extract_data <- function(raster)
	    {
	        data_vec <- extract(raster, matrix(data=c(lon,lat), ncol=2, nrow=1))
	        return(data_vec)
 	    }


	    #Extract land data points for each variable  (Use backslash so bash code ignores dollar signs (for R use only in this case)
	    var_data <- list(LAI   = extract_data(gridinfo),
	    phase     = extract_data(data_brick\$phen_phase),
	    clabile   = extract_data(data_brick\$casapool_clabile),
	    cplant    = rbind(extract_data(data_brick\$casapool_cplant_leaf), extract_data(data_brick\$casapool_cplant_wood), extract_data(data_brick\$casapool_cplant_root)),
 	    clitter   = rbind(extract_data(data_brick\$casapool_clitter_metb), extract_data(data_brick\$casapool_clitter_str), extract_data(data_brick\$casapool_clitter_cwd)),
 	    csoil     = rbind(extract_data(data_brick\$casapool_csoil_mic), extract_data(data_brick\$casapool_csoil_slow), extract_data(data_brick\$casapool_csoil_pass)),
	    nplant    = rbind(extract_data(data_brick\$casapool_nplant_leaf), extract_data(data_brick\$casapool_nplant_wood), extract_data(data_brick\$casapool_nplant_root)),
	    nlitter   = rbind(extract_data(data_brick\$casapool_nlitter_metb), extract_data(data_brick\$casapool_nlitter_str), extract_data(data_brick\$casapool_nlitter_cwd)),
	    nsoil     = rbind(extract_data(data_brick\$casapool_nsoil_mic), extract_data(data_brick\$casapool_nsoil_slow), extract_data(data_brick\$casapool_nsoil_pass)),
	    nsoilmin  = extract_data(data_brick\$casapool_nsoilmin),
	    pplant    = rbind(extract_data(data_brick\$casapool_pplant_leaf), extract_data(data_brick\$casapool_pplant_wood), extract_data(data_brick\$casapool_pplant_root)),
	    plitter   = rbind(extract_data(data_brick\$casapool_plitter_metb), extract_data(data_brick\$casapool_plitter_str), extract_data(data_brick\$casapool_plitter_cdw)),
	    psoil     = rbind(extract_data(data_brick\$casapool_psoil_mic), extract_data(data_brick\$casapool_psoil_slow), extract_data(data_brick\$casapool_psoil_pass)),
	    psoillab  = extract_data(data_brick\$casapool_psoillab),
	    psoilsorb = extract_data(data_brick\$casapool_psoilsorb),
	    psoilocc  = extract_data(data_brick\$casapool_psoilocc)
	    )


	    dims <- list()

	    #Define extra dimension for plant/soil/litter variables
	    dimplant  <- ncdim_def(name="pools_plant",units="", vals=1:3)
	    dimlitter <- ncdim_def(name="pools_litter",units="", vals=1:3)
	    dimsoil   <- ncdim_def(name="pools_soil",units="", vals=1:3)

	    #Name of modified restart file
	    file.remove("$init_restart_file")  #remove if already exists
	    file.copy(from = "${cable_outpath}/Initial_restart_file.nc", to = "$init_restart_file")

	    #Open restart file
	    ncid = nc_open("$init_restart_file", write=TRUE, readunlim=FALSE)

	    for(k in 1:length(var_data))
	    {

 	       # Define dimensions:
 	       if(length(var_data[[k]]) == length(ncid\$dim\$mp_patch\$vals)) {
        	    #mp only
           	 dims <- list(ncid\$dim\$mp_patch)
 	       } else if(nrow(var_data[[k]]) == 3 & grepl(pattern="plant", names(var_data)[k]))  {
	            #mp and plant
 	           dims <- list(ncid\$dim\$mp_patch, dimplant)
	        } else if(nrow(var_data[[k]]) == 3 & grepl(pattern="litter", names(var_data)[k])) {
	            #mp and litter
	            dims <- list(ncid\$dim\$mp_patch, dimlitter)
	        } else if(nrow(var_data[[k]]) == 3 & grepl(pattern="soil", names(var_data)[k]))  {
	            #mp and soil
	            dims <- list(ncid\$dim\$mp_patch, dimsoil)
 	   }


  	   # Define variable:
 	   hcvar = ncvar_def(name=vars[k], units='', dim=dims,
  	 		     missval=-9999,longname=vars[k])

  	   # Add variable and then variable data:
 	   ncid = ncvar_add(ncid,hcvar) # CRITICAL - use new netcdf file handle
 	   ncvar_put(ncid, vars[k], var_data[[k]])

 	   rm(hcvar)
 	}


	nc_close(ncid)

EOF

	#Run R script to create initial restart file
	Rscript ./create_initial_poolfile.R

fi



#############################
### Step 1: Stabilise NPP ###
#############################

#lai        #icycle    #casa outpath    #cable out path   #use GW?   #CO2 conc.  #spinup  #out casa netcdf
#./create_cable-nml_casaspin.sh -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $GWflag -n $CO2_con -e false -b true


if [[ "$step1" = true ]]; then


    iter=1  #loop counter


    while [[ $iter -lt 200 ]]   #while iter smaller than 20, keep going (using 20 so doesn't run all day...)
    do



        echo "Spinning up Step #1, iteration: ${iter} ----------------------------------###"


        #Set restart_in filename
        if [[ $iter -eq 1 ]]; then
            #first iteration, use initial file
            restart_in=$init_restart_file
        else
            #else use restart from previous year
            restart_in="${cable_outpath}/cable_restart_Step1_$(($iter-1)).nc"
        fi

        restart_out="cable_restart_Step1_${iter}.nc"



        #Check if finds restart file, if not abort code (something's not working?)
        if [[ ! -e "$restart_in" ]]; then
            echo "Cannot find restart file, code not working? ABORTING."
            exit
        fi



        #Read annual CO2 concentration
        #CO2_con=`sed -n "${year}p" < $CO2_file`
        #CO2_con=${CO2_con%.*}
        #echo "Annual CO2 concentration: ${CO2_con}"


        #Create namelist
        ./create_cable-nml_casaspin_Step1-3.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $GWflag -e false -b false -y $iter -i $restart_in -r $restart_out -k $casa_biome  #-n $CO2_con

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
            file_prev <- list.files(path="$cable_outpath", pattern=paste('cable_restart_Step1_$(($iter-1)).nc'), full.names=TRUE)    #pattern must match cable output file pattern
            file_current <- list.files(path="$cable_outpath", pattern=paste('cable_restart_Step1_${iter}.nc'), full.names=TRUE)

            if(length(file_prev)==0 | length(file_current)==0)
            {
                print("Looking for output file in $cable_outpath")
                stop("Cannot find output files Rscript Step1 !!!")
            }

            data <- lapply(list(file_prev,file_current), function(x) ncvar_get(nc_open(x), "CASA_Cplant"))

            #If difference <1%, create final restart file

            leaf_diff <- abs((data[[2]][1]-data[[1]][1])/data[[2]][1])
            wood_diff <- abs((data[[2]][2]-data[[1]][2])/data[[2]][2])
            root_diff <- abs((data[[2]][3]-data[[1]][3])/data[[2]][3])

	    difference <- sum(leaf_diff, wood_diff, root_diff, na.rm=TRUE)


            status <- paste("Difference in NPP:", round(difference*100, digits=2), "%, iter: ", $iter)

            print(status)

            if (difference < 0.01) {
                file.copy(from='${cable_outpath}/cable_restart_Step1_${iter}.nc', to='${cable_outpath}/cable_restart_Step1_final.nc')
            } else {
                write.csv(status, file="NPP_difference_iter_${iter}.csv")
            }
    
EOF

            #not pretty indenting but EOF statement has to be in first column...

            #Run R-script
            Rscript npp_check.R

        fi


        #If r-script creates restart file, exit loop and continue to Step 2
        if [[ -e "${cable_outpath}/cable_restart_Step1_final.nc" ]]; then
            echo "Step 1 finished, NPP varies < 1%"
            break
        fi


        #If NPP hasn't stabilised, check if R wrote text file (if did not, probably problem in code, abort)
        if [[ ! -e "NPP_difference_iter_${iter}.csv" && $iter -gt 1 ]]; then
            echo "Could not find R NPP difference file, problem in code? ABORT"
            exit
        fi


        iter=$(($iter+1))

    done #while


    echo "END Step1 of spin-up ----------------------------------###"


fi #Step1





############################
### Step 2: fast spin-up ###
############################


#./create_cable-nml_casaspin.sh -d $LAIflag -c $icycle -p $casa_outpath -g $GWflag -n $CO2_con
#-e true  # spin up (I think??)
#-t true  # set casaspinin option to true


#Loops through years once to create daily met, soil moisture and soil carbon outputs for each year for Step 3

if [[ "$step2" = true ]]; then

    ###################################
    ### 1st step: Create met inputs ###
    ###################################

    echo "Spinning up Step #2 met input ----------------------------------###"


    #Set step1 final restart file as starting point
    restart_in="${cable_outpath}/cable_restart_Step1_final.nc"

    restart_out="cable_restart_Step2.nc"


    #Check if finds restart file, if not abort code (something's not working?)
    if [[ ! -e "$restart_in" ]]; then
        echo "Cannot find restart file in Step 2, code not working? ABORTING."
        exit
    fi


    #Read annual CO2 concentration
    #CO2_con=`sed -n "${year}p" < $CO2_file`
    #CO2_con=${CO2_con%.*}
    #echo "Annual CO2 concentration: ${CO2_con}"


    #Create namelist										    #should spin==true?				     #casaspin in truei      -normal spinup true
    ./create_cable-nml_casaspin_Step1-3.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $GWflag -e TRUE -b FALSE -i $restart_in -t true -r $restart_out -k $casa_biome #-n $CO2_con

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

    restart_in="${cable_outpath}/cable_restart_Step2.nc"
    restart_out="CABLE_restart_final_Step2.nc"


    if [[ ! -d cable_finalpath ]]; then
        mkdir -p $cable_finalpath
    fi

    #List output met files from Step 2

    files=`find ${casa_outpath}/casa_dump_cnpspin*.nc`
    N=`find ${casa_outpath}/casa_dump_cnpspin*.nc | wc -l`


    #have to be non-intended or cable can't read the created file correctly !!!!!!
cat > ${casa_outpath}/fcnpspin.lst << EOF
$N
$files
EOF


    #Read CO2 concentration
    #CO2_con=`sed -n "${year}p" < $CO2_file`
    #CO2_con=${CO2_con%.*}
    #echo "Annual CO2 concentration: ${CO2_con}"


    ./create_cable-nml_casaspin_Step4.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $GWflag -e false -u true -i $restart_in -r $restart_out -k $casa_biome #-n $CO2_con

    ./$executable

    #save namelist
    cp cable.nml $cable_outpath/cable_namelist_step2_second_stage.nml


    echo "END Step2 of fast spin-up ----------------------------------###"


fi




###############################
### Step 3: Full simulation ###
###############################


#Determine which restart file to use depending on whether
#using fast spin-up
if [[ "$step2" = true ]]; then
     restart_init="${cable_outpath}/CABLE_restart_final_Step2.nc"
else
     restart_init="${cable_outpath}/cable_restart_Step1_final.nc"
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
            restart_in=$restart_init
        else
            #else use restart from previous year
            restart_in="${cable_outpath}/cable_restart_Step3_$(($iter-1)).nc"
        fi

        restart_out="cable_restart_Step3_${iter}.nc"



        #Check if finds restart file, if not abort code (something's not working?)
        if [[ ! -e "$restart_in" ]]; then
            echo "Cannot find restart file Step 3 spin-up, code not working? ABORTING."
            exit
        fi



        #Create namelist
            ./create_cable-nml_casaspin_Step1-3.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_outpath -g $GWflag -e false -b true -i $restart_in -r $restart_out -k $casa_biome  #-n $CO2_con
	

        #Run CABLE
        ./$executable

        #save namelist every year
        cp cable.nml $cable_outpath/cable_namelist_step3_iter_${iter}.nml



        #Check if NPP has reaches equilibrium (varies < 1%)
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
            file_prev <- list.files(path="$cable_outpath", pattern=paste('cable_restart_Step3_$(($iter-1)).nc'), full.names=TRUE)    #pattern must match cable output file pattern
            file_current <- list.files(path="$cable_outpath", pattern=paste('cable_restart_Step3_${iter}.nc'), full.names=TRUE)

	    file_time <- list.files(path="$cable_outpath", pattern=paste('cable_output_0.nc'), full.names=TRUE)

            if(length(file_prev)==0 | length(file_current)==0)
            {
                print("Looking for output file in $cable_outpath")
                stop("Cannot find output files Rscript Step3 !!!")
            }

            data <- lapply(list(file_prev,file_current), function(x) ncvar_get(nc_open(x), "CASA_Csoil"))

	    time_info <- ncvar_get(nc_open(file_time), "CASAGPP")
	    yrs <- round(length(time_info)/365, digits=0)


            #If difference <1%, create final restart file

            difference <- abs(data[[2]][3]-data[[1]][3]) / yrs

            status <- paste("Difference in passive Csoil:", round(difference*100, digits=4), "%, iter: ", $iter)
   
            print(status)

		
	    #Use 0.5 gC/m2 per 50 years as the threshold (i.e. 0.01 gC/m2 per yr)	

            if (difference < 0.01) {
                file.copy(from='${cable_outpath}/cable_restart_Step3_${iter}.nc', to='${cable_outpath}/cable_restart_Step3_final.nc')
            } else {
            status_all <- paste0(status, ", old Csoil: ", round(data[[1]][[3]], digits=5), ", new Csoil: ", round(data[[2]][[3]], digits=5))
                write.csv(status_all, file="Csoil_passive_difference_Step3_iter_${iter}.csv")
            }

EOF

            #not pretty indenting but EOF statement has to be on first column...

            #Run R-script
            Rscript csoil_check.R

        fi


        #If r-script creates restart file, exit loop and continue to Step 2
        if [[ -e "${cable_outpath}/cable_restart_Step3_final.nc" ]]; then
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

    restart_in="${cable_outpath}/CABLE_restart_Step3_final.nc"


    #Create name list (use spin-up TRUE to spin up soil moisture; also make sure to write CASA outputs to file)
    ./create_cable-nml_casaspin_Step4.sh -m $met_file -d $LAIflag -c $icycle -p $casa_outpath -j $cable_finalpath -g $GWflag -e true -i $restart_in -b true -z true -k $casa_biome #-n $CO2_con


    ./$executable


    echo "END Step3 final run FINISHED ----------------------------------###"

fi


