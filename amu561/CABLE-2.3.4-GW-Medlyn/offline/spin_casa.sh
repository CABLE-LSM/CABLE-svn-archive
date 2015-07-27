#!/bin/bash

#  spin_casa.sh
#  
#
#  Created by Anna Ukkola on 1/07/2015.
#



## Spins CASA for determined number of cycles
## Outputs cnipool file


### SET INPUT DIRECTORY
cd /srv/ccrc/data45/z3509830/CABLE_runs/Rainfall_assymmetry/Inputs/Met_inputs/
INDIR=`find * -maxdepth 1 -type d -name "STU100*"`

## Some options:
GWFLAG="TRUE"



cd ./../../CABLE-2.3.4-GW-Medlyn/offline/


#Loop through experiments
for D in $INDIR
do

echo "Experiment:"
echo $D


## name of out directory

OUTDIR="${D}_outs"

if [[ -d "./../../Outputs/"$OUTDIR ]]; then
    echo "output directory exists"
else
    mkdir ./../../Outputs/$OUTDIR
fi



#Site met file
metfile="Met_inputs/${D}/CABLE_met_input_${D}.nc"




########## Spin model ###########

#number of iterations
#maxI=2
#for ((I=1;I<=$maxI;I++)); #I in 1 2 3
#do

#echo "PRINTING index:"
#echo $I


#Set tolerance for spin-up
tol=$(echo "scale=6; 0.005" | bc)

#Initialise csoil and allocation variables for spin-up check
diff_cplant=$(echo "scale=6; 99999.0" | bc)
diff_csoil=$(echo "scale=6; 99999.0" | bc)

I=1




######## Run model if tolerance not achieved ########
while [[ $(echo "$diff_cplant > $tol" | bc -l) || $(echo "$diff_csoil > $tol" | bc -l) ]]
do


echo "PRINTING index:"
echo $I


## Set CNPpool file name
if [[ $I -eq 1 ]]; then
    pool_file="poolcnpInStubai.csv"
else
    pool_file="poolcnpOut_old_${D}.csv"
fi




### Create namelist file ###

cat > cable.nml << EOF
&cable
filename%met = "./../../Inputs/$metfile"
filename%out = "./../../Outputs/${OUTDIR}/out_cable.nc"
filename%log = "./../../Outputs/${OUTDIR}/log_cable.txt"
filename%restart_in  = ""
filename%restart_out = "./../../Outputs/${OUTDIR}/restart_out.nc"
filename%type    = "./../../Inputs/gridinfo_CSIRO_1x1.nc"
filename%veg    = "./../../Inputs/def_veg_params_Ticket2_with_medlyn.txt"
filename%soil    = "./../../Inputs/def_soil_params.txt"
filename%gw_elev = "./../../Inputs/$metfile"
vegparmnew = .TRUE.  ! using new format when true
soilparmnew = .TRUE.  ! using new format when true
spinup = .FALSE.  ! do we spin up the model?
delsoilM = 0.001   ! allowed variation in soil moisture for spin up
delsoilT = 0.01    ! allowed variation in soil temperature for spin up
output%restart = .FALSE.  ! should a restart file be created?
output%met = .FALSE.  ! input met data
output%flux = .FALSE.  ! convective, runoff, NEE
output%soil = .FALSE.  ! soil states
output%snow = .FALSE.  ! snow states
output%radiation = .FALSE.  ! net rad, albedo
output%carbon    = .FALSE.  ! NEE, GPP, NPP, stores
output%veg       = .FALSE.  ! vegetation states
output%params    = .FALSE.  ! input parameters used to produce run
output%balances  = .FALSE.  ! energy and water balances
output%averaging = "all"
check%ranges     = .FALSE.  ! variable ranges, input and output
check%energy_bal = .FALSE.  ! energy balance
check%mass_bal   = .FALSE.  ! water/mass balance
verbose = .FALSE. ! write details of every grid cell init and params to log?
leaps = .false. ! calculate timing with leap years?
logn = 88      ! log file number - declared in input module
fixedCO2 = 390.0   ! if not found in met file, in ppmv
spincasainput = .FALSE.    ! input required to spin casacnp offline
spincasa      = .FALSE.     ! spin casa before running the model if TRUE, and should be set to FALSE
l_casacnp     = .TRUE.  ! using casaCNP with CABLE
l_laiFeedbk   = .TRUE.  ! using prognostic LAI
l_vcmaxFeedbk = .FALSE.  ! using prognostic Vcmax
icycle = 1   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
casafile%cnpipool = "./../../Inputs/CASA_ins/$pool_file"
casafile%cnpbiome = "./../../Inputs/CASA_ins/pftlookup_csiro_v16_17tiles_Ticket2.csv"
casafile%cnpepool = "./../../Outputs/${OUTDIR}/poolcnpOut_${D}.csv"   ! end of run pool size
casafile%cnpmetout= "./../../Outputs/${OUTDIR}/casa_met_out.nc"  ! output daily met forcing for spinning casacnp
casafile%cnpmetin = ''          ! list of daily met files for spinning casacnp
casafile%phen     = "./../../Inputs/CASA_ins/modis_phenology_csiro.txt"
casafile%cnpflux  ="./../../Outputs/${OUTDIR}/cnpfluxOut_${D}.csv"
ncciy = 0  ! 0 for not using gswp; 4-digit year input for year of gswp met
gswpfile%rainf = "gswp/AustMerra3h${year}.nc"
gswpfile%snowf = "gswp/AustMerra3h${year}.nc"
gswpfile%LWdown= "gswp/AustMerra3h${year}.nc"
gswpfile%SWdown= "gswp/AustMerra3h${year}.nc"
gswpfile%PSurf = "gswp/AustMerra3h${year}.nc"
gswpfile%Qair  = "gswp/AustMerra3h${year}.nc"
gswpfile%Tair  = "gswp/AustMerra3h${year}.nc"
gswpfile%wind  = "gswp/AustMerra3h${year}.nc"
redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
wiltParam = 0.5
satuParam = 0.8
cable_user%GS_SWITCH = 'medlyn'             ! choices are:
! 1. leuning
! 2. medlyn
cable_user%FWSOIL_SWITCH = 'standard'        ! choices are:
! 1. standard
! 2. nonlinear
! 3. Lai and Ktaul 2000
cable_user%DIAG_SOIL_RESP = 'ON '
cable_user%LEAF_RESPIRATION = 'ON '
cable_user%RUN_DIAG_LEVEL= 'BASIC'        ! choices are:
! 1. BASIC
! 1. NONE
cable_user%CONSISTENCY_CHECK= .TRUE.      ! TRUE outputs combined fluxes at each timestep for comparisson to A control run
cable_user%CASA_DUMP_READ = .FALSE.      ! TRUE reads CASA forcing from netcdf format
cable_user%CASA_DUMP_WRITE = .FALSE.      ! TRUE outputs CASA forcing in netcdf format
cable_user%SSNOW_POTEV= 'HDM'      ! Humidity Deficit Method
cable_user%GW_MODEL = .${GWFLAG}.       !True means use the groundwater module, false means use default soil snow scheme
cable_user%alt_forcing = .FALSE.
cable_user%GSWP3       = .FALSE.
gw_params%EfoldMaxSatFrac = 6.0
gw_params%MaxHorzDrainRate = 0.01
gw_params%MaxSatFraction = 0.0  !Changed Anna
gw_params%EfoldHorzDrainRate = 1.0
&end

EOF





#run model
./cable





#Calculate difference in csoil and cplant between current and previous timestep


csoil_old=`awk -F "\"*,\"*" '{print $15}' ../../Inputs/CASA_ins/${pool_file}`
csoil_new=`awk -F "\"*,\"*" '{print $15}' ./../../Outputs/${OUTDIR}/poolcnpOut_${D}.csv`


if [[ $I -eq 1 ]]
then
    cplant_old=99999
else
    cleaf_old=`awk -F "\"*,\"*" '{print $11}' ../../Inputs/CASA_ins/cnpfluxOut_old_${D}.csv`
    cwood_old=`awk -F "\"*,\"*" '{print $12}' ../../Inputs/CASA_ins/cnpfluxOut_old_${D}.csv`
    croot_old=`awk -F "\"*,\"*" '{print $13}' ../../Inputs/CASA_ins/cnpfluxOut_old_${D}.csv`
    cplant_old=$cleaf_old+$cwood_old+$croot_old
fi

cleaf_new=`awk -F "\"*,\"*" '{print $11}' ../../Outputs/${OUTDIR}/cnpfluxOut_${D}.csv`
cwood_new=`awk -F "\"*,\"*" '{print $12}' ../../Outputs/${OUTDIR}/cnpfluxOut_${D}.csv`
croot_new=`awk -F "\"*,\"*" '{print $13}' ../../Outputs/${OUTDIR}/cnpfluxOut_${D}.csv`
cplant_new=$(echo "$cleaf_new+$cwood_new+$croot_new" | bc -l)


diff_cplant=$(echo "$cplant_new-$cplant_old" | bc -l)
diff_csoil=$(echo "$csoil_new-$csoil_old" | bc -l )


diff_csoil=$(echo "scale=6; a=$diff_csoil/1;if(0>a)a*=-1;a"|bc)
diff_cplant=$(echo "scale=6; a=$diff_cplant/1;if(0>a)a*=-1;a"|bc)



echo "PLANT C difference:"
echo $diff_cplant
echo "SOIL C difference:"
echo $diff_csoil




## Copy current pool file to casa inputs for next step
## Also rename file from previous iteration so it can be used for checking values
if [[ $I -gt 1 ]]
then
cp -rf ./../../Inputs/CASA_ins/poolcnpOut_old_${D}.csv ./../../Inputs/CASA_ins/poolcnpOut_old_${D}_prev.csv
cp -rf ./../../Inputs/CASA_ins/cnpfluxOut_old_${D}.csv ./../../Inputs/CASA_ins/cnpfluxOut_old_${D}_prev.csv
fi

cp -rf ./../../Outputs/${OUTDIR}/poolcnpOut_${D}.csv ./../../Inputs/CASA_ins/poolcnpOut_old_${D}.csv
cp -rf ./../../Outputs/${OUTDIR}/cnpfluxOut_${D}.csv ./../../Inputs/CASA_ins/cnpfluxOut_old_${D}.csv


I=$[I+1]



done # while loop, casa iterations

done # experiments $D



