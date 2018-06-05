#!/bin/bash

#PBS -m ae
#PBS -P dt6
#PBS -l walltime=27000
#PBS -l mem=7500MB
#PBS -l ncpus=1
#PBS -j oe
#PBS -q normal
#PBS -l wd
#PBS -l other=gdata1


module load netcdf/4.2.1.1
module load R

#  spin_casa.sh
#  
#
#  Created by Anna Ukkola on 1/07/2015.
#

##### SET SITE !!!!! #####

site="Konza" #"grassLTER"


## Spins CASA for determined number of cycles
## Outputs cnipool file

run_dir=`pwd`
data_dir="/g/data1/w35/amu561/Rainfall_assymmetry/"


### SET INPUT DIRECTORY
cd ${data_dir}/Inputs/Met_inputs/${site}


if [[ $run_dir =~ ([^,]+).*"CABLE-2.3.4-GW-Medlyn_spin"([^,]+) ]]; then
    INDIR=`find * -maxdepth 1 -type d -name "${site}${BASH_REMATCH[2]}"`
else

    #only for testing script
    INDIR=`find * -maxdepth 1 -type d -name "${site}200"`
fi




## Some options:
GWFLAG="TRUE"



cd $run_dir



# Create initial pool file #
cat > create_init_poolfile.R << EOF
library(raster)
library(ncdf4)
path <- "$data_dir"

dir.create(paste(path, "Inputs/CASA_ins", "$site", sep="/"))

#Read site lat/lon
site_file <- nc_open(paste(path, "/Inputs/Met_inputs/${site}/${site}100/${site}100_precip.nc", sep=""))
lat       <- ncvar_get(site_file, "latitude")
lon       <- ncvar_get(site_file, "longitude")
nc_close(site_file)

#Read in global pool file
pools <- read.csv(file=paste(path, "/Inputs/CASA_ins/cnppool1990_steady.csv", sep="/"), header=FALSE)
pools <- pools[,1:46]


colnames(pools) <- 
  c("ktau","npt","veg_iveg","soil_isoilm","casamet_isorder","casamet_lat","casamet_lon",
    "casamet_areacell","casamet_glai","casabiome_sla","phen_phase","casapool_clabile", 
    "casapool_cplant_leaf",  "casapool_cplant_wood", "casapool_cplant_root",
    "casapool_clitter_metb", "casapool_clitter_str", "casapool_clitter_cwd", 
    "casapool_csoil_mic",    "casapool_csoil_slow",  "casapool_csoil_pass",
    "casapool_nplant_leaf",  "casapool_nplant_wood", "casapool_nplant_root",
    "casapool_nlitter_metb", "casapool_nlitter_str", "casapool_nlitter_cwd",
    "casapool_nsoil_mic",    "casapool_nsoil_slow",  "casapool_nsoil_pass",
    "casapool_nsoilmin",
    "casapool_pplant_leaf",  "casapool_pplant_wood", "casapool_pplant_root",  
    "casapool_plitter_metb", "casapool_plitter_str", "casapool_plitter_cdw",
    "casapool_psoil_mic",    "casapool_psoil_slow",  "casapool_psoil_pass",
    "casapool_psoillab","casapool_psoilsorb","casapool_psoilocc","casabal_sumcbal","casabal_sumnbal","casabal_sumpbal")


#Read in gridinfo file for dimensions
gridinfo <- raster("$data_dir/Inputs/gridinfo_CSIRO_1x1_default.nc", varname="iveg")


data_brick <- brick(nrow=nrow(gridinfo), ncol=ncol(gridinfo), xmn=xmin(gridinfo), xmx=xmax(gridinfo), ymn=ymin(gridinfo), ymx=ymax(gridinfo), 
                    nl=(ncol(pools)), crs='+proj=longlat +datum=WGS84')

#Columns to convert to raster
lat_lon <- which(colnames(pools)=="casamet_lat" | colnames(pools)=="casamet_lon")

for(k in 1:nlayers(data_brick))
{
  data_raster <- raster(nrow=nrow(data_brick), ncol=ncol(data_brick))
  extent(data_raster) <- extent(data_brick)
  Cells <- cellFromXY(data_raster, pools[,rev(lat_lon)])
  data_raster[Cells] <- pools[,k]
  
  data_brick[[k]] <- data_raster
}


#Extract site values
site_values <- extract(data_brick, cbind(lon,lat))

#Replace lat/lon values with actual site coords
site_values[1,which(grepl("casamet_lon", colnames(site_values)))] <- lon
site_values[1,which(grepl("casamet_lat", colnames(site_values)))] <- lat

#Write to file
write.table(x=site_values, file=paste(path, "/Inputs/CASA_ins/${site}/poolcnp_init_${site}.csv", sep=""), 
          col.names=FALSE, row.names=FALSE, sep=",   ", fileEncoding= "UTF-8")

EOF


Rscript create_init_poolfile.R  #source R script to create init pool file



#Loop through experiments
for D in $INDIR
do

echo "Experiment:"
echo $D


## name of out directory

OUTDIR="${D}_outs"

if [[ -d "${data_dir}/Outputs/${site}/${OUTDIR}" ]]; then
    echo "output directory exists"
else
    mkdir -p ${data_dir}/Outputs/${site}/${OUTDIR}
fi



#Site met file
metfile="Met_inputs/${site}/${D}/CABLE_met_input_${D}.nc"




########## Spin model ###########

#number of iterations
#maxI=2
#for ((I=1;I<=$maxI;I++)); #I in 1 2 3
#do

#echo "PRINTING index:"
#echo $I


#Set tolerance for spin-up
tol=$(echo "scale=6; 0.05" | bc)

#Initialise csoil and allocation variables for spin-up check
diff_cplant=$(echo "scale=6; 99999.0" | bc)
diff_csoil=$(echo "scale=6; 99999.0" | bc)

I=1




######## Run model if tolerance not achieved ########
while [[ $(echo "$diff_cplant > $tol" | bc -l) -gt 0 || $(echo "$diff_csoil > $tol" | bc -l) -gt 0 ]]
do


echo "PRINTING index:"
echo $I


## Set CNPpool file name and restart in and out filenames
restart_out="restart_out_${I}.nc"
if [[ $I -eq 1 ]]; then
    pool_file="poolcnp_init_${site}.csv"
    restart_in=""
else
    restart_in="restart_out_${[I-1]}.nc"
    pool_file="poolcnpOut_old_${D}.csv"
fi





### Create namelist file ###

cat > cable.nml << EOF
&cable
filename%met = "${data_dir}/Inputs/$metfile"
filename%out = "${data_dir}/Outputs/${site}/${OUTDIR}/out_cable.nc"
filename%log = "${data_dir}/Outputs/${site}/${OUTDIR}/log_cable.txt"
filename%restart_in  = "${data_dir}/Outputs/${site}/${OUTDIR}/${restart_in}"
filename%restart_out = "${data_dir}/Outputs/${site}/${OUTDIR}/${restart_out}"
filename%type    = "${data_dir}/Inputs/gridinfo_CSIRO_1x1.nc"
filename%veg    = "${data_dir}/Inputs/def_veg_params_Ticket2_with_medlyn.txt"
filename%soil    = "${data_dir}/Inputs/def_soil_params.txt"
filename%gw_elev = ""
vegparmnew = .TRUE.  ! using new format when true
soilparmnew = .TRUE.  ! using new format when true
spinup = .FALSE.  ! do we spin up the model?
delsoilM = 0.001   ! allowed variation in soil moisture for spin up
delsoilT = 0.01    ! allowed variation in soil temperature for spin up
output%restart = .TRUE.  ! should a restart file be created?
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
check%ranges     = .TRUE.  ! variable ranges, input and output
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
casafile%cnpipool = "${data_dir}/Inputs/CASA_ins/${site}/${pool_file}"
casafile%cnpbiome = "${data_dir}/Inputs/CASA_ins/pftlookup_csiro_v16_17tiles_Ticket2.csv"
casafile%cnpepool = "${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv"   ! end of run pool size
casafile%cnpmetout= ""  ! output daily met forcing for spinning casacnp
casafile%cnpmetin = ""  ! list of daily met files for spinning casacnp
casafile%phen     = "${data_dir}/Inputs/CASA_ins/modis_phenology_csiro.txt"
casafile%cnpflux  ="${data_dir}/Outputs/${site}/${OUTDIR}/cnpfluxOut_${D}.csv"
output%CASA = .FALSE.   ! output CASA variables to netcdf?
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
! 3. Lai and Katul 2000
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
cable_user%L_newProfile = .TRUE. !Use new YP profile for root water extraction?
cable_user%or_evap = .FALSE.     !Use new Or soilE formulation?
&end

EOF





#run model
executable=`ls cable-r*`   #automatically find current cable version
./$executable




#Calculate difference in csoil and cplant between current and previous timestep


## C-soil ##
if [[ $I -eq 1 ]]
then
    csoil_old=99999
else
    cmic_old=`awk -F "\"*,\"*" '{print $19}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    cslow_old=`awk -F "\"*,\"*" '{print $20}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    cpass_old=`awk -F "\"*,\"*" '{print $21}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    csoil_old=$(echo "$cmic_old+$cslow_old+$cpass_old" | bc -l)
fi

cmic_new=`awk -F "\"*,\"*" '{print $19}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
cslow_new=`awk -F "\"*,\"*" '{print $20}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
cpass_new=`awk -F "\"*,\"*" '{print $21}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
csoil_new=$(echo "$cmic_new+$cslow_new+$cpass_new" | bc -l)


## C-plant ##
if [[ $I -eq 1 ]]
then
    cplant_old=99999
else
    cleaf_old=`awk -F "\"*,\"*" '{print $13}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    cwood_old=`awk -F "\"*,\"*" '{print $14}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    croot_old=`awk -F "\"*,\"*" '{print $15}' ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv`
    cplant_old=$(echo "$cleaf_old+$cwood_old+$croot_old" | bc -l)
fi

cleaf_new=`awk -F "\"*,\"*" '{print $13}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
cwood_new=`awk -F "\"*,\"*" '{print $14}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
croot_new=`awk -F "\"*,\"*" '{print $15}' ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv`
cplant_new=$(echo "$cleaf_new+$cwood_new+$croot_new" | bc -l)


echo "PRINTING old/new csoil"
echo $csoil_old
echo $csoil_new


diff_cplant=$(echo "$cplant_new-$cplant_old" | bc -l)
diff_csoil=$(echo "$csoil_new-$csoil_old" | bc -l)


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
cp -rf ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}_prev.csv
cp -rf ${data_dir}/Inputs/CASA_ins/${site}/cnpfluxOut_old_${D}.csv ${data_dir}/Inputs/CASA_ins/${site}/cnpfluxOut_old_${D}_prev.csv
fi

cp -rf ${data_dir}/Outputs/${site}/${OUTDIR}/poolcnpOut_${D}.csv ${data_dir}/Inputs/CASA_ins/${site}/poolcnpOut_old_${D}.csv
cp -rf ${data_dir}/Outputs/${site}/${OUTDIR}/cnpfluxOut_${D}.csv ${data_dir}/Inputs/CASA_ins/${site}/cnpfluxOut_old_${D}.csv


I=$[I+1]



done # while loop, casa iterations

done # experiments $D



