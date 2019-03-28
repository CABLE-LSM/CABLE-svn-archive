#!/bin/bash
# Argument = -t test -r server -p password -v

usage()
{
cat << EOF
usage: $0 options

This script setups up cable.nml and can run the model

OPTIONS:
   -h      Show this message
   -m      metfile.  not needed if year /= 0 (gswp run)
   -j      CABLE output path
   -o      Output file
   -l      log file 
   -i      restart file to be used as input (netcdf)
   -r      restart file output at end of run
   -v      vegetation data file
   -s      soil data file
   -y      year to run the model (gswp forcing)
   -e      true -> spin to equilibrium
   -g      Use groundwater scheme
   -a      averaging flag - daily,monthly,yearly
   -f      alternate (merra) forcing flag
   -c	   icycle option
   -p      CASA output path
   -d      true -> use LAI feedback
   -x      true -> use vxmax feedback
   -t      true -> spincasain option
   -u      true -> spincasa option
   -n      CO2 concentration
   -b      true -> write CASA output to netcdf
   -k      CASA biome parameter file
EOF
}

#set some defaults here?  or do with if;fi after getopts

while getopts "hm:j:o:l:i:r:v:s:y:e:g:a:f:c:p:d:x:t:u:n:b:k:" OPTION

do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         m)
             metfile=$OPTARG
             ;;
         j)
	     cableout_path=$OPTARG
	     ;;
         o)
             outfile=$OPTARG
             ;;
         l)
             logfile=$OPTARG
             ;;
         i)
             restart_infile=$OPTARG
             ;;
         r)
             restart_outfile=$OPTARG
             ;;
         v)
             vegfile=$OPTARG
             ;;
         s)
             soilfile=$OPTARG
             ;;
         y)
             year=$OPTARG
             ;;
         e)
             eqflag=$OPTARG
             ;;
         g)
             GWflag=$OPTARG
             ;;
         a)
             AVGflag=$OPTARG
             ;;
         f)
            ALTFORCEflag=$OPTARG
            ;;
	 c)
     	    icycle_flag=$OPTARG
	    ;;
	 p)
	    casaout_path=$OPTARG
	    ;;
	 d) 
	    LAIflag=$OPTARG
  	    ;;
	 x)
	    Vcmax_flag=$OPTARG
	    ;;
	 t)
	    spincasain=$OPTARG
 	    ;;
	 u)
	    spincasa=$OPTARG
	    ;;
	 n) 
	    CO2_var=$OPTARG
	    ;;
	 b)
	    casaout_flag=$OPTARG
	    ;;
	 k) 
            casabiome=$OPTARG
	    ;;
         ?)
             usage
             exit
             ;;
     esac
done

if [ -z $eqflag ]; then   #default to no spinup
  eqflag="false"
fi

if [ -z $eqflag ]; then   #default to no groundwater
  GWflag="false"
fi

if [ -z $year ]; then
  year=0
  echo "No year supplied. Setting ncicty=0"
fi

if [ -z $metfile ] && [ $year -eq 0 ]; then
  metfile="default_met.txt"
  echo "Using default met file: $metfile Hopefully this is correct as ncity=0"
fi

if [ -z $cableout_path ]; then
  cableout_path="./"
fi

if [ -z $outfile ]; then
  outfile="cable_output_${year}.nc"          #if gswp append year to log
fi

if [ -z $logfile ]; then
  logfile='cable_log.txt'
  if [ $year -ne 0 ]; then
      logfile="cable_log_${year}.txt"          #if gswp append year to log
  fi
  echo "Using default log file $logfile"
fi

if [ -z $restart_infile ]; then
  if [ $year -gt 1900 ]; then              #gswp run after first year
    let pyear=$year-1
    restart_infile="restart_${pyear}.nc"
  fi
  echo "Using default restart_in file $restart_infile"
fi

if [ -z $restart_outfile ]; then
  if [ $year -eq 0 ]; then
    restart_outfile="restart_out.nc"   #default for non-gswp run
  else
    restart_outfile="cable_restart_${year}.nc"
  fi
  echo "Using default restart_out file $restart_outfile"
fi

if [ -z $vegfile ]; then
  vegfile='def_veg_params_Ticket2.txt'
  echo "Using default veg file $vegfile"
fi

if [ -z $soilfile ]; then
  soilfile='def_soil_params.txt'
  echo "Using default soil file $soilfile"
fi

if [ -z $AVGflag ]; then
  AVGflag='all'
  echo "Setting the output averaging flag to all"
fi

if [ -z $ALTFORCEflag ]; then   #default to no spinup
  ALTFORCEflag="false"
fi

if [ -z $icycle_flag ]; then   #icycle
  icycle_flag=0
fi

if [ -z $casaout_path ]; then   #CASA output path
  casaout_path="./" 
fi

if [ -z $LAIflag ]; then     #Run CASA with LAI feedback?
  LAIflag="FALSE"
fi

if [ -z $Vcmax_flag ]; then  #Run CASA with Vcmax feedback?
  Vcmax_flag="FALSE"
fi

if [ -z $spincasain ]; then  #true for step2, output casa met?
  spincasain="FALSE"
fi

if [ -z $spincasa ]; then    #true for step 3, spin casa?
  spincasa="FALSE"
fi

if [ -z $CO2_var ]; then   #CO2 concentration
  CO2_var="370"
fi

if [ -z $casaout_flag ]; then
  casaout_flag="TRUE"
fi

if [ -z $casabiome ];then
  casabiome="pftlookup_csiro_v16_17tiles_Ticket2.csv"
fi

#Set CASA-CNP flag (should CASA be used or not?)
if [[ $icycle_flag -eq 0 ]]; then
    l_casacnp="FALSE"
else
    l_casacnp="TRUE"
fi


in_path="./CABLE-AUX/"


touch $(pwd)/cable.nml

cat > $(pwd)/cable.nml << EOF
&cable
   filename%met = "$metfile"
   filename%out = "${cableout_path}/${outfile}"
   filename%log = "${cableout_path}/${logfile}"
   filename%restart_in  = "$restart_infile" 
   filename%restart_out = "${cableout_path}/${restart_outfile}"
   filename%type    = "${in_path}/offline/gridinfo_CSIRO_1x1.nc"
   filename%veg    = "${in_path}/core/biogeophys/$vegfile"
   filename%soil    = "${in_path}/core/biogeophys/${soilfile}"
   vegparmnew = .TRUE.  ! using new format when true
   soilparmnew = .TRUE.  ! using new format when true
   spinup = .${eqflag}.  ! do we spin up the model?
   delsoilM = 0.01   ! allowed variation in soil moisture for spin up
   delsoilT = 0.1    ! allowed variation in soil temperature for spin up
   output%restart = .TRUE.  ! should a restart file be created?
   output%met = .TRUE.  ! input met data
   output%flux = .TRUE.  ! convective, runoff, NEE
   output%soil = .TRUE.  ! soil states
   output%snow = .TRUE.  ! snow states
   output%radiation = .TRUE.  ! net rad, albedo
   output%carbon    = .TRUE.  ! NEE, GPP, NPP, stores
   output%veg       = .TRUE.  ! vegetation states
   output%params    = .TRUE.  ! input parameters used to produce run
   output%balances  = .TRUE.  ! energy and water balances
   output%averaging = "${AVGflag}"
   check%ranges     = .FALSE.  ! variable ranges, input and output
   check%energy_bal = .TRUE.  ! energy balance
   check%mass_bal   = .TRUE.  ! water/mass balance
   verbose = .FALSE. ! write details of every grid cell init and params to log?
   leaps = .FALSE. ! calculate timing with leap years?
   logn = 88      ! log file number - declared in input module
   fixedCO2 = ${CO2_var}.0   ! if not found in met file, in ppmv
   spincasainput = .${spincasain}.    ! input required to spin casacnp offline
   spincasa      = .${spincasa}.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.
   l_casacnp     =  .${l_casacnp}. ! using casaCNP with CABLE
   l_laiFeedbk   = .${LAIflag}.  ! using prognostic LAI
   l_vcmaxFeedbk = .${Vcmax_flag}.  ! using prognostic Vcmax
   icycle = ${icycle_flag}   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
   casafile%cnpipool=''
   casafile%cnpbiome='${in_path}/core/biogeochem/${casabiome}'
   casafile%cnpepool='${casaout_path}/poolcnpOut.csv'    ! end of run pool size
   casafile%cnpmetout=''         ! output daily met forcing for spinning casacnp
   casafile%cnpmetin=''          ! list of daily met files for spinning casacnp
   casafile%phen='${in_path}/core/biogeochem/modis_phenology_csiro.txt'
   casafile%cnpflux='${casaout_path}/cnpfluxOut.csv'
   casafile%dump_cnpspin = "${casaout_path}/casa_dump_cnpspin_${year}.nc" !Step 2 output
   casafile%cnpspin = "${casaout_path}/fcnpspin.lst"             !Step 3 input
   output%CASA = .${casaout_flag}. ! output to netcdf file?
   ncciy = 0 ! 0 for not using gswp; 4-digit year input for year of gswp met
   redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
   wiltParam = 0.5
   satuParam = 0.8
   cable_user%FWSOIL_SWITCH = 'standard'        ! choices are: 
                                                ! 1. standard
                                                ! 2. non-linear extrapolation
                                                ! 3. Lai and Katul 2000
                                                ! 4. Haverd2013
   cable_user%DIAG_SOIL_RESP = 'ON ' 
   cable_user%LEAF_RESPIRATION = 'ON ' 
   cable_user%RUN_DIAG_LEVEL= 'BASIC'           ! choices are:
                                                ! 1. BASIC
                                                ! 1. NONE
   cable_user%CONSISTENCY_CHECK= .FALSE.      ! TRUE outputs combined fluxes at each timestep for comparisson to A control run 
   cable_user%CASA_DUMP_READ = .FALSE.      ! TRUE reads CASA forcing from netcdf format
   cable_user%CASA_DUMP_WRITE = .FALSE.      ! TRUE outputs CASA forcing in netcdf format
   cable_user%SSNOW_POTEV= 'HDM'      ! Humidity Deficit Method
   cable_user%L_newProfile = .TRUE.    !Y-P's new root water uptake?
   cable_user%GS_SWITCH = 'M'

&end

EOF

