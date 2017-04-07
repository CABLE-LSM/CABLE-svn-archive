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
   -i      CABLE restart file to be used as input (netcdf)
   -r      CABLE restart file output at end of run
   -v      vegetation data file
   -s      soil data file
   -y      year to run the model (gswp forcing)
   -e      true -> spin to equilibrium
   -g      Start CASA from zero
   -a      CASA restart file to be used as input (netcdf)
   -f      CASA restart file output at end of run
   -c	   icycle option
   -p      CASA output path
   -d      true -> use LAI feedback
   -x      true -> use vxmax feedback
   -t      true -> spincasain option
   -u      true -> spincasa option
   -n      CO2 concentration
   -b      true -> write CASA output to netcdf
   -k      CASA biome parameter file
   -z      Should all variables be outputted?
EOF
}

#set some defaults here?  or do with if;fi after getopts

while getopts "hm:j:o:l:i:r:v:s:y:e:g:a:f:c:p:d:x:t:u:n:b:k:z:" OPTION

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
             casa_zero=$OPTARG
             ;;
         a)
             restart_in_casa=$OPTARG
             ;;
         f)
             restart_out_casa=$OPTARG
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
         z)
             all_vars=$OPTARG
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

if [ -z $casa_zero ]; then   #default to no groundwater
  casa_zero="false"
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
  vegfile='def_veg_params_zr_clitt_correct.txt'
  echo "Using default veg file $vegfile"
fi

if [ -z $soilfile ]; then
  soilfile='def_soil_params.txt'
  echo "Using default soil file $soilfile"
fi

if [ -z $restart_in_casa ]; then
  restart_in_casa='casa_restart_in.nc'
  echo "Using default CASA restart input file $restart_in_casa"
fi

if [ -z $restart_out_casa ]; then   #default to no spinup
  restart_out_casa="casa_restart_out.nc"
  echo "Using default CASA restart output file $restart_in_casa"
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
  CO2_var="390"
fi

if [ -z $casaout_flag ]; then
  casaout_flag="FALSE"
fi

if [ -z $casabiome ];then
  casabiome="pftlookup_csiro_v16_17tiles_Ticket2.csv"
fi
if [ -z $all_vars ]; then
  all_vars="FALSE"
fi

#Set CASA-CNP flag (should CASA be used or not?)
if [[ $icycle_flag -eq 0 ]]; then
    l_casacnp="FALSE"
else
    l_casacnp="TRUE"
fi


in_path="./Inputs"


touch $(pwd)/cable.nml

cat > $(pwd)/cable.nml << EOF
&cable
   filename%met = "${in_path}/$metfile"
   filename%out = "${cableout_path}/${outfile}"
   filename%log = "${cableout_path}/${logfile}"
   filename%restart_in  = "$restart_infile" 
   filename%restart_out = "${cableout_path}/${restart_outfile}"
   filename%type    = "${in_path}/gridinfo_CSIRO_1x1.nc"
   filename%veg    = "${in_path}/$vegfile"
   filename%soil    = "${in_path}/${soilfile}"
   filename%gw_elev = "${in_path}/GSWP3_elevation_slope_mean_stddev_site.nc"
   vegparmnew = .TRUE.  ! using new format when true
   soilparmnew = .TRUE.  ! using new format when true
   spinup = .${eqflag}.  ! do we spin up the model?
   delsoilM = 0.01   ! allowed variation in soil moisture for spin up
   delsoilT = 0.1    ! allowed variation in soil temperature for spin up
   output%restart = .TRUE.  ! should a restart file be created?
   output%met = .${all_vars}.  ! input met data
   output%flux = .${all_vars}.  ! convective, runoff, NEE
   output%soil = .${all_vars}.  ! soil states
   output%snow = .${all_vars}.  ! snow states
   output%radiation = .${all_vars}.  ! net rad, albedo
   output%carbon    = .${all_vars}.  ! NEE, GPP, NPP, stores
   output%veg       = .${all_vars}.  ! vegetation states
   output%params    = .${all_vars}.  ! input parameters used to produce run
   output%balances  = .${all_vars}.  ! energy and water balances
   output%averaging = "all"
   check%ranges     = .FALSE.  ! variable ranges, input and output
   check%energy_bal = .TRUE.  ! energy balance
   check%mass_bal   = .TRUE.  ! water/mass balance
   verbose = .FALSE. ! write details of every grid cell init and params to log?
   leaps = .FALSE. ! calculate timing with leap years?
   logn = 88      ! log file number - declared in input module
   fixedCO2 = ${CO2_var}.0   ! if not found in met file, in ppmv
   l_casacnp     =  .${l_casacnp}. ! using casaCNP with CABLE
   l_laiFeedbk   = .${LAIflag}.  ! using prognostic LAI
   l_vcmaxFeedbk = .${Vcmax_flag}.  ! using prognostic Vcmax
   icycle = ${icycle_flag}   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
   cable_user%CASA_fromZero = .${casa_zero}.
   casafile%cnpipool='${restart_in_casa}'
   casafile%cnpbiome='${in_path}/CASA_inputs/${casabiome}'
   casafile%cnpepool='${casaout_path}/${restart_out_casa}'    ! end of run pool size
   casafile%cnpmetout='casamet.nc'         ! output daily met forcing for spinning casacnp
   casafile%cnpmetin=''          ! list of daily met files for spinning casacnp
   casafile%phen='${in_path}/CASA_inputs/modis_phenology_csiro.txt'
   casafile%cnpflux='${casaout_path}/cnpfluxOut.csv'
   spincasainput = .${spincasain}.    ! input required to spin casacnp offline
   spincasa      = .${spincasa}.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.
   cable_user%CASA_DUMP_READ = .${spincasa}.      ! TRUE reads CASA forcing from netcdf format
   cable_user%CASA_DUMP_WRITE = .${spincasain}.   ! TRUE outputs CASA forcing in netcdf format
   cable_user%CASA_NREP = 0 !number of times to repeat CASA forcing
   !cable_user%CASA_SPIN_STARTYEAR = 1819        ! default = 1950
   !cable_user%CASA_SPIN_ENDYEAR   = 1850        ! default = 1960
   !casafile%dump_cnpspin = "${casaout_path}/casa_dump_cnpspin_${year}.nc" !Step 2 output
   !casafile%cnpspin = "${casaout_path}/fcnpspin.lst"             !Step 3 input
   output%casa = .${casaout_flag}. ! output CASA variables?
   ncciy = 0 ! 0 for not using gswp; 4-digit year input for year of gswp met
   gswpfile%rainf = "./gswp/Rainf/GSWP3.BC.Rainf.3hrMap.${year}.nc"
   gswpfile%snowf = "./gswp/Snowf/GSWP3.BC.Snowf.3hrMap.${year}.nc"
   gswpfile%LWdown= "./gswp/LWdown/GSWP3.BC.LWdown.3hrMap.${year}.nc"
   gswpfile%SWdown= "./gswp/SWdown/GSWP3.BC.SWdown.3hrMap.${year}.nc" 
   gswpfile%PSurf = "./gswp/PSurf/GSWP3.BC.PSurf.3hrMap.${year}.nc"
   gswpfile%Qair  = "./gswp/Qair/GSWP3.BC.Qair.3hrMap.${year}.nc"
   gswpfile%Tair  = "./gswp/Tair/GSWP3.BC.Tair.3hrMap.${year}.nc"
   gswpfile%wind  = "./gswp/Wind/GSWP3.BC.Wind.3hrMap.${year}.nc"
   gswpfile%mask  = "./surface_data/GSWP3_landmask.nc"
   redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
   cable_user%YearStart = 0
   cable_user%YearEnd = 0
   cable_user%MetType='site'
   cable_user%FWSOIL_SWITCH = 'standard'        ! choices are:
                                                ! 1. standard
                                                ! 2. non-linear extrapolation
                                                ! 3. Lai and Katul 2000  (check spelling if fails...)
                                                ! 4. Haverd2013
   cable_user%DIAG_SOIL_RESP = 'ON ' 
   cable_user%LEAF_RESPIRATION = 'ON ' 
   cable_user%RUN_DIAG_LEVEL= 'BASIC'           ! choices are:
                                                ! 1. BASIC
                                                ! 1. NONE
   cable_user%CONSISTENCY_CHECK= .FALSE.      ! TRUE outputs combined fluxes at each timestep for comparisson to A control run 
   cable_user%SSNOW_POTEV= 'HDM'      ! Humidity Deficit Method
   cable_user%GW_MODEL = .TRUE.       !True means use the groundwater module, false means use default soil snow scheme
   cable_user%alt_forcing = .FALSE.
   cable_user%GSWP3 = .FALSE.
   cable_user%GS_SWITCH = 'medlyn'
   cable_user%or_evap = .TRUE.         !Use Or soilE formulation?
   cable_user%litter = .false.         !Use litter scheme to suppress Esoil?
   cable_user%L_NEW_ROUGHNESS_SOIL = .false.
   cable_user%CLIMATE_fromZero=.${casa_zero}.
   cable_user%CALL_CLIMATE = .TRUE.
   cable_user%soil_struc='default'
   wiltParam = 0.5
   satuParam = 0.8
   gw_params%MaxSatFraction     = -1.0
   gw_params%MaxHorzDrainRate   = 4.0e-5! 1e-4!Hi->1e-2 ; Lo->5e-5 ; med->1e-4; med2->1e-3; Hi2->1e-3
   gw_params%EfoldHorzDrainRate = 1.0 !Hi->5.0 ; Lo->1.0; med->2.0 -> med2->3.0; Hi2 ->5.0
   cable_user%sync_netcdf_file = .FALSE.  !Should code write output upon crashing? Slows down code if TRUE
&end

EOF

