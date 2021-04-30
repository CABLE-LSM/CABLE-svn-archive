#!/bin/basha
# Argument = -t test -r server -p password -v

usage()
{
cat << EOF
usage: $0 options

This script setups up cable.nml and can run the model

OPTIONS:
   -h      Show this message
   -m      met file path
   -o      Output file
   -l      log file
   -i      restart file to be used as input (netcdf)
   -r      restart file output at end of run
   -t      typefile = gridinfo file
   -v      vegetation data file
   -s      soil data file
   -y      year to run the model (gswp forcing)
   -e      true -> spin to equilibrium
   -a      averaging flag - daily,monthly,yearly
   -f      alternate (merra) forcing flag
   -c	   CO2 concentration
EOF
}

#set some defaults here?  or do with if;fi after getopts

metdir="/g/data/w35/mm3972/model/cable/src/CABLE-AUX/offline/"

while getopts “hm:r:o:l:i:r:t:v:s:y:e:g:a:f:c:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         m)
             metfile_path=$OPTARG
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
         t)
             typefile=$OPTARG
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
             CO2_con=$OPTARG
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

if [ -z $GWflag ]; then   #default to no groundwater
  GWflag="false"
fi

if [ -z $year ]; then
  year=0
  echo "No year supplied. Setting ncicty=0"
fi

if [ -z $metfile_path ]; then
  metfile_path=""
  echo "Using default met file: $metfile Hopefully this is correct as ncity=0"
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
  if [ $year -gt 1926 ]; then              #gswp run after first year
    let pyear=$year-1
    restart_infile="restart_${pyear}.nc"
  fi
  echo "Using default restart_in file $restart_infile"
fi

if [ -z $restart_outfile ]; then
  if [ $year -eq 0 ]; then
    restart_outfile="restart_out.nc"   #default for non-gswp run
  else
    restart_outfile="restart_${year}.nc"
  fi
  echo "Using default restart_out file $restart_outfile"
fi

if [ -z $typefile ]; then
  typefile='gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc'
  echo "Using default veg file $typefile"
fi

if [ -z $vegfile ]; then
  vegfile='def_veg_params_zr_clitt_albedo_fix.txt'
  echo "Using default veg file $vegfile"
fi

if [ -z $soilfile ]; then
  soilfile='def_soil_params.txt'
  echo "Using default soil file $soilfile"
fi

if [ -z $AVGflag ]; then
  AVGflag='monthly'
  echo "Setting the output averaging flag to monthly"
fi

if [ -z $ALTFORCEflag ]; then   #default to no spinup
  ALTFORCEflag="false"
fi

if [ -z $CO2_con ];then
  CO2_con=370.0
fi

#if [ $(pwd)/cable.nml ];then
#  rm $(pwd)/cable.nml
#fi

touch $(pwd)/cable_on.nml

cat > $(pwd)/cable_on.nml << EOF
&cable
   filename%met = ""
   filename%out = "${outfile}"
   filename%log = "${logfile}"
   filename%restart_in  = "${restart_infile}"
   filename%restart_out = "${restart_outfile}"
   filename%type = "/g/data/w35/amu561/Steven_CABLE_runs/gridinfo_AWRA_matched_veg_inputs.nc"
   filename%veg  = "/g/data/w35/amu561/Steven_CABLE_runs/CABLE_source/CABLE-AUX/core/biogeophys/def_veg_params_zr_clitt_albedo_fix.txt"
   filename%gw_elev ="/g/data/w35/amu561/Steven_CABLE_runs/gridinfo_AWRA_matched_veg_inputs.nc"
   filename%soil = "/g/data/w35/amu561/Steven_CABLE_runs/CABLE_source/CABLE-AUX/core/biogeophys/def_soil_params.txt"
   vegparmnew = .TRUE.  ! using new format when true
   soilparmnew = .TRUE.  ! using new format when true
   spinup = .FALSE.  ! do we spin up the model?
   delsoilM = 0.001   ! allowed variation in soil moisture for spin up
   delsoilT = 0.01    ! allowed variation in soil temperature for spin up
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
   verbose = .TRUE. ! write details of every grid cell init and params to log?
   leaps = .TRUE. ! calculate timing with leap years?
   logn = 88      ! log file number - declared in input module
   fixedCO2 = 350.0   ! MMY fixed co2 to compare with LIS-CABLE if not found in met file, in ppmv
   spincasa      = .FALSE.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.
   l_casacnp     = .FALSE.  ! using casaCNP with CABLE
   l_laiFeedbk   = .FALSE.  ! using prognostic LAI
   l_vcmaxFeedbk = .FALSE.  ! using prognostic Vcmax
   icycle = 0   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
   casafile%cnpipool=' *** SET PATH IN cable.nml *** '
   casafile%cnpbiome=' *** SET PATH IN cable.nml *** '
   casafile%cnpepool='poolcnpOut.csv'    ! end of run pool size
   casafile%cnpmetout='casamet.nc'                ! output daily met forcing for spinning casacnp
   casafile%cnpmetin=''          ! list of daily met files for spinning casacnp
   casafile%phen=' *** SET PATH IN cable.nml *** '
   casafile%cnpflux='cnpfluxOut.csv'
   ncciy = $year ! 0 for not using gswp; 4-digit year input for year of gswp met
   gswpfile%rainf = "${metfile_path}/Rainf/GSWP3.BC.Rainf.3hrMap.${year}.nc"
   gswpfile%snowf = "${metfile_path}/Snowf/GSWP3.BC.Snowf.3hrMap.${year}.nc"
   gswpfile%LWdown= "${metfile_path}/LWdown/GSWP3.BC.LWdown.3hrMap.${year}.nc"
   gswpfile%SWdown= "${metfile_path}/SWdown/GSWP3.BC.SWdown.3hrMap.${year}.nc"
   gswpfile%PSurf =  "/g/data/w35/amu561/Steven_CABLE_runs/gridinfo_AWRA_matched_veg_inputs.nc" 
   gswpfile%Qair  = "${metfile_path}/Qair/GSWP3.BC.Qair.3hrMap.${year}.nc"
   gswpfile%Tair  = "${metfile_path}/Tair/GSWP3.BC.Tair.3hrMap.${year}.nc"
   gswpfile%wind  = "${metfile_path}/Wind/GSWP3.BC.Wind.3hrMap.${year}.nc"
   gswpfile%mask  = "/g/data/w35/amu561/Steven_CABLE_runs/gridinfo_AWRA_matched_veg_inputs.nc"
   redistrb = .FALSE.  ! Turn on/off the hydraulic redistribution
   wiltParam = 0.5
   satuParam = 0.8
   cable_user%FWSOIL_SWITCH = 'Haverd2013'        ! choices are:
                                                 ! 1. standard
                                                 ! 2. non-linear extrapolation
                                                 ! 3. Lai and Katul 2000
                                                 ! 4. Haverd2013
   cable_user%DIAG_SOIL_RESP = 'ON '
   cable_user%LEAF_RESPIRATION = 'ON '
   cable_user%RUN_DIAG_LEVEL= 'BASIC'        ! choices are:
                                                 ! 1. BASIC
                                                 ! 1. NONE
   cable_user%CONSISTENCY_CHECK = .FALSE.      ! TRUE outputs combined fluxes at each timestep for comparisson to A control run
   cable_user%CASA_DUMP_READ = .FALSE.        ! TRUE reads CASA forcing from netcdf format
   cable_user%CABLE_RUNTIME_COUPLED = .FALSE.
   cable_user%CASA_DUMP_WRITE = .FALSE.       ! TRUE outputs CASA forcing in netcdf format
   cable_user%SSNOW_POTEV = 'HDM'              ! Humidity Deficit Method
   cable_user%GW_MODEL = .FALSE.          ! True means use the groundwater module, false means use default soil snow scheme
   cable_user%GSWP3 = .TRUE.
   cable_user%MetType = 'gswp3'
   cable_user%GS_SWITCH = 'medlyn'
   cable_user%or_evap = .FALSE.
   cable_user%soil_struc = 'default'
   cable_user%litter = .FALSE.
   cable_user%call_climate = .false.
   cable_user%POPLUC = .false.
   cable_user%soil_thermal_fix = .true.
   gw_params%MaxSatFraction = 3000.0
   gw_params%MaxHorzDrainRate = 2e-5
   gw_params%EfoldHorzDrainRate = 1.0
   gw_params%hkrz = 0.1
   gw_params%zdepth = 2.0
   gw_params%ssgw_ice_switch = .true.

&end

EOF
