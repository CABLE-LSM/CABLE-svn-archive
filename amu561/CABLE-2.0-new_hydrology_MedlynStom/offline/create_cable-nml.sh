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
   -o      Output file
   -l      log file 
   -i      restart file to be used as input (netcdf)
   -r      restart file output at end of run
   -v      vegetation data file
   -s      soil data file
   -y      year to run the model (gswp forcing)
   -e      true -> spin to equilibrium
   -a      averaging flag - daily,monthly,yearly
   -f      alternate (merra) forcing flag
EOF
}

#set some defaults here?  or do with if;fi after getopts

while getopts “hm:r:o:l:i:r:v:s:y:e:g:a:f:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         m)
             metfile=$OPTARG
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

if [ -z $metfile] && [ $year -eq 0 ]; then
  metfile="default_met.txt"
  echo "Using default met file: $metfile Hopefully this is correct as ncity=0"
fi

if [ -z $outfile ]; then
  outfile="out_cable.nc"          #if gswp append year to log
fi

if [ -z $logfile ]; then
  logfile='log_cable.txt'
  if [ $year -ne 0 ]; then
      logfile="cable_log_${year}.txt"          #if gswp append year to log
  fi
  echo "Using default log file $logfile"
fi

if [ -z $restart_infile ]; then
#  restart_infile="restart_2009.nc"
  if [ $year -gt 1979 ]; then              #gswp run after first year
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

if [ -z $vegfile ]; then
  vegfile='def_veg_params.txt'
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

if [ -f cable.nml ]; then
  mv cable.nml old.cable.nml
fi

touch cable.nml

cat > cable.nml << EOF
&cable
   filename%met = "$metfile"
   filename%out = "$outfile"
   filename%log = "$logfile"
   filename%restart_in  = "$restart_infile" 
   filename%restart_out = "$restart_outfile"
   filename%type    = "./../../Inputs/CABLE_GSWP3_HGSD_DRT_Surface_Data.nc"
   filename%veg    = "./../../Inputs/$vegfile"
   filename%soil    = "./../../Inputs/$soilfile"
   filename%gw_elev = "$metfile"
   vegparmnew = .TRUE.  ! using new format when true
   soilparmnew = .FALSE.  ! using new format when true     !Anna changed to FALSE
   spinup = .${eqflag}.  ! do we spin up the model?
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
   verbose = .FALSE. ! write details of every grid cell init and params to log?
   leaps = .false. ! calculate timing with leap years?
   logn = 88      ! log file number - declared in input module
   fixedCO2 = 350.0   ! if not found in met file, in ppmv
   spincasainput = .FALSE.    ! input required to spin casacnp offline
   spincasa      = .FALSE.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.
   l_casacnp     = .FALSE.  ! using casaCNP with CABLE 
   l_laiFeedbk   = .FALSE.  ! using prognostic LAI
   l_vcmaxFeedbk = .FALSE.  ! using prognostic Vcmax
   icycle = 0   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P
   casafile%cnpipool=' *** SET PATH IN cable.nml *** '
   casafile%cnpbiome=' *** SET PATH IN cable.nml *** '
   casafile%cnpepool=''    ! end of run pool size
   casafile%cnpmetout=''                ! output daily met forcing for spinning casacnp
   casafile%cnpmetin=''          ! list of daily met files for spinning casacnp
   casafile%phen=' *** SET PATH IN cable.nml *** '
   casafile%cnpflux='cnpfluxOut.csv'
   ncciy = $year ! 0 for not using gswp; 4-digit year input for year of gswp met
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
   cable_user%GS_SWITCH = 'leuning'             ! choices are:
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
   cable_user%GW_MODEL = .${GWflag}.       !True means use the groundwater module, false means use default soil snow scheme
   cable_user%alt_forcing = .${ALTFORCEflag}.
   cable_user%GSWP3       = .FALSE.
   gw_params%EfoldMaxSatFrac = 6.0
   gw_params%MaxHorzDrainRate = 0.01
   gw_params%MaxSatFraction = 0.0  !Changed Anna
   gw_params%EfoldHorzDrainRate = 1.0
&end

EOF



