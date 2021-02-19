#!/bin/bash
#
# Explor / Pearcey
# https://slurm.schedmd.com/sbatch.html
# Name - 8 letters and digits
#SBATCH -J harvar10
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.out
# Explor partitions (sinfo): std (2x16, parallel), sky (2x16, parallel, AVX512), hf (2x4, serial),
#                            P100 (2x16, GPU), GTX (2x16, GPU), ivy (2x8, parallel), k20 (2x8, GPU)
#SBATCH -p sky
# Nodes / tasks
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --ntasks-per-node=4
# Check memory on *nix with /usr/bin/time -v ./prog
# time (day-hh:mm:ss) / memory (optional, units K,M,G,T)
#SBATCH -t 00:09:59
#SBATCH --mem=4G
# notify: Valid type values are NONE,BEGIN,END,FAIL,REQUEUE,ALL,STAGE_OUT,TIME_LIMIT,TIME_LIMIT_90/80/50,ARRAY_TASKS
#SBATCH --mail-type=FAIL,STAGE_OUT,TIME_LIMIT
#SBATCH --mail-user=matthias.cuntz@inrae.fr
#
# Gadi
# https://opus.nci.org.au/display/Help/How+to+submit+a+job
#PBS -N plume_test
#PBS -P x45
# express / normal / copyq (2x24, cascadelake)
# expressbw / normalbw (2x14, broadwell) / normalsl (2x16, skylake)- ex-Raijin nodes
#PBS -q normal
#PBS -l walltime=17:59:59
#PBS -l mem=192GB
#PBS -l ncpus=32
# #PBS -l jobfs=1GB
#PBS -l storage=gdata/x45
#PBS -l software=netCDF:MPI:Intel:GNU
#PBS -r y
#PBS -l wd
#PBS -j oe
#PBS -S /bin/bash
#PBS -M juergen.knauer@csiro.au

# cuntz@explor, cuntz@mcinra, moc801@gadi cuntz@gadi
# kna016@pearcey knauer@pearcey, jk8585@gadi knauer@gadi
# not yet vxh599@gadi nor hav014@pearcey
system=knauer@gadi

# MPI run or single processor run
# nproc should fit with job tasks
dompi=1   # 0: normal run: ./cable
# 1: MPI run: mpiexec -n ${nproc} ./cable_mpi
nproc=32   # Number of cores for MPI runs
# must be same as above: SBATCH -n nproc or PBS -l ncpus=nproc
ignore_mpi_err=0 # 0/1: 1: continue even if mpi run failed


# --------------------------------------------------------------------
#
# Full Cable run with biomass spinup, POP, land-use change, etc.
#
# This script uses CRU-JRA forcing.
#
# Global meteo and land-use change data can be used with a mask giving land points.
# In step 0, the land mask can be extracted for a single point, an area
# or a number of random points can be chosen.
# Alternatively, single site met, LUH2, forcing and mask can be extracted from the global data sets.
#
# The run sequence is as follows:
#   1. Create a climate restart file using Cable's default vegetation distribution.
#   2. First phase of spinup with static land use, fixed atmospheric CO2 and
#      N deposition from 1700, and 30 years of repeated meteorology.
#      Zero initial biomass stocks, POP and climate.
#   3. Bring biomass stocks into equilibrium, restricting labile P and mineral N pools.
#      Repeat a and b several times.
#      a) Start from restart files.
#      b) Use dump files for the biophysics. Quasi-equilibrium of soil and litter pools
#         using an analytic steady-state solution.
#   4. Same as 3 but without any restriction on labile P and mineral N pools.
#      Repeat a and b several times.
#   5. Second phase of spinup with dynamic land use, atmospheric CO2 and N deposition.
#      a) Dynamic land use from 1580 to 1699, using still fixed atmospheric CO2 and
#         N deposition from 1700, and 30 years of repeated meteorology.
#      b) Run from 1700 to 1899 with dynmic land use, varying atmospheric CO2 and N deposition,
#         but still with 30 years of repeated meteorology.
#   6. Final historical run, everything dynamic from 1900 to 2017.
#   7. Final future run, everything dynamic (not all met types).
#
# Written, Matthias Cuntz, Aug 2019, following the run scripts and namelists provided by Vanessa Haverd
#
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Sequence switches
#
imeteo=1        # 0: Use global meteo, land use and mask
                # 1: Use local mask, but global meteo and land use (doextractsite=1)
                # 2: Use local meteo, land use and mask (doextractsite=2)
# Step 0
doextractsite=1 # 0: Do not extract local meteo, land use nor mask
                # 1: Do extract only mask at specific site/region (imeteo=1)
                # 2: Do extract meteo, land use and mask at specific site/region (imeteo=2)
                #    Does not work with randompoints /= 0 but with latlon
#experiment=Test_parallel  # assigned later!
purge_restart=1 # 1/0: Do/Do not delete all restart files (for a completely new run, e.g. if settings changed)
Test_parallel=1 # 1/0: Test run or not?
randompoints=0   # <0: use -1*randompoints random grid points from ${LandMaskFilePath}/${experiment}_points.csv if existing
                 # 0:  use latlon
                 # >0: generate and use randompoints random grid points from GlobalLandMaskFile
# lat,lon  or  latmin,latmax,lonmin,lonmax   # must have . in numbers otherwise indexes taken
# latlon=42.536875,-72.172602
#latlon=-34.5,-33.5,149.5,156.5
#latlon=-3.0,4.5,15.5,25.5
latlon=-44.0,-10.0,110.0,155.0  # Australia
# Step 1
doclimate=0     # 1/0: Do/Do not create climate restart file
# Step 2
dofromzero=0    # 1/0: Do/Do not first spinup phase from zero biomass stocks
# Step 3
doequi1=0       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with restricted P and N pools
nequi1=3        #      number of times to repeat steps in doequi1
# Step 4
doequi2=0       # 1/0: Do/Do not bring biomass stocks into quasi-equilibrium with unrestricted P and N pools
nequi2=3        #      number of times to repeat steps in doequi2
# Step 5a
doiniluc=0      # 1/0: Do/Do not spinup with dynamic land use (5a)
# Step 5b
doinidyn=0      # 1/0: Do/Do not full dynamic spinup from 1700 to 1899 (5b)
# Step 6
dofinal=0       # 1/0: Do/Do not final run from 1900 to 2017
# Step 7
dofuture=0      # 1/0: Do/Do not future runs (plume only)

# --------------------------------------------------------------------
# Other switches
#
# MetType
mettype='plume'     # 'cru', 'plume', 'bios'
metmodel='hadgem2'  # 'hadgem2', 'ipsl' (only used if mettype='plume')
RCP='8.5'           # 'hist', '2.6', '4.5', '6.0', '8.5' (no future runs if RCP='hist')

# Cable
doc13o2=0             # 1/0: Do/Do not calculate 13C
c13o2_simple_disc=0   # 1/0: simple or full 13C leaf discrimination
explicit_gm=0         # 1/0: explicit (finite) or implicit mesophyll conductance
coordinate_photosyn=0 # 1/0: Do/Do not coordinate photosynthesis
acclimate_photosyn=0  # 1/0: Do/Do not acclimate photosynthesis
use_LUTgm=1           # 1/0: Do/Do not use lookup table for parameter conversion accounting for gm? (only used if explicit_gm=1)
coord=F               # T/F: version of photosyn. optimisation (optimised(F) or forced (T))
Rubisco_params='Bernacchi_2002'   # 'Bernacchi_2002' or 'Walker_2013'

# --------------------------------------------------------------------
# Setup
#

set -e

trap cleanup 1 2 3 6

pid=$$
isdir="${PWD}"
prog=$0
pprog=$(basename ${prog})
pdir=$(dirname ${prog})
tmp=${TMPDIR:-"/tmp"}

system=$(echo ${system} | tr A-Z a-z)
sys=${system#*@}
user=${system%@*}

#
# Special things on specific computer system such as loading modules
#
export mpiexecdir=
if [[ "${sys}" == "explor" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    # # INTELMPI - load mpi module first, otherwise intel module will not prepend LD_LIBRARY_PATH
    # module load intelmpi/2018.5.274
    # module load intel/2018.5
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    # export mpiexecdir=/soft/env/soft/all/intel/2018.3/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/
    # # INTEL / OpenMPI - load mpi module first, otherwise intel module will not pre-pend LD_LIBRARY_PATH
    # module load openmpi/3.0.0/intel18
    # module load intel/2018.5
    # export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-ifort2018.0/lib
    # export mpiexecdir=/opt/soft/hf/openmpi-3.0.0-intel18/bin/
    # GNU / OpenMPI
    module load gcc/6.3.0
    module load openmpi/3.0.1/gcc/6.3.0
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/local/lib:${HOME}/local/netcdf-fortran-4.4.4-gfortran63/lib
    export mpiexecdir=/opt/soft/hf/openmpi/3.0.1/gcc/6.3.0/bin/
    if [[ ${doextractsite} -ge 1 ]] ; then module load python/intel/2019/3 ; fi
elif [[ "${sys}" == "mcinra" ]] ; then
    # # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-gfortran"
    # export mpiexecdir=/usr/local/openmpi-3.1.4-gfortran/bin/
    # # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-ifort"
    export mpiexecdir=/usr/local/openmpi-3.1.5-ifort/bin/
elif [[ "${sys}" == "pearcey" ]] ; then
    # prog is slurm_script
    pdir=${isdir}
    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module unload intel-mpi/5.0.1.035
    module add netcdf/4.3.3.1 openmpi/1.8.8
elif [[ "${sys}" == "raijin" ]] ; then
    module del intel-cc intel-fc
    module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
    module add netcdf/4.3.3.1
elif [[ "${sys}" == "gadi" ]] ; then
    pdir=${isdir}
    . /etc/bashrc
    module purge
    module load intel-compiler/2019.5.281
    module load intel-mpi/2019.5.281
    module load netcdf/4.6.3
    module load hdf5/1.10.5
    if [[ ${doextractsite} -ge 1 ]] ; then
        module load python3/3.7.4
        export PYTHONPATH=${PYTHONPATH}:/g/data/x45/python/lib/python3.7/site-packages
    fi
    if [[ ${randompoints} -eq 0 ]] ; then
        module load nco/4.9.2  # needed for cropping outputs
    fi
    export mpiexecdir=/apps/intel-mpi/2019.5.281/intel64/bin/
fi
if [[ ! -z ${mpiexecdir} ]] ; then export mpiexecdir="${mpiexecdir}/" ; fi

#
# Directories of things
#

## Experiment naming for gm_acclim_coord experiment
experiment=test_${mettype}
if [[ ("${system}" == "jk8585@gadi" || "${system}" == "knauer@gadi") && "${mettype}" == "plume" ]] ; then
    if [[ ${explicit_gm} -eq 1 ]] ; then 
        if [[ ${coordinate_photosyn} -eq 1 ]] ; then 
            if [[ ${acclimate_photosyn} -eq 1 ]] ; then
                experiment=exp_coord${coord}_PSacclim
            else
	        experiment=exp_coord${coord}
            fi
        else
	    if [[ ${acclimate_photosyn} -eq 1 ]] ; then
	        experiment=exp_nocoord_PSacclim
	    else
	        experiment=exp_nocoord
	    fi	
	fi
    else
	if [[ ${coordinate_photosyn} -eq 1 ]] ; then
	    if [[ ${acclimate_photosyn} -eq 1 ]] ; then
		experiment=imp_coord${coord}_PSacclim
	    else
		experiment=imp_coord${coord}
	    fi
	else
	    if [[ ${acclimate_photosyn} -eq 1 ]] ; then
		experiment=imp_nocoord_PSacclim
	    else
		experiment=imp_nocoord
	    fi
        fi
    fi

    if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
        experiment=${experiment}_B2002
    elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
        experiment=${experiment}_W2013
    else
	echo "unknwn option for 'Rubisco_params', exiting..."
	exit 1
    fi
    
    if [[ ${Test_parallel} -eq 1 ]] ; then
	experiment=TEST_${experiment}
    fi
fi


#
# Relative directories must be relative to the directory of this script,
#   not relative to the directory from which this script is launched (if different)
#   nor relative to the run path.
#
if [[ "${system}" == "cuntz@explor" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/home/oqx29/zzy20/prog/cable/single_sites/${experiment}"
    cablehome="/home/oqx29/zzy20/prog/cable"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi"
    else
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablehome}/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    # GlobalLandMaskFile="/home/oqx29/zzy20/data/crujra/daily_1deg/glob_ipsl_1x1.nc"
    # Global CRU
    GlobalMetPath="/home/oqx29/zzy20/data/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/home/oqx29/zzy20/data/LUH2_v3_1deg/"
elif [[ "${system}" == "cuntz@mcinra" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/Users/cuntz/prog/vanessa/cable/single_sites/${experiment}"
    cablehome="/Users/cuntz/prog/vanessa/cable"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-ifort"
        # exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi-gfortran"
    else
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="${cablehome}/CABLE-AUX"
    # Global Mask, CRU, LUC
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    GlobalMetPath=
    GlobalTransitionFilePath=
elif [[ "${system}" == "moc801@gadi" || "${system}" == "cuntz@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    # sitepath="/home/801/moc801/prog/cable/single_sites/${experiment}"
    sitepath="/scratch/x45/moc801/cable/speed1000"
    cablehome="/home/801/moc801/prog/cable"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable-mpi"
    else
        exe="${cablehome}/branches/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    # Global CRU
    GlobalMetPath="/g/data/x45/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2018/1deg/EXTRACT"
elif [[ "${system}" == "kna016@pearcey" || "${system}" == "knauer@pearcey" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_run/parallel_runs/${experiment}"
    cablehome="/OSM/CBR/OA_GLOBALCABLE/work/Juergen/CABLE_code"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/NESP2pt9_BLAZE/offline/cable-mpi"
    else
        exe="${cablehome}/NESP2pt9_BLAZE/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/OSM/CBR/OA_GLOBALCABLE/work/Vanessa/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    # Global CRU
    GlobalMetPath="/OSM/CBR/OA_GLOBALCABLE/work/CRU-JRA55/crujra/daily_1deg"
    # Global LUC
    GlobalTransitionFilePath="/OSM/CBR/OA_GLOBALCABLE/work/LUH2/v3/1deg"
elif [[ "${system}" == "jk8585@gadi" || "${system}" == "knauer@gadi" ]] ; then
    # Run directory: runpath="${sitepath}/run_xxx"
    sitepath="/g/data/x45/jk8585/global_runs/gm_acclim_coord/${experiment}"
    workpath="/home/599/jk8585/CABLE_run/gm_acclim_coord/global_runs"
    cablehome="/home/599/jk8585/CABLE_code"
    # Cable executable
    if [[ ${dompi} -eq 1 ]] ; then
        exe="${cablehome}/gm_acclim_coord/offline/cable-mpi"
	#exe="${cablehome}/NESP2pt9_BLAZE/offline/cable-mpi"
    else
	#exe="${cablehome}/NESP2pt9_BLAZE/offline/cable"
        exe="${cablehome}/gm_acclim_coord/offline/cable"
    fi
    # CABLE-AUX directory (uses offline/gridinfo_CSIRO_1x1.nc and offline/modis_phenology_csiro.txt)
    aux="/g/data/x45/CABLE-AUX"
    # Global Mask
    GlobalLandMaskFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"
    SurfaceFile="${aux}/offline/gridinfo_CSIRO_1x1.nc"   # note that SurfaceFile does not need subsetting
    # Global Met
    if [[ "${mettype}" == 'plume' ]] ; then
	GlobalMetPath="/g/data/x45/ipbes/${metmodel}/1deg"
        CO2Path="/g/data/x45/ipbes/co2"
	NdepPath="/g/data/x45/ipbes/ndep"
	GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
	GlobalLandMaskFile="/g/data/x45/ipbes/masks/gridinfo_CSIRO_1x1.nc"
    elif [[ "${mettype}" == 'bios' ]] ; then
	GlobalLandMaskFile="/g/data/x45/BIOS3_forcing/acttest9/acttest9" # no file extension
        GlobalMetPath="/g/data/x45/BIOS3_forcing/acttest9/met/"  # last slash is needed
	ParamPath="/g/data/x45/BIOS3_forcing/acttest9/params/"
    elif [[ "${mettype}" == 'cru' ]] ; then
	GlobalLandMaskFile="/g/data/x45/ipbes/masks/glob_ipsl_1x1.nc"
        GlobalMetPath="/g/data/x45/CRUJRA2020/daily_1deg"
    fi
    # Global LUC
    GlobalTransitionFilePath="/g/data/x45/LUH2/GCB_2019/1deg/EXTRACT"
else
    echo "System not known."
    exit 1
fi
# Run directory
runpath="${sitepath}/run"

# Cable parameters
namelistpath="$(dirname ${workpath})/namelists"
filename_veg="$(dirname ${workpath})/params/def_veg_params.txt"
filename_soil="$(dirname ${workpath})/params/def_soil_params.txt"
casafile_cnpbiome="$(dirname ${workpath})/params/pftlookup.csv"
# Other scripts
ScriptsPath="${workpath}/scripts"
# Mask
LandMaskFile="${sitepath}/mask/${experiment}_landmask.nc"
# Met forcing
MetPath="${sitepath}/met/cru_jra_1deg"
ClimateFile="${sitepath}/run/${mettype}_climate_rst.nc"
if [[ ${doclimate} -eq 0 && ! -f ${ClimateFile} ]] ; then
   ClimateFile="/g/data/x45/ipbes/cable_climate/ipsl_climate_rst_glob_1deg.nc"
fi 

# LUC
TransitionFilePath="${sitepath}/LUH2/v3/1deg"
# 13C
filename_d13c_atm="$(dirname ${workpath})/params/graven_et_al_gmd_2017-table_s1-delta_13c-1700-2025.txt"

# --------------------------------------------------------------------
# Start Script
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Helper functions
#
# usage of script
function usage()
{
    printf "${pprog} [-h]\n"
    printf "Runs Cable on a single grid cell with spinup, POP, land-use change, etc.\n"
    printf "Behaviour of the script is controlled by switches at the top of the script (ca. line 101ff).\n"
    printf "\n"
    printf "Options\n"
    printf "    -h    Prints this help screen.\n"
}
#
# cleanup at end or at trap
function cleanup()
{
    \rm -f ${tmp}/*.${pid}*
    exit 1
}
#
# absolute path
function abspath()
{
    idir=${PWD}
    cd ${1}
    odir=${PWD}
    cd ${idir}
    echo "${odir}"
}
#
# filename with absolute path
function absfile()
{
    f=$(basename ${1})
    d=$(dirname ${1})
    d=$(abspath ${d})
    echo "${d}/${f}"
}
#
# rm file if present
function irm()
{
    for i in "$@" ; do
        if [[ -f ${i} ]] ; then rm ${i} ; fi
    done
}
#
# sed command from comma-separated list: var1=str1,var2=str2
#   no , but = allowed in str
function csed()
{
    com=""
    for i in $(echo ${1} | tr ',' ' ') ; do
        v=$(echo ${i} | cut -d '=' -f 1)
        s=$(echo ${i} | cut -d '=' -f 2-)
        # com="${com} -e '/^[[:blank:]]*${v}[[:blank:]]*=/s/${v}[[:blank:]]*=.*/${v} = ${s}/'"
        # com="${com} -e 's/${v}[[:blank:]]*=.*/${v} = ${s}/'"
        com="${com} -e s|${v}[[:blank:]]*=.*|${v}=${s}|"
    done
    printf "%s" "${com}"
}
#
# get list of variables in file with nco (cdo omits the dimension variables)
function ncvarlist()
{
    out=$(ncks --trd -m ${1} | grep -E ': type' | cut -f 1 -d ' ' | sed 's/://' | sort)
    echo ${out}
}
#
# check if first argument string is in rest of arguments
# returns first argument if present otherwise empty string (check with -z)
function isin()
{
    tofind=${1}
    shift
    out=""
    for i in $@ ; do
        if [[ ${i} == ${tofind} ]] ; then out=${i} ; fi
    done
    echo ${out}
}
#
# returns argument to extract lat and lon with ncks
function nckslatlon()
{
    vars=$(ncvarlist ${1})
    if [[ -z $(isin latitude ${vars}) ]] ; then ilat='lat' ; else ilat='latitude' ; fi
    if [[ -z $(isin longitude ${vars}) ]] ; then ilon='lon' ; else ilon='longitude' ; fi
    if [[ -z $(echo ${2} | cut -f 3 -d ',') || -z $(echo ${2} | cut -f 4 -d ',') ]] ; then
        iilat=$(echo ${2} | cut -f 1 -d ',')
        iilon=$(echo ${2} | cut -f 2 -d ',')
        echo "-d ${ilat},${iilat} -d ${ilon},${iilon}"
    else
        iilat1=$(echo ${2} | cut -f 1 -d ',')
        iilat2=$(echo ${2} | cut -f 2 -d ',')
        iilon1=$(echo ${2} | cut -f 3 -d ',')
        iilon2=$(echo ${2} | cut -f 4 -d ',')
        echo "-d ${ilat},${iilat1},${iilat2} -d ${ilon},${iilon1},${iilon2}"
    fi
}
#
# copy files adding first argument to filenames
function copyid()
{
    rid=${1}
    shift 1
    for i in $@ ; do
        cp ${i} ${i%.*}_${rid}.${i##*.}
    done
}
#
# rename files adding first argument to filenames
function renameid()
{
    rid=${1}
    shift 1
    for i in $@ ; do
        mv ${i} ${i%.*}_${rid}.${i##*.}
    done
}

# --------------------------------------------------------------------------------------------------
# Preparation
#
# Get options
while getopts "h" option ; do
    case ${option} in
        h) usage; exit;;
        *) printf "Error ${pprog}: unimplemented option.\n\n" 1>&2;  usage 1>&2; exit 1;;
    esac
done
shift $((${OPTIND} - 1))

#
# get directories
pdir=$(abspath ${pdir})
cd ${pdir}
adir=$(abspath ${aux})
exe=$(absfile ${exe})
mkdir -p ${runpath}
rdir=$(abspath ${runpath})
ndir=$(abspath ${namelistpath})
sdir=$(abspath ${ScriptsPath})

#
# prepare run directory
cd ${rdir}
mkdir -p logs
mkdir -p outputs
mkdir -p restart
ln -sf ${adir}
ln -sf ${exe}
iexe=$(basename ${exe})
cd ${pdir}

#
# set stacksize to unlimited if permitted, otherwise to 15 bit if possible
set +e
ulimit -s unlimited 2> /dev/null || ulimit -s 32768
set -e

# --------------------------------------------------------------------
# Info
#
t1=$(date +%s)
printf "Started at %s\n" "$(date)"

printf "\nSetup\n"
printf "    Serial / Parallel\n"
printf "        dompi=${dompi}\n"
printf "            nproc=${nproc}\n"
printf "            ignore_mpi_err=${ignore_mpi_err}\n"
printf "\n"
printf "    Sequence\n"
printf "        imeteo=${imeteo}\n"
printf "        doextractsite=${doextractsite}\n"
printf "            experiment=${experiment}\n"
printf "            randompoints=${randompoints}\n"
printf "            latlon=${latlon}\n"
printf "        doclimate=${doclimate}\n"
printf "        dofromzero=${dofromzero}\n"
printf "        doequi1=${doequi1}\n"
printf "            nequi1=${nequi1}\n"
printf "        doequi2=${doequi2}\n"
printf "            nequi2=${nequi2}\n"
printf "        doiniluc=${doiniluc}\n"
printf "        doinidyn=${doinidyn}\n"
printf "        dofinal=${dofinal}\n"
printf "        dofuture=${dofuture}\n"
printf "\n"
printf "    Options\n"
printf "        mettype=${mettype}\n"
printf "        RCP=${RCP}\n"
printf "        doc13o2=${doc13o2}\n"
printf "        c13o2_simple_disc=${c13o2_simple_disc}\n"
printf "        explicit_gm=${explicit_gm}\n"
printf "        use_LUTgm=${use_LUTgm}\n"
printf "        coordinate_photosyn=${coordinate_photosyn}\n"
printf "        acclimate_photosyn=${acclimate_photosyn}\n"
printf "        Rubisco_params=${Rubisco_params}\n"
printf "\n"
printf "    Directories\n"
printf "        sitepath=${sitepath}\n"
printf "        cablehome=${cablehome}\n"
printf "        exe=${exe}\n"
printf "        aux=${aux}\n"
printf "        GlobalLandMaskFile=${GlobalLandMaskFile}\n"
printf "        GlobalMetPath=${GlobalMetPath}\n"
printf "        GlobalTransitionFilePath=${GlobalTransitionFilePath}\n"
printf "        runpath=${runpath}\n"
printf "        namelistpath=${namelistpath}\n"
printf "        filename_veg=${filename_veg}\n"
printf "        filename_soil=${filename_soil}\n"
printf "        casafile_cnpbiome=${casafile_cnpbiome}\n"
printf "        LandMaskFile=${LandMaskFile}\n"
printf "        SurfaceFile=${SurfaceFile}\n"
printf "        MetPath=${MetPath}\n"
printf "        ClimateFile=${ClimateFile}\n"
printf "        TransitionFilePath=${TransitionFilePath}\n"
printf "        filename_d13c_atm=${filename_d13c_atm}\n"
printf "\n"

# --------------------------------------------------------------------
# Sequence
#

# 0. Extract meteo, land use and mask for one specific site from global files
if [[ ${doextractsite} -ge 2 ]] ; then
    # set +e
    # xcdo=$(which cdo 2> /dev/null)
    # set -e
    # if [[ -z ${xcdo} ]] ; then module load cdo ; fi
    set +e
    xnco=$(which ncks 2> /dev/null)
    set -e
    if [[ -z ${xnco} ]] ; then module load nco ; fi
fi
if [[ ${doextractsite} -eq 1 ]] ; then
    echo "0. Set local mask"
    cd ${pdir}
    # mask
    LandMaskFilePath=$(dirname ${LandMaskFile})
    mkdir -p ${LandMaskFilePath}
    LandMaskFile=$(absfile ${LandMaskFile})
    # echo $(basename ${LandMaskFile})
    # generate random points if ${randompoints} > 0
    if [[ ${randompoints} -ne 0 ]] ; then
        dogeneraterandom=1
        rpoints=${randompoints}
        if [[ ${randompoints} -lt 0 ]] ; then
            if [[ -f ${LandMaskFilePath}/${experiment}_points.csv ]] ; then dogeneraterandom = 0 ; fi
            rpoints=$(( ${randompoints} * -1 ))
        fi
        # generate random points
        com=$(csed "basepath=\"${sitepath}\"")
        com=${com}$(csed "gridinfo_file=\"${GlobalLandMaskFile}\"")
        com=${com}$(csed "outname=\"${LandMaskFilePath}/${experiment}_points.csv\"")
        sed ${com} ${sdir}/generate_latlonlist.py > ${LandMaskFilePath}/generate_latlonlist.py
        python3 ${LandMaskFilePath}/generate_latlonlist.py ${rpoints}

        # set mask to generated random points
        com=$(csed "path=\"${LandMaskFilePath}\"")
        com=${com}$(csed "maskfname=\"${LandMaskFile}\"")
        com=${com}$(csed "latlonfile=\"${LandMaskFilePath}/${experiment}_points.csv\"")
        com=${com}$(csed "gridinfo_file=\"${GlobalLandMaskFile}\"")
	if [[ "${mettype}" == "bios" ]] ; then
           com=${com}$(csed "res=0.25")
	fi    
        sed ${com} ${sdir}/create_landmask.py > ${LandMaskFilePath}/create_landmask.py
        sed -i -e "s|from lnutils.*|sys.path.insert(1,'${sdir}'); from lnutils import latlon2ixjy|" ${LandMaskFilePath}/create_landmask.py
        python3 ${LandMaskFilePath}/create_landmask.py
    else
        # # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${GlobalLandMaskFile} ${LandMaskFile}
        # ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}
        com=$(csed "path=\"${LandMaskFilePath}\"")
        com=${com}$(csed "maskfname=\"${LandMaskFile}\"")
        com=${com}$(csed "gridinfo_file=\"${GlobalLandMaskFile}\"")
	if [[ "${mettype}" == "bios" ]] ; then
           com=${com}$(csed "res=0.25")
	fi 
        sed ${com} ${sdir}/create_landmask.py > ${LandMaskFilePath}/create_landmask.py
        sed -i -e "s|from lnutils.*|sys.path.insert(1,'${sdir}'); from lnutils import latlon2ixjy|" ${LandMaskFilePath}/create_landmask.py
        python3 ${LandMaskFilePath}/create_landmask.py ${latlon}
    fi
fi
if [[ ${doextractsite} -eq 2 ]] ; then
    echo "0. Extract local meteo and mask"
    cd ${pdir}
    # meteorology
    met_list="pre pres dlwrf dswrf spfh tmax tmin ugrd vgrd ndep"
    mkdir -p ${MetPath}
    MetPath=$(abspath ${MetPath})
    for mm in ${met_list} ; do
        echo ${GlobalMetPath}/${mm}
        mkdir -p ${MetPath}/${mm}
        for nc in ${GlobalMetPath}/${mm}/*.nc ; do
            ff=$(basename ${nc})
            echo "    ${ff}"
            # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${nc} ${MetPath}/${mm}/${ff}
            ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${MetPath}/${mm}/${ff}
        done
    done
    echo ${GlobalMetPath}/co2/*.csv
    mkdir -p ${MetPath}/co2
    cp ${GlobalMetPath}/co2/*.csv ${MetPath}/co2/

    # land use
    mkdir -p ${TransitionFilePath}
    TransitionFilePath=$(abspath ${TransitionFilePath})
    echo ${GlobalTransitionFilePath}
    for nc in ${GlobalTransitionFilePath}/*.nc ; do
        ff=$(basename ${nc})
        echo "    ${ff}"
        # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${nc} ${TransitionFilePath}/${ff}
        ncks -O $(nckslatlon ${nc} ${latlon}) ${nc} ${TransitionFilePath}/${ff}
    done

    # mask
    LandMaskFilePath=$(dirname ${LandMaskFile})
    mkdir -p ${LandMaskFilePath}
    LandMaskFile=$(absfile ${LandMaskFile})
    echo $(basename ${LandMaskFile})
    # cdo -s -f nc4 -z zip sellonlatbox,-72.5,-72.0,42.5,43.0 ${GlobalLandMaskFile} ${LandMaskFile}
    ncks -O $(nckslatlon ${GlobalLandMaskFile} ${latlon}) ${GlobalLandMaskFile} ${LandMaskFile}
fi

# Choose meteo, land use and mask directories and files
if [[ ${imeteo} -eq 0 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${GlobalLandMaskFile})
elif [[ ${imeteo} -eq 1 ]] ; then
    MetPath=$(abspath ${GlobalMetPath})
    TransitionFilePath=$(abspath ${GlobalTransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
elif [[ ${imeteo} -eq 2 ]] ; then
    MetPath=$(abspath ${MetPath})
    TransitionFilePath=$(abspath ${TransitionFilePath})
    LandMaskFile=$(absfile ${LandMaskFile})
else
    printf "Error ${pprog}: imeteo option unknown: ${imeteo}.\n\n"
    exit 1
fi
# absolute paths of other parameter files
ClimateFile=$(absfile ${ClimateFile})
filename_veg=$(absfile ${filename_veg})
filename_soil=$(absfile ${filename_soil})
casafile_cnpbiome=$(absfile ${casafile_cnpbiome})
filename_d13c_atm=$(absfile ${filename_d13c_atm})


# global changes to cable.nml file (for all steps below, but can be overwritten in later steps)
com=$(csed "filename%veg=\"${filename_veg}\",filename%soil=\"${filename_soil}\"")
com=${com}$(csed "filename%type=\"${SurfaceFile}\"")
com=${com}$(csed "filename%out=\"outputs/${mettype}_out_cable.nc\"")
com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
com=${com}$(csed "filename%restart_out=\"restart/${mettype}_cable_rst.nc\"")
com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
com=${com}$(csed "casafile%cnpepool=\"restart/${mettype}_casa_rst.nc\"")
com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
com=${com}$(csed "cable_user%POP_restart_out=\"restart/pop_${mettype}_ini.nc\"")
com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
com=${com}$(csed "cable_user%climate_restart_in=\"restart/${mettype}_climate_rst.nc\"")
com=${com}$(csed "cable_user%climate_restart_out=\"restart/${mettype}_climate_rst.nc\"")
com=${com}$(csed "cable_user%RunIden=\"${mettype}\"")
com=${com}$(csed "cable_user%MetType=\"${mettype}\"")
com=${com}$(csed "cable_user%CALL_POP=.TRUE.")
com=${com}$(csed "output%averaging=\"monthly\"")
com=${com}$(csed "output%grid=\"mask\"")
com=${com}$(csed "leaps=.FALSE.")
com=${com}$(csed "cable_user%SOIL_STRUC=\"default\"")
com=${com}$(csed "cable_user%Rubisco_parameters=\"${Rubisco_params}\"")
if [[ ${explicit_gm} -eq 1 ]] ; then
    com=${com}$(csed "cable_user%explicit_gm=.true.")
    if [[ ${use_LUTgm} -eq 1 ]] ; then
	if [[ "${Rubisco_params}" == "Bernacchi_2002" ]] ; then
	    com=${com}$(csed "cable_user%gm_LUT_file=\"$(dirname ${workpath})/params/gm_LUT_351x3601x7_1pt8245_Bernacchi2002.nc\"")
	elif [[ "${Rubisco_params}" == "Walker_2013" ]] ; then
	    com=${com}$(csed "cable_user%gm_LUT_file=\"$(dirname ${workpath})/params/gm_LUT_351x3601x7_1pt8245_Walker2013.nc\"")
	fi    
    else
       com=${com}$(csed "cable_user%gm_LUT_file=\"\"")
    fi
else
    com=${com}$(csed "cable_user%explicit_gm=.false.")
fi
if [[ ${coordinate_photosyn} -eq 1 ]] ; then
    com=${com}$(csed "cable_user%coordinate_photosyn=.true.")
else
    com=${com}$(csed "cable_user%coordinate_photosyn=.false.")
fi
if [[ ${acclimate_photosyn} -eq 1 ]] ; then
    com=${com}$(csed "cable_user%acclimate_photosyn=.true.")
else
    com=${com}$(csed "cable_user%acclimate_photosyn=.false.")
fi

sed ${com} ${ndir}/cable.nml > ${ndir}/cable_global_${experiment}.nml



## global plume.nml file
if [[ "${mettype}" == 'plume' ]] ; then
    com=$(csed "BasePath=\"${MetPath}\"")
    if [[ "${metmodel}" == 'ipsl' ]] ; then
        com=${com}$(csed "Forcing=\"ipsl-cm5a-lr\"")
    elif [[ "${metmodel}" == 'hadgem2' ]] ; then
	com=${com}$(csed "Forcing=\"hadgem2-es\"")
    fi
    
    sed ${com} ${ndir}/plume.nml > ${ndir}/plume_${metmodel}.nml
fi


## global bios.nml file
if [[ "${mettype}" == 'bios' ]] ; then
    com=$(csed "met_path=\"${MetPath}/\",param_path=\"${ParamPath}\"")
    com=${com}$(csed "landmaskflt_file=\"${GlobalLandMaskFile}.flt\"")
    com=${com}$(csed "landmaskhdr_file=\"${GlobalLandMaskFile}.hdr\"")
    com=${com}$(csed "rain_file=\"1900010120201231_rain_b2003.bin\"")
    com=${com}$(csed "swdown_file=\"1900010120201231_rad_b2003.bin\"")
    com=${com}$(csed "tairmax_file=\"1900010120201231_tmax_noclim_b2003.bin\"")
    com=${com}$(csed "tairmin_file=\"1900010120201231_tmin_noclim_b2003.bin\"")
    com=${com}$(csed "wind_file=\"1900010120201231_windspeed_ms_b2003.bin\"")
    com=${com}$(csed "vp0900_file=\"1900010120201231_vph09_b2003.bin\"")
    com=${com}$(csed "vp1500_file=\"1900010120201231_vph15_b2003.bin\"")
    com=${com}$(csed "co2_file=\"1750_2020_globalCO2_time_series.bin\"")
    sed ${com} ${ndir}/bios.nml > ${ndir}/bios_global.nml
fi




# delete all restart files if required
if [[ ${purge_restart} -eq 1 ]] ; then
    rm -f ${rdir}/restart/*
fi


# 1. Create climate restart file
if [[ ${doclimate} -eq 1 ]] ; then
    echo "1. Create climate restart file"
    rid="climate_restart"
    # Met forcing
    if [[ "${mettype}" == 'cru' ]] ; then
        irm ${rdir}/cru.nml
        com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
        com=${com}$(csed "Run=\"S0_TRENDY\"")
        sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    elif [[ "${mettype}" == 'plume' ]] ; then
        irm ${rdir}/plume.nml
	com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	com=${com}$(csed "Run=\"spinup\"")
	com=${com}$(csed "RCP=\"hist\"")
	com=${com}$(csed "CO2=\"static1850\"")
	com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	com=${com}$(csed "NDEP=\"static1850\"")
        com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
    elif [[ "${mettype}" == 'bios' ]] ; then
        irm ${rdir}/bios.nml
	com=$(csed "Run=\"spinup\"")
	sed ${com} ${ndir}/bios_global.nml > ${rdir}/bios.nml
    fi	
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1700")
    com=${com}$(csed "YearEnd=2017")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.true.")
    com=${com}$(csed "cable_user%YearStart=1860")
    com=${com}$(csed "cable_user%YearEnd=1889")
    #com=${com}$(csed "cable_user%YearEnd=1925")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
    com=${com}$(csed "cable_user%limit_labile=.true.")
    com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"\"")
    com=${com}$(csed "cable_user%POP_fromZero=.true.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"\"")
    com=${com}$(csed "cable_user%POPLUC=.false.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")	
    # do not calculate 13C because there is no 13C in the climate restart file
    com=${com}$(csed "cable_user%c13o2=.false.")
    sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} cable.nml ${mettype}.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} ${mettype}_climate_rst.nc
    cp ${mettype}_climate_rst.nc ${ClimateFile}
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    cd ..
    cd ${pdir}
else
  cp ${ClimateFile} ${mettype}_climate_rst.nc
fi 


# 2. First spinup phase from zero biomass
if [[ ${dofromzero} -eq 1 ]] ; then
    echo "2. First spinup from zero biomass"
    rid="zero_biomass"
    # Met forcing
    if [[ "${mettype}" == 'cru' ]] ; then
        irm ${rdir}/cru.nml
        com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
        com=${com}$(csed "Run=\"S0_TRENDY\"")
        sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    elif [[ "${mettype}" == 'plume' ]] ; then
        irm ${rdir}/plume.nml
	com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	com=${com}$(csed "Run=\"spinup\"")
	com=${com}$(csed "RCP=\"hist\"")
	com=${com}$(csed "CO2=\"static1850\"")
	com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	com=${com}$(csed "NDEP=\"static1850\"")
        com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
    fi 
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=1700")
    com=${com}$(csed "YearEnd=2017")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.true.")
    com=${com}$(csed "cable_user%YearStart=1860")
    com=${com}$(csed "cable_user%YearEnd=1889")
    #com=${com}$(csed "cable_user%YearEnd=1861")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
    com=${com}$(csed "cable_user%limit_labile=.true.")
    com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"\"")
    com=${com}$(csed "cable_user%POP_fromZero=.true.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} cable.nml ${mettype}.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 3. Biomass into quasi-equilibrium with unrestricted N and P pools
if [[ ${doequi1} -eq 1 ]] ; then
    echo "3. Bring biomass into quasi-equilibrium with unrestricted N and P pools"
    for ((iequi1=1; iequi1<=${nequi1}; iequi1++)) ; do
        rid="spinup_limit_labile"
        # rid="spinup_limit_labile${iequi}"
        if [[ 1 -eq 1 ]] ; then
            # 3a. 30 year run starting from restart files
            echo "   3a. 30 year spinup from accumulated biomass; iequi1=${iequi1}/${nequi1}"
            # Met forcing
            if [[ "${mettype}" == 'cru' ]] ; then
	        irm ${rdir}/cru.nml
                com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
                com=${com}$(csed "Run=\"S0_TRENDY\"")
                sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
	    elif [[ "${mettype}" == 'plume' ]] ; then
                irm ${rdir}/plume.nml
	        com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	        com=${com}$(csed "Run=\"spinup\"")
	        com=${com}$(csed "RCP=\"hist\"")
	        com=${com}$(csed "CO2=\"static1850\"")
	        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	        com=${com}$(csed "NDEP=\"static1850\"")
                com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	        sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
            fi 	
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            #com=${com}$(csed "cable_user%YearEnd=1841")
            com=${com}$(csed "icycle=2")
            com=${com}$(csed "spincasa=.false.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
            com=${com}$(csed "cable_user%limit_labile=.true.")
            com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
                com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
            irm logs/log_cable.txt logs/log_out_cable.txt
            if [[ ${dompi} -eq 1 ]] ; then
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
                ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
            else
                ./${iexe} > logs/log_out_cable.txt
            fi
            # save output
            renameid ${rid} cable.nml ${mettype}.nml LUC.nml
            mv *_${rid}.nml restart/
            cd logs
            renameid ${rid} log_cable.txt log_out_cable.txt
            cd ../restart
            copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
            if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
            cd ../outputs
            renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
            if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
            cd ..
            cd ${pdir}
        fi
        #
        if [[ 1 -eq 1 ]] ; then
            # 3b. analytic quasi-equilibrium of biomass pools
            echo "   3b. Analytic solution of biomass pools"
            rid="spinup_analytic_limit_labile"
            # rid="spin_casa_limit_labile${iequi}"
            # Met forcing
	    if [[ "${mettype}" == 'cru' ]] ; then
                irm ${rdir}/cru.nml
                com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
                com=${com}$(csed "Run=\"S0_TRENDY\"")
                sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
	    elif [[ "${mettype}" == 'plume' ]] ; then
                irm ${rdir}/plume.nml
	        com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	        com=${com}$(csed "Run=\"spinup\"")
	        com=${com}$(csed "RCP=\"hist\"")
	        com=${com}$(csed "CO2=\"static1850\"")
	        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	        com=${com}$(csed "NDEP=\"static1850\"")
                com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	        sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
            fi
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            #com=${com}$(csed "cable_user%YearEnd=1841")
            com=${com}$(csed "icycle=12")
            com=${com}$(csed "spincasa=.true.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
            #com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1841")
            com=${com}$(csed "cable_user%limit_labile=.true.")
            com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
                com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
            irm logs/log_cable.txt logs/log_out_cable.txt
            if [[ ${dompi} -eq 1 ]] ; then
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
                ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
            else
                ./${iexe} > logs/log_out_cable.txt
            fi
            # save output
            renameid ${rid} cable.nml ${mettype}.nml LUC.nml
            mv *_${rid}.nml restart/
            cd logs
            renameid ${rid} log_cable.txt log_out_cable.txt
            cd ../restart
            copyid ${rid} pop_${mettype}_ini.nc ${mettype}_casa_rst.nc
            if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
            if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
                cd ../outputs
                renameid ${rid} ${mettype}_out_casa.nc
                if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
                cd ..
            fi
            cd ${pdir}
        fi
    done
fi


# 4. Biomass into quasi-equilibrium with restricted N and P pools
if [[ ${doequi2} -eq 1 ]] ; then
    echo "4. Bring biomass into quasi-equilibrium with restricted N and P pools"
    for ((iequi2=1; iequi2<=${nequi2}; iequi2++)) ; do
        rid="spinup"
        # rid="spinup${iequi}"
        if [[ 1 -eq 1 ]] ; then
            # 4a. 30 year run starting from restart files
            echo "   4a. 30 year spinup from accumulated biomass; iequi2=${iequi2}/${nequi2}"
            # Met forcing
	    if [[ "${mettype}" == 'cru' ]] ; then
                irm ${rdir}/cru.nml
                com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
                com=${com}$(csed "Run=\"S0_TRENDY\"")
                sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
	    elif [[ "${mettype}" == 'plume' ]] ; then
                irm ${rdir}/plume.nml
	        com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	        com=${com}$(csed "Run=\"spinup\"")
	        com=${com}$(csed "RCP=\"hist\"")
	        com=${com}$(csed "CO2=\"static1850\"")
	        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	        com=${com}$(csed "NDEP=\"static1850\"")
                com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	        sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
            fi
	    # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            #com=${com}$(csed "cable_user%YearEnd=1841")
            com=${com}$(csed "icycle=2")
            com=${com}$(csed "spincasa=.false.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.true.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1860")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1869")
            com=${com}$(csed "cable_user%limit_labile=.false.")
            com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
                com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
            irm logs/log_cable.txt logs/log_out_cable.txt
            if [[ ${dompi} -eq 1 ]] ; then
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
                ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
            else
                ./${iexe} > logs/log_out_cable.txt
            fi
            # save output
            renameid ${rid} cable.nml ${mettype}.nml LUC.nml
            mv *_${rid}.nml restart/
            cd logs
            renameid ${rid} log_cable.txt log_out_cable.txt
            cd ../restart
            copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
            if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
            cd ../outputs
            renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
            if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
            cd ..
            cd ${pdir}
        fi
        #
        if [[ 1 -eq 1 ]] ; then
            # 4b. analytic quasi-equilibrium of biomass pools
            echo "   4b. Analytic solution of biomass pools"
            rid="spinup_analytic"
            # rid="spin_casa${iequi}"
            # Met forcing
	    if [[ "${mettype}" == 'cru' ]] ; then
                irm ${rdir}/cru.nml
                com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
                com=${com}$(csed "Run=\"S0_TRENDY\"")
                sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
            elif [[ "${mettype}" == 'plume' ]] ; then
                irm ${rdir}/plume.nml
	        com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	        com=${com}$(csed "Run=\"spinup\"")
	        com=${com}$(csed "RCP=\"hist\"")
	        com=${com}$(csed "CO2=\"static1850\"")
	        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	        com=${com}$(csed "NDEP=\"static1850\"")
                com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	        sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
            fi
            # LUC
            irm ${rdir}/LUC.nml
            com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
            com=${com}$(csed "YearStart=1700")
            com=${com}$(csed "YearEnd=2017")
            sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
            # CABLE
            irm ${rdir}/cable.nml
            com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
            com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
            com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
            com=${com}$(csed "cable_user%YearStart=1840")
            com=${com}$(csed "cable_user%YearEnd=1859")
            #com=${com}$(csed "cable_user%YearEnd=1841")
            com=${com}$(csed "icycle=12")
            com=${com}$(csed "spincasa=.true.")
            com=${com}$(csed "cable_user%CASA_fromZero=.false.")
            com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
            com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
            com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
            com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")
            com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
            #com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1841")
            com=${com}$(csed "cable_user%limit_labile=.false.")
            com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
            com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
            com=${com}$(csed "cable_user%POP_fromZero=.false.")
            com=${com}$(csed "cable_user%POP_out=\"ini\"")
            com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
            com=${com}$(csed "cable_user%POPLUC=.true.")
            com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
            com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
            com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
            if [[ ${doc13o2} -eq 1 ]] ; then
                com=${com}$(csed "cable_user%c13o2=.true.")
                com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
                com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
                if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
                else
                    com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
                fi
            fi
            sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
            # run model
            cd ${rdir}
            irm logs/log_cable.txt logs/log_out_cable.txt
            if [[ ${dompi} -eq 1 ]] ; then
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
                ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
                if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
            else
                ./${iexe} > logs/log_out_cable.txt
            fi
            # save output
            renameid ${rid} cable.nml ${mettype}.nml LUC.nml
            mv *_${rid}.nml restart/
            cd logs
            renameid ${rid} log_cable.txt log_out_cable.txt
            cd ../restart
            copyid ${rid} pop_${mettype}_ini.nc ${mettype}_casa_rst.nc
            if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
            if [[ ${dompi} -eq 0 ]] ; then # no output only restart if MPI
                cd ../outputs
                renameid ${rid} ${mettype}_out_casa.nc
                if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
                cd ..
            fi
            cd ${pdir}
        fi
    done
fi

# 5a. First dynamic land use
if [[ ${doiniluc} -eq 1 ]] ; then
    echo "5a. First dynamic land use"
    # Met forcing
    if [[ "${mettype}" == 'cru' ]] ; then
	rid="1580_1699"
        YearStart=1580
        YearEnd=1699
        irm ${rdir}/cru.nml
        com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
        com=${com}$(csed "Run=\"S0_TRENDY\"")
        sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    elif [[ "${mettype}" == 'plume' ]] ; then
	rid="1580_1849"
	YearStart=1580
        YearEnd=1849
	irm ${rdir}/plume.nml
	com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	com=${com}$(csed "Run=\"spinup\"")
	com=${com}$(csed "RCP=\"hist\"")
	com=${com}$(csed "CO2=\"static1850\"")
        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	com=${com}$(csed "NDEP=\"static1850\"")
        com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
    fi
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=${YearStart}")
    com=${com}$(csed "YearEnd=${YearEnd}")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
    com=${com}$(csed "cable_user%YearStart=${YearStart}")
    com=${com}$(csed "cable_user%YearEnd=${YearEnd}")
    com=${com}$(csed "icycle=12")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.true.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"annually\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1840")  # not used here I assume
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
    com=${com}$(csed "cable_user%limit_labile=.false.")
    com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
    com=${com}$(csed "cable_user%POP_fromZero=.false.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"init\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} cable.nml ${mettype}.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_${mettype}_ini.nc ${mettype}_casa_rst.nc ${mettype}_LUC_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
    # cd ../outputs
    # renameid ${rid} ${mettype}_out_LUC.nc
    # cd ..
    cd ${pdir}
fi


# 5b. Second full dynamic spinup
if [[ ${doinidyn} -eq 1 ]] ; then
    echo "5b. Full dynamic spinup"
    # Met forcing
    if [[ "${mettype}" == 'cru' ]] ; then
	rid="1700_1899"
        YearStart=1700
        YearEnd=1899
        irm ${rdir}/cru.nml
        com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
        com=${com}$(csed "Run=\"S1_TRENDY\"")
        sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
    elif [[ "${mettype}" == 'plume' ]] ; then
	rid="1850_1900"
        YearStart=1850
        YearEnd=1900
        irm ${rdir}/plume.nml
	com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	com=${com}$(csed "Run=\"1850_1900\"")
	com=${com}$(csed "RCP=\"hist\"")
	com=${com}$(csed "CO2=\"varying\"")
        com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	com=${com}$(csed "NDEP=\"varying\"")
        com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
    fi
    # LUC
    irm ${rdir}/LUC.nml
    com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
    com=${com}$(csed "YearStart=${YearStart}")
    com=${com}$(csed "YearEnd=${YearEnd}")
    sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    # CABLE
    irm ${rdir}/cable.nml
    com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
    com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
    com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
    com=${com}$(csed "cable_user%YearStart=${YearStart}")
    com=${com}$(csed "cable_user%YearEnd=${YearEnd}")
    #com=${com}$(csed "cable_user%YearEnd=1851")
    com=${com}$(csed "icycle=2")
    com=${com}$(csed "spincasa=.false.")
    com=${com}$(csed "cable_user%CASA_fromZero=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
    com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
    com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
    #MCTEST com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"daily\"")
    #MCTEST com=${com}$(csed "output%averaging=\"all\"")
    com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1000")  # not used here I assume
    com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1010")
    com=${com}$(csed "cable_user%limit_labile=.false.")
    com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
    com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
    com=${com}$(csed "cable_user%POP_fromZero=.false.")
    com=${com}$(csed "cable_user%POP_out=\"ini\"")
    com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
    com=${com}$(csed "cable_user%POPLUC=.true.")
    com=${com}$(csed "cable_user%POPLUC_RunType=\"restart\"")
    com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
    com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
    if [[ ${doc13o2} -eq 1 ]] ; then
        com=${com}$(csed "cable_user%c13o2=.true.")
        com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
        com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_in_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
        com=${com}$(csed "cable_user%c13o2_restart_out_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
        if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
        else
            com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
        fi
    fi
    sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
    # run model
    cd ${rdir}
    irm logs/log_cable.txt logs/log_out_cable.txt
    if [[ ${dompi} -eq 1 ]] ; then
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
        ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
        if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
    else
        ./${iexe} > logs/log_out_cable.txt
    fi
    # save output
    renameid ${rid} cable.nml ${mettype}.nml LUC.nml
    mv *_${rid}.nml restart/
    cd logs
    renameid ${rid} log_cable.txt log_out_cable.txt
    cd ../restart
    copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc ${mettype}_LUC_rst.nc
    if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
    cd ../outputs
    renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
    if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
    cd ..
    cd ${pdir}
fi


# 6. Final centennial run
if [[ ${dofinal} -eq 1 ]] ; then
    echo "6. Final centennial run"

    YearStartPrev=$YearStart
    YearEndPrev=$YearEnd
    if [[ "${mettype}" == 'cru' ]] ; then
        YearStart=1900
        YearEnd=2017
	chunksize=300
    elif [[ "${mettype}" == 'plume' ]] ; then
        YearStart=1901
        YearEnd=2005
	chunksize=300
    fi

    nextchunk=1
    eyear=$YearStart
    while (( ${nextchunk} > 0 )) ; do
        syear=$(echo "${eyear} + 1" | bc)
        eyear=$(echo "${eyear} + ${chunksize}" | bc)
	if [[ ${syear} -eq $(echo "${YearStart} + 1" | bc) ]] ; then
	    syear=$(echo "${syear} - 1" | bc)
	    eyear=$(echo "${eyear} - 1" | bc)
	fi
        if [[ ${eyear} -gt ${YearEnd} ]] ; then
	    eyear=$YearEnd
	    nextchunk=0
	fi
        rid="${syear}_${eyear}"

        if [[ -f ${rdir}/restart/${mettype}_cable_rst_${rid}.nc ]] ; then
            echo "skipping" $rid
	    continue
	fi
	
        echo "  simulating years" $syear "-" $eyear 

        # Reset restart files
        if [[ ${syear} -eq ${YearStart} ]] ; then
	    YearStartRestart=$YearStartPrev
	    YearEndRestart=$YearEndPrev
	else
	    YearStartRestart=$(echo "${syear} - ${chunksize}" | bc)
 	    YearEndRestart=$(echo "${syear} - 1" | bc)
	fi
	resid="${YearStartRestart}_${YearEndRestart}"
	cp ${rdir}/restart/${mettype}_cable_rst_${resid}.nc ${rdir}/restart/${mettype}_cable_rst.nc
        cp ${rdir}/restart/${mettype}_casa_rst_${resid}.nc ${rdir}/restart/${mettype}_casa_rst.nc
        cp ${rdir}/restart/${mettype}_climate_rst_${resid}.nc ${rdir}/restart/${mettype}_climate_rst.nc
        cp ${rdir}/restart/${mettype}_LUC_rst_${resid}.nc ${rdir}/restart/${mettype}_LUC_rst.nc
	cp ${rdir}/restart/pop_${mettype}_ini_${resid}.nc ${rdir}/restart/pop_${mettype}_ini.nc
	    
	# Met forcing
	if [[ "${mettype}" == 'cru' ]] ; then
            irm ${rdir}/cru.nml
            com=$(csed "BasePath=\"${MetPath}\",MetPath=\"${MetPath}\",LandMaskFile=\"${LandMaskFile}\"")
            com=${com}$(csed "Run=\"S2_TRENDY\"")
            sed ${com} ${ndir}/cru.nml > ${rdir}/cru.nml
	elif [[ "${mettype}" == 'plume' ]] ; then
            irm ${rdir}/plume.nml
	    com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	    com=${com}$(csed "Run=\"1901_2005\"")
	    com=${com}$(csed "RCP=\"hist\"")
	    com=${com}$(csed "CO2=\"varying\"")
            com=${com}$(csed "CO2file=\"${CO2Path}/co2_1850_2005_hist.dat\"")
	    com=${com}$(csed "NDEP=\"varying\"")
            com=${com}$(csed "NDEPfile=\"${NdepPath}/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc\"")
	    sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml
	fi
	
	# LUC
	irm ${rdir}/LUC.nml
	com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
	com=${com}$(csed "YearStart=${YearStart}")
	com=${com}$(csed "YearEnd=${YearEnd}")
	sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml
    
        # CABLE
        irm ${rdir}/cable.nml
	com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
	com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
	com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
	com=${com}$(csed "cable_user%YearStart=${syear}")
	com=${com}$(csed "cable_user%YearEnd=${eyear}")
	#com=${com}$(csed "cable_user%YearEnd=1902")
	com=${com}$(csed "icycle=2")
	com=${com}$(csed "spincasa=.false.")
	com=${com}$(csed "cable_user%CASA_fromZero=.false.")
	com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
	com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
	com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
	com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1850")
	com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
	com=${com}$(csed "cable_user%limit_labile=.false.")
	com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
	com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
	com=${com}$(csed "cable_user%POP_fromZero=.false.")
	com=${com}$(csed "cable_user%POP_out=\"ini\"")
	com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
	com=${com}$(csed "cable_user%POPLUC=.true.")
	com=${com}$(csed "cable_user%POPLUC_RunType=\"restart\"")
	com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
	com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
	com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
	if [[ ${doc13o2} -eq 1 ]] ; then
            com=${com}$(csed "cable_user%c13o2=.true.")
            com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
            com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_in_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
            com=${com}$(csed "cable_user%c13o2_restart_out_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
            if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
		com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
            else
		com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
            fi
	fi
	sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml
	# run model
	cd ${rdir}
	irm logs/log_cable.txt logs/log_out_cable.txt
	if [[ ${dompi} -eq 1 ]] ; then
            if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
            if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	else
            ./${iexe} > logs/log_out_cable.txt
	fi
        # save output
        renameid ${rid} cable.nml ${mettype}.nml LUC.nml
        mv *_${rid}.nml restart/
        cd logs
        renameid ${rid} log_cable.txt log_out_cable.txt
        cd ../restart
        copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc ${mettype}_LUC_rst.nc
        if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ${mettype}_c13o2_luc_rst.nc ; fi
        cd ../outputs
        renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc ${mettype}_out_LUC.nc
        if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
        cd ..
        cd ${pdir}
    done	
fi


# 7. Future run (not all Met types)
if [[ ${dofuture} -eq 1 && "${mettype}" == 'plume' && "${RCP}" != 'hist' ]] ; then
    echo "7. Future run using climate projections"
    if [[ ${dofinal} -eq 1 ]] ; then 
      YearStartPrev=$syear
      YearEndPrev=$eyear
    else
      YearStartPrev=1901
      YearEndPrev=2005
    fi

    rcpnd=$(echo ${RCP} | sed "s!\.!!")
    rcpd=$(echo ${RCP} | sed "s!\.!p!")
    YearStart=2006
    YearEnd=2099
    chunksize=300

    nextchunk=1
    eyear=$YearStart
    while (( ${nextchunk} > 0 )) ; do
        syear=$(echo "${eyear} + 1" | bc)
        eyear=$(echo "${eyear} + ${chunksize}" | bc)
	if [[ ${syear} -eq $(echo "${YearStart} + 1" | bc) ]] ; then
	    syear=$(echo "${syear} - 1" | bc)
	    eyear=$(echo "${eyear} - 1" | bc)
	fi
        if [[ ${eyear} -gt ${YearEnd} ]] ; then
	    eyear=$YearEnd
	    nextchunk=0
	fi
        rid="${syear}_${eyear}_${rcpd}"

        if [[ -f restart/${mettype}_cable_rst_${rid}.nc ]] ; then
            echo "skipping" $rid
	    continue
	fi
	
        echo "  simulating years" $syear "-" $eyear

        # Reset restart files
        if [[ ${syear} -eq ${YearStart} ]] ; then
	    YearStartRestart=$YearStartPrev
	    YearEndRestart=$YearEndPrev
	    resid="${YearStartRestart}_${YearEndRestart}"
	else
	    YearStartRestart=$(echo "${syear} - ${chunksize}" | bc)
 	    YearEndRestart=$(echo "${syear} - 1" | bc)
	    resid="${YearStartRestart}_${YearEndRestart}_${rcpd}"
	fi				    
	cp ${rdir}/restart/${mettype}_cable_rst_${resid}.nc ${rdir}/restart/${mettype}_cable_rst.nc
        cp ${rdir}/restart/${mettype}_casa_rst_${resid}.nc ${rdir}/restart/${mettype}_casa_rst.nc
        cp ${rdir}/restart/${mettype}_climate_rst_${resid}.nc ${rdir}/restart/${mettype}_climate_rst.nc
        #cp ${rdir}/restart/${mettype}_LUC_rst_${resid}.nc ${rdir}/restart/${mettype}_LUC_rst.nc   # static LUC for future runs
	cp ${rdir}/restart/pop_${mettype}_ini_${resid}.nc ${rdir}/restart/pop_${mettype}_ini.nc
	
	# Met forcing
	irm ${rdir}/plume.nml
	com=$(csed "LandMaskFile=\"${LandMaskFile}\"")
	com=${com}$(csed "Run=\"2006_2099\"")
	com=${com}$(csed "RCP=\"${RCP}\"")
	com=${com}$(csed "CO2=\"varying\"")
	com=${com}$(csed "CO2file=\"${CO2Path}/co2_2006_2100_rcp${rcpnd}.dat\"")
	com=${com}$(csed "NDEP=\"varying\"")
	com=${com}$(csed "NDEPfile=\"${NdepPath}/RCP${rcpnd}/ndep_total_2000-2109_1.0x1.0_FD.nc\"")
	sed ${com} ${ndir}/plume_${metmodel}.nml > ${rdir}/plume.nml

	# LUC
	irm ${rdir}/LUC.nml
	com=$(csed "TransitionFilePath=\"${TransitionFilePath}\",ClimateFile=\"${ClimateFile}\"")
	com=${com}$(csed "YearStart=${YearStart}")
	com=${com}$(csed "YearEnd=${YearEnd}")
	sed ${com} ${ndir}/LUC.nml > ${rdir}/LUC.nml	

	# CABLE
	irm ${rdir}/cable.nml
	com=$(csed "casafile%cnpbiome=\"${casafile_cnpbiome}\"")
	com=${com}$(csed "filename%restart_in=\"restart/${mettype}_cable_rst.nc\"")
	com=${com}$(csed "cable_user%CLIMATE_fromZero=.false.")
	com=${com}$(csed "cable_user%YearStart=${syear}")
	com=${com}$(csed "cable_user%YearEnd=${eyear}")
	#com=${com}$(csed "cable_user%YearEnd=2008")
	com=${com}$(csed "icycle=2")
	com=${com}$(csed "spincasa=.false.")
	com=${com}$(csed "cable_user%CASA_fromZero=.false.")
	com=${com}$(csed "cable_user%CASA_DUMP_READ=.false.")
	com=${com}$(csed "cable_user%CASA_DUMP_WRITE=.false.")
	com=${com}$(csed "cable_user%CASA_OUT_FREQ=\"monthly\"")
	com=${com}$(csed "cable_user%CASA_SPIN_STARTYEAR=1850")
	com=${com}$(csed "cable_user%CASA_SPIN_ENDYEAR=1859")
	com=${com}$(csed "cable_user%limit_labile=.false.")
	com=${com}$(csed "casafile%out=\"outputs/${mettype}_out_casa.nc\"")
	com=${com}$(csed "casafile%cnpipool=\"restart/${mettype}_casa_rst.nc\"")
	com=${com}$(csed "cable_user%POP_fromZero=.false.")
	com=${com}$(csed "cable_user%POP_out=\"ini\"")
	com=${com}$(csed "cable_user%POP_restart_in=\"restart/pop_${mettype}_ini.nc\"")
	com=${com}$(csed "cable_user%POPLUC=.true.")
	com=${com}$(csed "cable_user%POPLUC_RunType=\"static\"")
	com=${com}$(csed "cable_user%LUC_outfile=\"outputs/${mettype}_out_LUC.nc\"")
	com=${com}$(csed "cable_user%LUC_restart_in=\"restart/${mettype}_LUC_rst.nc\"")
	com=${com}$(csed "cable_user%LUC_restart_out=\"restart/${mettype}_LUC_rst.nc\"")
	if [[ ${doc13o2} -eq 1 ]] ; then
	    com=${com}$(csed "cable_user%c13o2=.true.")
	    com=${com}$(csed "cable_user%c13o2_delta_atm_file=\"${filename_d13c_atm}\"")
	    com=${com}$(csed "cable_user%c13o2_outfile=\"outputs/${mettype}_out_casa_c13o2.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_in_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_out_flux=\"restart/${mettype}_c13o2_flux_rst.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_in_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_out_pools=\"restart/${mettype}_c13o2_pools_rst.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_in_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
	    com=${com}$(csed "cable_user%c13o2_restart_out_luc=\"restart/${mettype}_c13o2_luc_rst.nc\"")
	    if [[ ${c13o2_simple_disc} -eq 1 ]] ; then
		com=${com}$(csed "cable_user%c13o2_simple_disc=.true.")
	    else
		com=${com}$(csed "cable_user%c13o2_simple_disc=.false.")
	    fi
	fi
	sed ${com} ${ndir}/cable_global_${experiment}.nml > ${rdir}/cable.nml

	# run model
	cd ${rdir}
	irm logs/log_cable.txt logs/log_out_cable.txt
	if [[ ${dompi} -eq 1 ]] ; then
            if [[ ${ignore_mpi_err} -eq 1 ]] ; then set +e ; fi
            ${mpiexecdir}mpiexec -n ${nproc} ./${iexe} > logs/log_out_cable.txt
            if [[ ${ignore_mpi_err} -eq 1 ]] ; then set -e ; fi
	else
            ./${iexe} > logs/log_out_cable.txt
	fi
	# save output
	renameid ${rid} cable.nml ${mettype}.nml LUC.nml
	mv *_${rid}.nml restart/
	cd logs
	renameid ${rid} log_cable.txt log_out_cable.txt
	cd ../restart
	copyid ${rid} pop_${mettype}_ini.nc ${mettype}_climate_rst.nc ${mettype}_casa_rst.nc ${mettype}_cable_rst.nc
	if [[ ${doc13o2} -eq 1 ]] ; then copyid ${rid} ${mettype}_c13o2_flux_rst.nc ${mettype}_c13o2_pools_rst.nc ; fi
	cd ../outputs
	renameid ${rid} ${mettype}_out_cable.nc ${mettype}_out_casa.nc
	if [[ ${doc13o2} -eq 1 ]] ; then renameid ${rid} ${mettype}_out_casa_c13o2.nc ; fi
	cd ..
	cd ${pdir}
    done
fi


## if latlon is specified, crop outputs to that region
if [[ ${randompoints} -eq 0 ]] ; then
   echo cropping outputs to specified latlon region
   
   lat1=$(echo ${latlon} | cut -d "," -f 1)
   lat2=$(echo ${latlon} | cut -d "," -f 2)
   lon1=$(echo ${latlon} | cut -d "," -f 3)
   lon2=$(echo ${latlon} | cut -d "," -f 4)

   cd ${rdir}/outputs
   files=$(ls *_cable_*.nc *_LUC_*.nc)  # only for cable and LUC files because casa is not mapped!
   for file in ${files} ; do
      ncks -O -d y,${lat1},${lat2} -d x,${lon1},${lon2} ${file} ${file}  # pretty slow, but working
      ## cdo does not work for 5-dim variables:
      #filename=$(basename -s .nc ${file})
      #cdo -O sellonlatbox,${lon1},${lon2},${lat1},${lat2} ${filename}.nc ${filename}_cropped.nc
      #mv -f ${filename}_cropped.nc ${filename}.nc
   done

   cd ${pdir}
fi
   

# --------------------------------------------------------------------
# Finish
#
cd ${isdir}

t2=$(date +%s)
dt=$((t2-t1))
printf "\n"
if [[ ${dt} -lt 60 ]] ; then
    printf "Finished at %s   in %i seconds.\n" "$(date)" ${dt}
else
    dm=$(echo "(${t2}-${t1})/60." | bc -l)
    printf "Finished at %s   in %.2f minutes.\n" "$(date)" ${dm}
fi

exit 0

# grep -Ei '(filename%restart_in|cable_user%CLIMATE_fromZero|cable_user%YearStart|cable_user%YearEnd|icycle|spincasa|cable_user%CASA_fromZero|cable_user%CASA_DUMP_READ|cable_user%CASA_DUMP_WRITE|cable_user%CASA_OUT_FREQ|cable_user%CASA_SPIN_STARTYEAR|cable_user%CASA_SPIN_ENDYEAR|cable_user%limit_labile|casafile%out|casafile%cnpipool|cable_user%POP_fromZero|cable_user%POP_out|cable_user%POP_restart_in|cable_user%POPLUC|cable_user%POPLUC_RunType|cable_user%LUC_outfile|cable_user%LUC_restart_in|cable_user%LUC_restart_out)' ../driver_files/cable.nml.1900_2017 | pbcopy
