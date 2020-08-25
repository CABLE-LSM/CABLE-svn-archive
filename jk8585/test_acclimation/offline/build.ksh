#!/bin/ksh

export dosvn=1 # 1/0: do/do not check svn

# so that script can be called by 'bash build.ksh' if no ksh installed
if [ -z ${PS3} ] ; then # https://stackoverflow.com/questions/3327013/how-to-determine-the-current-shell-im-working-on
    eval 'function print(){ printf "$@\n"; }'
fi

known_hosts()
{
    if [ -z ${PS3} ] ; then
	kh=(pear jigg ccrc mael mcin vm_o gadi)
    else
        set -A kh pear jigg ccrc mael mcin vm_o gadi
    fi
}


host_mael()
{
   export NCDIR='/share/apps/netcdf/intel/4.1.3/lib'
   export NCMOD='/share/apps/netcdf/intel/4.1.3/include'
   export FC=ifort
   export CFLAGS='-O2 -fp-model precise -fpp'
   export CFLAGS='-O3 -fp-model precise  -ipo --parallel '   
   #export LDFLAGS='-L/share/apps/intel/Composer/lib/intel64 -L/share/apps/netcdf/intel/4.1.3/lib -O2'
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -traceback -g -fp-model precise -fpp' 
      export LDFLAGS='-L/share/apps/intel/Composer/lib/intel64 -L/share/apps/netcdf/intel/4.1.3/lib'
   fi
   export LD='-lnetcdf -lnetcdff'
   build_build
   cd ../
   build_status
}


host_ccrc()
{
   export NCDIR='/usr/local/netcdf/intel/4.1.3/lib'
   export NCMOD='/usr/local/netcdf/intel/4.1.3/include'
   export FC=ifort
   export CFLAGS='-O2 -fp-model precise -fpp'   #-traceback
   if [[ $1 = 'debug' ]]; then
      export CFLAGS='-O0 -debug -g -ftrapuv -CB -check bounds -diag-enable warn -fpp'
# -diag-enable sc2 -diag-enable sc-single-file
   fi
   export LD='-lnetcdf -lnetcdff'
   export LDFLAGS='-L/usr/local/intel/Compiler/11.1/lib/intel64 -L//usr/local/netcdf/intel/4.1.3/lib -O2'
   if [[ $1 = 'debug' ]]; then
      export LDFLAGS='-L/usr/local/intel/Compiler/11.1/lib/intel64 -L//usr/local/netcdf/intel/4.1.3/lib -O0 -debug -g -ftrapuv -diag-enable warn'
# -diag-enable sc2 -diag-enable sc-single-file
   fi
   build_build
   cd ../
   build_status
}


## jiggle
host_jigg()
{
   export NCDIR='/usr/local/lib'
   export NCMOD='/usr/local/include'
   export FC=gfortran
   export CFLAGS='-O2 -x f95-cpp-input'
   export LD='-lnetcdf -lnetcdff'
   export LDFLAGS='-L/usr/local/lib -O2'
   build_build
   cd ../
   build_status
}


## pearcey.hpsc.csiro.au 
host_pear()
{
   # if [ ! -z ${PS3} ] ; then
   . /apps/modules/Modules/default/init/ksh
   # fi

   module del intel-cc intel-fc
   module add intel-cc/16.0.1.150 intel-fc/16.0.1.150
   module add netcdf/4.3.3.1

   export NCDIR=$NETCDF_ROOT'/lib/'
   export NCMOD=$NETCDF_ROOT'/include/'
   export FC='ifort'
   export  CFLAGS='-O0 -fp-model precise -fpe0 -fpp -g -debug -traceback -fp-stack-check -no-ftz -ftrapuv -check all,noarg_temp_created -C '
   #export CFLAGS='-O0 -fpe=0 -fpe-all=0 -fpp -g -debug -traceback -fp-stack-check -no-ftz -ftrapuv -check bounds 
   #export CFLAGS='-O2 -fp-model precise -fpp'
   export CFLAGS="${CFLAGS} -D__CRU2017__"
   export CFLAGS="${CFLAGS} -D__NETCDF3__"
   export LDFLAGS='-g -L'$NCDIR' -O0'
   export MFLAGS='-j 8'
   export LD='-lnetcdf -lnetcdff'
   export dosvn=0
   export MFLAGS='-j 8'
   build_build
   cd ../
   build_status
}


# MatthiasCuntz@INRA
host_mcin()
{
    idebug=0
    iintel=1
    ignu=0
    inag=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
	if [[ "${1}" == "debug" ]] ; then
	    idebug=1
	    shift 1
	elif [[ "${1}" == "ifort" || "${1}" == "intel" ]] ; then
	    iintel=1
	    ignu=0
	    inag=0
	    shift 1
	elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            iintel=0
	    ignu=1
	    inag=0
	    shift 1
        elif [[ "${1}" == "nagfor" || "${1}" == "nag" ]] ; then
            iintel=0
	    ignu=0
	    inag=1
            shift 1
	else
	    echo "Error: command line option not known: " ${1}
	    exit 1
	fi
    done
    if [[ ${iintel} -eq 1 ]] ;  then
        # INTEL
        /opt/intel/compilers_and_libraries/mac/bin/compilervars.sh intel64
        export FC=ifort
        # release
	export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382 -xHost"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
	    export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv -init=arrays,snan
        fi
	# export CFLAGS="${CFLAGS} -mtune=corei7"
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=''
        export NCROOT='/usr/local/netcdf-fortran-4.4.5-ifort'
    elif [[ ${ignu} -eq 1 ]] ;  then
        # GFORTRAN
        export FC=gfortran
        # release
        export CFLAGS="-cpp -O3 -Wno-aggressive-loop-optimizations -ffree-form -ffixed-line-length-132"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-cpp -O -g -pedantic-errors -Wall -W -Wno-maybe-uninitialized -ffree-form -ffixed-line-length-132 -fbacktrace -ffpe-trap=zero,overflow -finit-real=nan" #  -ffpe-trap=zero,overflow,underflow
        fi
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__"
	export LD=''
        export NCROOT='/usr/local/netcdf-fortran-4.4.5-gfortran'
    elif [[ ${inag} -eq 1 ]] ;  then
        # NAG
        export FC=nagfor
        # release
        export CFLAGS="-O4"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-C -C=dangling -g -nan -O0 -strict95 -gline"
        fi
	export CFLAGS="${CFLAGS} -fpp -colour -unsharedf95 -kind=byte -ideclient -ieee=full -free -not_openmp"
	# export CFLAGS="${CFLAGS} -march=native"
	export CFLAGS="${CFLAGS} -D__NAG__"
        export LD='-ideclient -unsharedrts'
        export NCROOT='/usr/local/netcdf-fortran-4.4.5-nagfor'
    fi
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    # All compilers
    export NCCROOT='/usr/local'
    export NCCLIB=${NCCROOT}'/lib'
    export NCLIB=${NCROOT}'/lib'
    export NCMOD=${NCROOT}'/include'
    export LDFLAGS="-L${NCLIB} -L${NCCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz"
    export dosvn=0
    export MFLAGS='-j 8'
    build_build
    cd ../
    build_status
}


# MatthiasCuntz@Explor
host_vm_o()
{
    idebug=0
    iintel=1
    ignu=0
    np=$#
    for ((i=0; i<${np}; i++)) ; do
	if [[ "${1}" == "debug" ]] ; then
	    idebug=1
	    shift 1
	elif [[ "${1}" == "ifort" || "${1}" == "intel" ]] ; then
	    iintel=1
	    ignu=0
	    shift 1
	elif [[ "${1}" == "gfortran" || "${1}" == "gnu" ]] ; then
            iintel=0
	    ignu=1
	    shift 1
	else
	    echo "Error: command line option not known: " ${1}
	    exit 1
	fi
    done
    if [[ ${iintel} -eq 1 ]] ;  then
        # INTEL
        module load intel/2018.5
        export FC=ifort
        # release
        export CFLAGS="-O3 -fpp -nofixed -assume byterecl -fp-model precise -ip -xHost -diag-disable=10382"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-check all,noarg_temp_created -warn all -g -debug -traceback -fp-stack-check -O0 -debug -fpp -nofixed -assume byterecl -fp-model precise -ip -xHost -diag-disable=10382"
        fi
	# export CFLAGS="${CFLAGS} -march=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -march=core-avx2"     # std / hf
	# export CFLAGS="${CFLAGS} -mtune=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -march=skylake-avx512 # sky
	# export CFLAGS="${CFLAGS} -march=ivybridge"     # ivy / k20
	# export CFLAGS="${CFLAGS} -march=avx"           # ivy / k20
	# export CFLAGS="${CFLAGS} -mtune=ivybridge"     # ivy / k20
	export CFLAGS="${CFLAGS} -D__INTEL__ -D__INTEL_COMPILER__"
        export LD=''
        export NCROOT='/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-ifort2018.0'
    else
        # GFORTRAN
        module load gcc/6.3.0
        export FC=gfortran
        # release
        export CFLAGS="-O3 -Wno-aggressive-loop-optimizations -cpp -ffree-form -ffixed-line-length-132"
        if [[ ${idebug} -eq 1 ]] ; then
            # debug
            export CFLAGS="-pedantic-errors -Wall -W -O -g -Wno-maybe-uninitialized -cpp -ffree-form -ffixed-line-length-132"
        fi
	# export CFLAGS="${CFLAGS} -march=broadwell"     # std / hf
	# export CFLAGS="${CFLAGS} -mavx2"               # std / hf
	# export CFLAGS="${CFLAGS} -march=skylake-avx512 # sky
	# export CFLAGS="${CFLAGS} -march=ivybridge"     # ivy / k20
	# export CFLAGS="${CFLAGS} -mavx"                # ivy / k20
	export CFLAGS="${CFLAGS} -D__GFORTRAN__ -D__gFortran__" # -pg # gprof --line ./cable gmon.out
        export LD=''
        export NCROOT='/home/oqx29/zzy20/local/netcdf-fortran-4.4.4-gfortran63'
    fi
    # export CFLAGS="${CFLAGS} -D__C13DEBUG__"
    export CFLAGS="${CFLAGS} -D__CRU2017__"
    export CFLAGS="${CFLAGS} -D__NETCDF3__"

    # All compilers
    export NCCROOT='/home/oqx29/zzy20/local'
    export NCCLIB=${NCCROOT}'/lib'
    export NCLIB=${NCROOT}'/lib'
    export NCMOD=${NCROOT}'/include'
    export LDFLAGS="-L${NCCLIB} -L${NCLIB} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz" # -pg
    export dosvn=0
    export MFLAGS='-j 8'
    build_build
    cd ../
    build_status
}


## gadi.nci.org.au
host_gadi()
{
   if [ -z ${PS3} ] ; then
      . /etc/bashrc
   else
      . /etc/kshrc
   fi
   module purge
   module add intel-compiler/2019.5.281
   module add netcdf/4.6.3

   export FC='ifort'
   export NCDIR=${NETCDF_ROOT}'/lib/Intel'
   export NCMOD=${NETCDF_ROOT}'/include/Intel'
   #if [[ ${1} == 'debug' ]]; then
       # export CFLAGS='-O0 -fpp -traceback -g -fp-model precise -ftz -fpe0'
       #export CFLAGS="-fpp -O0 -debug extended -traceback -g -check all,noarg_temp_created -warn all -fp-stack-check -nofixed -assume byterecl -fp-model precise -diag-disable=10382 -fpe0" # -fpe-all=0 -no-ftz -ftrapuv"
       #export LDFLAGS='-O0'
   #else
       # export CFLAGS='-O2 -fpp -fp-model precise'
       export CFLAGS="-fpp -O3 -nofixed -assume byterecl -fp-model precise -ip -diag-disable=10382"
       export LDFLAGS='-O3'
       # export CFLAGS="${CFLAGS} -xCORE-AVX2 -axSKYLAKE-AVX512,CASCADELAKE" # given in user training: does not work
       export CFLAGS="${CFLAGS} -xCASCADELAKE" # or -xCORE-AVX512 # queues: express / normal
       # export CFLAGS="${CFLAGS} -xBROADWELL"  # or -xCORE-AVX512; queues: expressbw / normalbw
       # export CFLAGS="${CFLAGS} -xSKYLAKE"    # or -xSKYLAKE-AVX512 depends on performance; queues: normalsl
   #fi
   export CFLAGS="${CFLAGS} -D__CRU2017__"
   export CFLAGS="${CFLAGS} -D__NETCDF3__"
   export LDFLAGS='-L'${NCDIR}' '${LDFLAGS}
   export LD='-lnetcdf -lnetcdff'
   export MFLAGS='-j 8'
   build_build
   cd ../
   build_status
}


## unknown machine, user entering options stdout 
host_read()
{
   print "\n\tWhat is the ROOT path of your NetCDF library" \
         "and .mod file. "
   print "\tRemember these have to be created by the same " \
         "Fortran compiler you" 
   print "\twant to use to build CABLE. e.g./usr/local/intel"
   read NCDF_ROOT
   
   print "\n\tWhat is the path, relative to the above ROOT, of " \
         "your NetCDF library." 
   print "\n\tPress enter for default [lib]."
   read NCDF_DIR
   if [[ $NCDF_DIR == '' ]]; then
      export NCDIR=$NCDF_ROOT/'lib'
   else   
      export NCDIR=$NCDF_ROOT/$NCDF_DIR
   fi

   print "\n\tWhat is the path, relative to the above ROOT, of " \
         "your NetCDF .mod file."
   print "\n\tPress enter for default [include]."
   read NCDF_MOD
   if [[ $NCDF_MOD == '' ]]; then
      export NCMOD=$NCDF_ROOT/'include'
   else   
      export NCMOD=$NCDF_ROOT/$NCDF_MOD
   fi

   print "\n\tWhat is the Fortran compiler you wish to use."
   print "\te.g. ifort, gfortran"
   
   print "\n\tPress enter for default [ifort]."
   read FCRESPONSE 
   if [[ $FCRESPONSE == '' ]]; then
      export FC='ifort'
   else   
      export FC=$FCRESPONSE
   fi

   print "\n\tWhat are the approriate compiler options"
   print "\te.g.(ifort) -O2 -fp-model precise "
   print "\n\tPress enter for default [-O2 -fp-model precise]."
   read CFLAGRESPONSE 
   if [[ $CFLAGRESPONSE == '' ]]; then
      export CFLAGS='-O2 -fp-model precise'
   else   
      export CFLAGS=$CFLAGRESPONSE
   fi

   iflags='-L'$NCDIR' -O2'
   export LDFLAGS=$iflags

   print "\n\tWhat are the approriate libraries to link"
   print "\te.g.(most systems) -lnetcdf "
   print "\n\tPress enter for default [-lnetcdf]."
   read LDRESPONSE 
   if [[ $LDRESPONSE == '' ]]; then
      export LD='-lnetcdf'
   else   
      export LD=$LDRESPONSE
   fi
}


host_write()
{
   print '#!/bin/ksh' > junk
   print '' >> junk
   print 'known_hosts()' >> junk
   print '{' >> junk
   print '   set -A kh' ${kh[*]} $HOST_MACH >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
   print '## '$HOST_COMM >> junk
   print 'host_'$HOST_MACH'()' >> junk
   print '{' >> junk
   print '   export NCDIR='"'"$NCDIR"'" >> junk
   print '   export NCMOD='"'"$NCMOD"'" >> junk
   print '   export FC='$FC >> junk
   print '   export CFLAGS='"'"$CFLAGS"'" >> junk
   print '   export LD='"'"$LD"'" >> junk
   print '   export LDFLAGS='"'"$LDFLAGS"'" >> junk
   print '   build_build' >> junk
   print '   cd ../' >> junk
   print '   build_status' >> junk
   print '}' >> junk
   print '' >> junk
   print '' >> junk
}


clean_build()
{
      print '\ncleaning up\n'
      print '\n\tPress Enter too continue buiding, Control-C to abort now.\n'
      read dummy
      rm -fr .tmp
}


not_recognized()
{  
   print "\n\n\tThis is not a recognized host for which we " \
         "know the location of the" 
   print "\tnetcdf distribution and correct compiler switches."

   print "\n\tPlease enter these details as prompted, and the " \
         "script will be " 
   print "\tupdated accordingly. " 
   print "\n\tIf this is a common machine for CABLE users, " \
         "please email"
   print "\n\t\t cable_help@nf.nci.org.au "  
   print "\n\talong with your new build.ksh so that we can " \
         "update the script "
   print "\tfor all users. "
   print "\n\tTo enter compile options for this build press " \
         "enter, otherwise " 
   print "\tControl-C to abort script."           
   
   host_read

   print "\n\tPlease supply a comment include the new build " \
         "script." 
   print "\n\tGenerally the host URL e.g. gadi.nci.org.au "
   read HOST_COMM
   
   build_build
}


do_i_no_u()
{
    if [ -z ${PS3} ] ; then
       kmax=${#kh[*]}
       k=0
   else
       integer kmax=${#kh[*]}
       integer k=0
   fi
   typeset -f subr

   while [[ $k -lt $kmax ]]; do
      if [[ $HOST_MACH = ${kh[$k]} ]];then
         echo 'Host recognized as' $HOST_MACH
         subr=host_${kh[$k]}
         $subr $*
      fi
      (( k = k + 1 ))
   done
}


build_status()
{
   if [[ -f .tmp/cable ]]; then
        mv .tmp/cable .
        print '\nBUILD OK\n'
   else
      print '\nOooops. Something went wrong\n'        
      print '\nKnown build issues:\n'        
      print '\nSome systems require additional library. \n'        
      print '\nEdit Makefile_offline; add -lnetcdff to LD = ...\n'        
   fi
   exit
}

      
i_do_now()
{
      cd ../
      host_write
      tail -n +7 build.ksh > build.ksh.tmp
      cat junk build.ksh.tmp > build.ksh.new
      mv build.ksh.new build.ksh
      chmod u+x build.ksh 
      rm -f build.ksh.tmp build.ksh.new junk 
      build_status
}


build_build()
{
   if [[ ${dosvn} -eq 1 ]] ; then
       # write file for consumption by Fortran code
       # get SVN revision number 
       CABLE_REV=`svn info | grep Revis |cut -c 11-18`
       if [[ $CABLE_REV = "" ]]; then
           echo "this is not an svn checkout"
           CABLE_REV=0
           echo "setting CABLE revision number to " $CABLE_REV 
       fi         
       echo $CABLE_REV > ~/.cable_rev
       # get SVN status 
       CABLE_STAT=`svn status`
       echo $CABLE_STAT >> ~/.cable_rev
   fi
 
   if [[ ! -d .tmp ]]; then
      mkdir .tmp
   fi
   
   if [[ -f cable ]]; then
      print '\ncable executable exists. copying to dated backup file\n' 
      mv cable cable.`date +%d.%m.%y`
   fi
   
   # directories contain source code
   PHYS="../core/biogeophys"
   UTIL="../core/utils"
   DRV="."
   CASA="../core/biogeochem"
   BLAZE="../core/blaze"
   
   /bin/cp -p $PHYS/*90  ./.tmp
   /bin/cp -p $UTIL/*90  ./.tmp
   /bin/cp -p $DRV/*90   ./.tmp
   /bin/cp -p $CASA/*90  ./.tmp
   /bin/cp -p $BLAZE/*90 ./.tmp
    
   /bin/cp -p Makefile_offline ./.tmp
   
   cd .tmp/
   make -f Makefile_offline ${MFLAGS}
}


###########################################
## build.ksh - MAIN SCRIPT STARTS HERE   ##
###########################################

if [[ $1 = 'clean' ]]; then
    clean_build
    shift 1
fi

known_hosts
HOST_MACH=`uname -n | cut -c 1-4 | tr - _`
do_i_no_u $*

# only ksh because host_write writes ksh script
not_recognized
i_do_now