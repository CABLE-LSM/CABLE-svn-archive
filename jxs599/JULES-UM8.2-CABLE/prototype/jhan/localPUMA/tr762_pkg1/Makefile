###############################################################################
##
## This is the top-level makefile for building a JULES executable.
##
## This makefile is maintained to be compliant with the GNU make utility ONLY, 
## although it will most probably work with other make utilities.  The GNU make
## utility is available free-of-charge from http://www.gnu.org/software/make/, 
## but most Unix/Linux systems have it lurking somewhere anyway.  (Type
## 'make -v' or 'which make' to check.)
##
## To build, you will need to set the following three variables to point to 
## the directory in which you are building JULES (most probably the one 
## containing this makefile) and the location of the pre-compiled netCDF 
## Fortran interface.  The netCDF interface library is available free-of-charge
## from Unidata at http://www.unidata.ucar.edu/software/netcdf/.  Note that
## this makefile does NOT build the netCDF library, that must be done 
## beforehand as a separate process.
##
## CDF_LIB_PATH should point to the directory containing the netCDF library
## archives (e.g. libnetcdff.a and libnetcdf.a).
## CDF_MOD_PATH should point to the directory containing the netCDF Fortran 90
## module file (e.g. netcdf.mod).
##
## There are two options for the build that can be passed to the makefile from
## the command line.  The variable COMPILER allows the user to tell the make 
## which compiler they are using and the makefile will pick up compiler-specific
## variables from the separate text files, e.g. Makefile.comp.sun, that come 
## with this distribution. 
##
## The variable BUILD allows the user to maintain several JULES libraries and 
## executables built with different compiler flags.  This is particularly useful
## when making any modifications to the JULES source code, as it allows the user
## to switch to a debugging version of the model without having to recompile
## all of the source code or edit any makefiles.  This makefile is set up to 
## build a "normal" version, a debugging version and an optimised/fast version.
## Currently, there are 3 allowed values for BUILD: run, debug, fast.
## 
## These two options can be specified from the command line, e.g.,
## 
##    make BUILD=run COMPILER=sun
##
## KNOWN BUGS.
## * Using Intel compiler on Linux, GNU make sometimes returns root "/" from 
##   $(PWD) command.  Get around this by specifying JULESDIR=/home/jules on 
##   the command line.
## * Using Intel compiler on Linux, GNU make does not move the .mod files to 
##   MOD_DIR after each sub-make completes.  Get around this by adding 
##   -module $(JULESDIR)/MODS to FF_* variables in Makefile.comp.intel.
##
###############################################################################
JULESDIR=$(PWD)


## The user is advised not to modify code below this line. ##


# We carry FPP_VARS all the way through so we can add useful bits to it as we go
FPP_VARS=

######################################
## Make some decisions about the    ##
## netCDF library.                  ##
######################################
CDFDUMMY=false
OBJ_CDF_DUMMY=
ifeq ($(CDFDUMMY),true)
DIR_CDF_DUMMY=$(JULESDIR)/utils/netcdf_dummy
CDF_MOD_PATH=$(DIR_CDF_DUMMY)
CDF_LIB_PATH=$(DIR_CDF_DUMMY)
OBJ_CDF_DUMMY=$(DIR_CDF_DUMMY)/jules_netcdf_dummy.o
CLEAN_DEPENDENCIES += clean_cdfdummy
# Add CDF_DUMMY as a preprocessor var
FPP_VARS += NCDF_DUMMY
endif

ifndef CDF_LIB_PATH
$(error You need to set the variable CDF_LIB_PATH either in the JULES Makefile or from the command line)
endif

ifndef CDF_MOD_PATH
$(error You need to set the variable CDF_MOD_PATH either in the JULES Makefile or from the command line)
endif

#######################################
## Set up variables for drhook dummy ##
#######################################
DIR_DRHOOK_DUMMY=$(JULESDIR)/utils/drhook_dummy
DRHOOK_OBJECTS=$(DIR_DRHOOK_DUMMY)/parkind1.o \
               $(DIR_DRHOOK_DUMMY)/yomhook.o
CLEAN_DEPENDENCIES += clean_drhook_dummy

################################################################
## Set up variables for mpi dummy library                     ##
##                                                            ##
## Compiling with MPI is only supported with FCM make         ##
## http://www.metoffice.gov.uk/research/collaboration/fcm     ##
################################################################
DIR_MPI_DUMMY=$(JULESDIR)/utils/mpi_dummy
MPI_OBJECTS=$(DIR_MPI_DUMMY)/mpi.o
CLEAN_DEPENDENCIES += clean_mpi_dummy
FPP_VARS += MPI_DUMMY

######################################
## Rules for making library archive ##
## members from source files.       ##
######################################
include $(JULESDIR)/Makefile.common.mk

######################################
## Choose compiler.                 ##
######################################
## Set default values for these variables else you can get 
## some nasty "rm -f *" commands issued during a clean.
FC=f90
CC=cc
LIB_INC=-L
LIB_PRE=-l
LIB_FPRE=lib
LIB_FSUF=.a
LINKER=$(FC)
MOD_INC=-M
MOD_FSUF=.mod

COMPILER=
ifndef COMPILER
$(error You need to set the variable COMPILER either in the JULES Makefile or from the command line)
endif

ifeq ($(COMPILER),sun)
include $(JULESDIR)/Makefile.comp.sun
FPP_VARS += COMPILER_SUN
endif

ifeq ($(COMPILER),g95)
include $(JULESDIR)/Makefile.comp.g95
FPP_VARS += COMPILER_G95
endif

ifeq ($(COMPILER),gfortran)
include $(JULESDIR)/Makefile.comp.gfortran
FPP_VARS += COMPILER_GFORTRAN
endif

ifeq ($(COMPILER),nag)
include $(JULESDIR)/Makefile.comp.nag
FPP_VARS += COMPILER_NAG
endif

ifeq ($(COMPILER),intel)
include $(JULESDIR)/Makefile.comp.intel
FPP_VARS += COMPILER_INTEL
endif

ifeq ($(COMPILER),pgf)
include $(JULESDIR)/Makefile.comp.pgf
FPP_VARS += COMPILER_PGF
endif

ifeq ($(COMPILER),misc)
include $(JULESDIR)/Makefile.comp.misc
FPP_VARS += COMPILER_MISC
endif

ifeq ($(COMPILER)$(FC),misc@@)
$(error You can only use COMPILER misc if you have set the variables in Makefile.comp.misc)
endif


######################################
## Construct an number of useful    ##
## variables based on user/compiler ##
## options.                         ##
######################################
EXENAME=jules$(LIB_SPECIAL).exe
JULESLIBNOM=$(LIB_FPRE)jules$(LIB_SPECIAL)$(LIB_FSUF)
JULESLIB=$(JULESDIR)/$(JULESLIBNOM)
ARC=$(JULESLIB)

MAINOBJNOM=jules.o
MAINOBJ=$(JULESDIR)/$(MAINOBJNOM)
HEADMAKE=$(JULESDIR)/Makefile $(JULESDIR)/Makefile.common.mk
MOD_DIR=$(JULESDIR)/mods
INCLUDE_DIR=$(JULESDIR)/includes

LIB_PATH  = $(LIB_INC)$(JULESDIR) 
LIB_PATH += $(LIB_INC)$(CDF_LIB_PATH) 

LIB_LIBS  = $(LIB_PRE)jules$(LIB_SPECIAL)
ifeq ($(CDFDUMMY),false)
LIB_LIBS += $(LIB_PRE)netcdff $(LIB_PRE)netcdf 
endif

MOD_PATH  = $(MOD_PUT)$(MOD_DIR)
MOD_PATH += $(MOD_INC)$(MOD_DIR)
MOD_PATH += $(MOD_INC)$(CDF_MOD_PATH)
MOD_PATH += $(MOD_INC)$(DIR_DRHOOK_DUMMY)

CFLAGS=
FFLAGS=$(FF_RUN)
LDFLAGS=$(LIB_PATH) $(LIB_LIBS)

######################################
## Start of compilation/library     ##
## options.                         ##
######################################
# Options are run, debug, fast, condor.
BUILD=run
LIB_SPECIAL=

ifeq ($(BUILD),debug)
  LIB_SPECIAL=_debug
  FFLAGS=$(FF_DBG)
  LINKER=$(FC)
endif

ifeq ($(BUILD),fast)
  LIB_SPECIAL=_fast
  FFLAGS=$(FF_FAST)
  LINKER=$(FC)
endif

ifeq ($(BUILD),condor)
  LIB_SPECIAL=_condor
  FFLAGS=$(FF_CON)
  LINKER=condor_compile $(FC)
endif

################################################
## Build pre-processor variable definitions   ##
################################################
FPP_VARS += SCMA BL_DIAG_HACK
# Add the preproccesor define flag to each specifed var
FPP_DEFS = $(foreach FPP_VAR,$(FPP_VARS),$(FPP_FDEF)$(FPP_VAR))

#################################################
## Set up the paths for pre-processor includes ##
#################################################
FPP_INC_PATH = $(FPP_INC)$(INCLUDE_DIR)/shared 
FPP_INC_PATH += $(FPP_INC)$(INCLUDE_DIR)/standalone 
FPP_INC_PATH += $(FPP_INC)$(DIR_MPI_DUMMY)


#######################################
## Make these variables available to ##
## the sub-makes.                    ##
#######################################
export FC
export FPP
export FPP_DEFS
export FPP_INC_PATH
export MOD_FSUF
export MOD_DIR
export CFLAGS
export FFLAGS
export MOD_PATH
export LIB_PATH
export LIB_LIBS
export JULESDIR
export ARC

######################################
## Variables for executing the      ##
## recursive build and clean.       ##
######################################
SRC_DIRS = src/params/standalone                                \
           src/params/imogen                                    \
           src/control/shared                                   \
           src/science/params                                   \
           src/util                                             \
           src/util/cube                                        \
           src/util/grid                                        \
           src/io/file_handling/core/drivers/ascii              \
           src/io/file_handling/core/drivers/ncdf               \
           src/io/file_handling/core                            \
           src/io/file_handling/gridded                         \
           src/io/file_handling/timeseries                      \
           src/control/standalone/var                           \
           src/control/imogen/var                               \
           src/io/model_interface                               \
           src/io/input                                         \
           src/io/input/time_varying/interpolation              \
           src/io/input/time_varying                            \
           src/io/output                                        \
           src/control/standalone/parallel                      \
           src/io/dump                                          \
           src/initialisation/standalone/ancillaries            \
           src/initialisation/standalone/grid                   \
           src/initialisation/standalone/initial_conditions     \
           src/initialisation/standalone/params                 \
           src/control/standalone/spinup                        \
           src/control/standalone/update                        \
           src/initialisation/standalone                        \
           src/control/imogen                                   \
           src/control/standalone                               \
           src/science/radiation                                \
           src/science/snow                                     \
           src/science/soil                                     \
           src/science/surface                                  \
           src/science/vegetation

CLEAN_SRC_DIRS=$(SRC_DIRS:%=clean_%)
MAKE_SRC_DIRS=$(SRC_DIRS:%=make_%)

CLEAN_DEPENDENCIES += $(CLEAN_SRC_DIRS)

######################################
## Dependencies for building library##
## and executable.                  ##
######################################

all : $(EXENAME)
.PHONY : all

$(EXENAME) : $(JULESLIB) $(MAINOBJ) $(OBJ_CDF_DUMMY) mpi drhook $(HEADMAKE)
	$(LINKER) $(FFLAGS) -o $(EXENAME) $(MAINOBJ) \
	    $(OBJ_CDF_DUMMY) \
      $(DRHOOK_OBJECTS) \
      $(MPI_OBJECTS)    \
	    $(LIB_PATH) \
	    $(MOD_PATH) \
	    $(LIB_LIBS)

$(JULESLIB) : $(OBJ_CDF_DUMMY) mpi drhook $(MAKE_SRC_DIRS)
	touch $(JULESLIB)

$(MAINOBJ) : $(JULESLIB)
	$(EXTRACT) $(JULESLIB) $(MAINOBJNOM)

$(OBJ_CDF_DUMMY) : $(HEADMAKE)
	@$(MAKE) -C $(DIR_CDF_DUMMY)
	
mpi : $(HEADMAKE)
	@$(MAKE) -C $(DIR_MPI_DUMMY)

drhook : $(HEADMAKE)
	@$(MAKE) -C $(DIR_DRHOOK_DUMMY)

$(MAKE_SRC_DIRS) : make_% : $(HEADMAKE)
	@$(MAKE) -C $*

######################################
## Dependencies for cleaning.       ##
######################################
.PHONY : clean clean_cdfdummy clean_mpi_dummy clean_drhook_dummy
clean: $(CLEAN_DEPENDENCIES)
	$(RM) $(JULESLIB) $(EXENAME) $(wildcard $(MOD_DIR)/*$(MOD_FSUF)) $(MAINOBJ)
$(CLEAN_SRC_DIRS) : clean_% :
	@$(MAKE) clean -C $*
clean_cdfdummy :
	@$(MAKE) clean -C $(DIR_CDF_DUMMY)
clean_mpi_dummy :
	@$(MAKE) clean -C $(DIR_MPI_DUMMY)
clean_drhook_dummy :
	@$(MAKE) clean -C $(DIR_DRHOOK_DUMMY)

## End of file.
