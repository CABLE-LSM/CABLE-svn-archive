!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: 
!
! Contact: 
!
! History: 
!
!
! ==============================================================================
! Uses:
! 
! CALLs:
!
! input  file: 
!
! output file:
!==============================================================================

MODULE cable_namelists_mod

   ! declare vars for switches (default .FALSE.) re-declared thru namelist
   LOGICAL ::                     &
      vegparmnew,       & ! using new format input file (BP dec 2007)
      spinup,           & ! model spinup to soil state equilibrium?
      spinConv,         & ! has spinup converged?
      spincasainput,    & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                    ! FALSE: READ input to spin CASA-CNP 
      spincasa,         & ! TRUE: CASA-CNP Will spin mloop times,
                                    ! FALSE: no spin up
      l_casacnp,        & ! using CASA-CNP with CABLE
      l_laiFeedbk,      & ! using prognostic LAI
      l_vcmaxFeedbk       ! using prognostic Vcmax
   
CONTAINS

SUBROUTINE cable_namelists
!   USE cable_def_types_mod
   USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
                                   verbose, fixedCO2,output,check,patchout,    &
                                   soilparmnew
   USE cable_common_module,  ONLY: cable_user,     &
                                   filename, redistrb,          & 
                                   wiltParam, satuParam,    &
                                   calcsoilalbedo
!   USE cable_data_module,    ONLY: driver_type, point2constants
!   USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
!                                   get_met_data,close_met_file
!   USE cable_output_module,  ONLY: create_restart,open_output_file,            &
!                                   write_output,close_output_file
!   USE cable_cbm_module
!   
!   USE cable_diag_module
!   
   ! modules related to CASA-CNP
   USE casadimension,       ONLY: icycle 
   USE casavariable,        ONLY: casafile!, casa_biome, casa_pool, casa_flux,  &
!                                  casa_met, casa_balance
!   USE phenvariable,        ONLY: phen_variable
!
!   IMPLICIT NONE
!   
   ! CABLE namelist: model configuration, runtime/user switches 
   CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml' 
!   
!   ! timing variables 
!   INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
!   
!   
!   ! CABLE parameters
!   TYPE (soil_parameter_type) :: soil ! soil parameters	
!   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
!   TYPE (driver_type)    :: C         ! constants used locally  
!   
!   TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
!   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
!   
!   ! CASA-CNP variables 
!   TYPE (casa_biome)     :: casabiome
!   TYPE (casa_pool)      :: casapool
!   TYPE (casa_flux)      :: casaflux
!   TYPE (casa_met)       :: casamet
!   TYPE (casa_balance)   :: casabal
!   TYPE (phen_variable)  :: phen 
  
   REAL              :: &  
      delsoilM,         & ! allowed variation in soil moisture for spin up
      delsoilT            ! allowed variation in soil temperature for spin up
!  
!   ! temporary storage for soil moisture/temp. in spin up mode
!   REAL, ALLOCATABLE, DIMENSION(:,:)  :: & 
!      soilMtemp,                         &   
!      soilTtemp      
!
!   !___ unique unit/file identifiers for cable_diag: arbitrarily 5 here 
!   INTEGER, SAVE :: iDiagZero=0, iDiag1=0, iDiag2=0, iDiag3=0, iDiag4=0
!!
   ! switches etc defined thru namelist (by default cable.nml)
   NAMELIST/CABLE/                  &
                  filename,         & ! TYPE, containing input filenames 
                  vegparmnew,       & ! use new soil param. method
                  soilparmnew,      & ! use new soil param. method
                  calcsoilalbedo,   & ! albedo considers soil color Ticket #27
                  spinup,           & ! spinup model (soil) to steady state 
                  delsoilM,delsoilT,& ! 
                  output,           &
                  patchout,         &
                  check,            &
                  verbose,          &
                  leaps,            &
                  logn,             &
                  fixedCO2,         &
                  spincasainput,    &
                  spincasa,         &
                  l_casacnp,        &
                  l_laiFeedbk,      &
                  l_vcmaxFeedbk,    &
                  icycle,           &
                  casafile,         &
                  ncciy,            &
                  gswpfile,         &
                  redistrb,         &
                  wiltParam,        &
                  satuParam,        &
                  cable_user           ! additional USER switches 

   ! END header

   ! Open, read and close the namelist file.
   OPEN( 10, FILE = CABLE_NAMELIST )
      READ( 10, NML=CABLE )   !where NML=CABLE defined above
   CLOSE(10)

END SUBROUTINE cable_namelists

END MODULE cable_namelists_mod





