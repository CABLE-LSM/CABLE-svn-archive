module cable_namelists_module

contains

! read namelists and assign to these passed vars, eventually retire
! these vars, or at least declare in this module AND use from here
subroutine read_namelists( filename, vegparmnew, soilparmnew, calcsoilalbedo,  &
                           spinup, delsoilM, delsoilT, output, patchout,       &
                           check, verbose, leaps, logn, fixedCO2,              &
                           spincasainput, spincasa, l_casacnp, l_laiFeedbk,    &
                           l_vcmaxFeedbk, icycle, casafile, ncciy, gswpfile,   &
                           redistrb, wiltParam, satuParam, cable_user )  

   USE cable_IO_vars_module, ONLY: gswp_type,&
                                   output_inclusion_type, checks_type
   USE cable_common_module,  ONLY: kbl_user_switches,     &
                                   filenames_type
    
   USE casavariable,        ONLY: casafiles_type

   implicit none
 
   INTEGER            :: icycle
   TYPE(casafiles_type) :: casafile
   TYPE(checks_type) :: check ! what types of checks to perform
   REAL :: fixedCO2 ! CO2 level if CO2air not in met file
   TYPE(output_inclusion_type)      :: output 
   TYPE(output_inclusion_type)      :: patchout ! do we want patch-specific info

   TYPE(gswp_type)      :: gswpfile
   
   LOGICAL :: leaps   ! use leap year timing?

   INTEGER ::                                                                  &
      ncciy,      & ! year number (& switch) for gswp run
      logn          ! log file unit number
   
   LOGICAL ::                                                                  &
      verbose,    & ! print init and param details of all grid cells?
      soilparmnew   ! read IGBP new soil map. Q.Zhang @ 12/20/2010
  
   LOGICAL :: calcsoilalbedo 
   ! hydraulic_redistribution switch _soilsnow module
   LOGICAL ::                                                                  &
      redistrb ! Turn on/off the hydraulic redistribution
   
   ! hydraulic_redistribution parameters _soilsnow module
   REAL :: wiltParam, satuParam
   
   TYPE(filenames_type) :: filename
   
   TYPE(kbl_user_switches) :: cable_user

   REAL              ::                                                       &  
      delsoilM,        & 
      delsoilT         
  
   LOGICAL ::                                                                 & 
      vegparmnew, spinup, spinConv, spincasainput, spincasa,                  &
      l_casacnp, l_laiFeedbk, l_vcmaxFeedbk
   
   ! CABLE namelist: model configuration, runtime/user switches 
   CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml' 
   CHARACTER(LEN=200), PARAMETER :: NEWCABLE_NAMELIST='newcable.nml' 
   !jhan: Do not need namelist here, can use any namelists and
   ! re-assign to these vars to pass back   
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

   CHARACTER(LEN=200) ::                                                       &
      MetForcing,          & ! Driving Input Data per timestep - site data (obs) 
      CoarseSurfaceInfo,   & ! Defaults to this data if superceded by site data
      VegParameters,       & ! Vegetation parameters file 
      SoilParameters!,      & ! Soil parameters file 
      
   NAMELIST/cable_param_files/      &
      VegParameters,       & ! Vegetation parameters file 
      SoilParameters!,      & ! Soil parameters file 
      
   NAMELIST/offline_driver_files/   &
      MetForcing,          & ! Driving Input Data per timestep - site data (obs) 
      CoarseSurfaceInfo!,   & ! Defaults to this data if superceded by site data

   ! Open, read and close the namelist file.
   OPEN( 10, FILE = CABLE_NAMELIST )
      READ( 10, NML=CABLE )   !where NML=CABLE defined above
   CLOSE(10)
   
   OPEN( 10, FILE = NEWCABLE_NAMELIST )
      READ( 10, NML=cable_param_files)   !where NML=CABLE defined above
      READ( 10, NML=offline_driver_files)   !where NML=CABLE defined above
   CLOSE(10)

      filename%met   = MetForcing
      filename%type  = CoarseSurfaceInfo
      filename%veg   = VegParameters
      filename%soil  = SoilParameters

End subroutine read_namelists

End module cable_namelists_module
