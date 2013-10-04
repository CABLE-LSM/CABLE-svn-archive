
!===COPYRIGHT==================================================================
! The source codes are part of the australian 
! Community Atmosphere Biosphere Land Exchange (CABLE) model. 
! Please register online at xxx and sign the agreement before use 
! contact: whox@xxxx.yyy about registration user agreement
!==============================================================================


!==============================================================================
! Name: cable_driver
! Purpose: offline driver for CABLE model
! CALLed from: executed PROGRAM 
! MODULEs used:   define_dimensions
!                 define_types
!                 io_variables
!                 cable_common_module
!                 input_module
!                 output_module
!                 cbm_module
!                 casadimension
!                 casavariable
! 
! CALLs:       open_met_file
!              load_parameters
!              open_output_file
!              spincasacnp
!              get_met_data
!              casa_feedback
!              cbm
!              bgcdriver
!              sumcflux
!              write_output
!              casa_poolout
!              casa_fluxout
!              create_restart
!              close_met_file
!              close_output_file
!              prepareFiles
!
! Major contribution: land surface modeling team, CSIRO, Aspendale
!
! input  file: ￼[SiteName].nc
!              ￼poolcnpIn[SiteName].csv -- for CASA-CNP only
!              ￼gridinfo_CSIRO_1x1.nc
!              ￼def_veg_params.txt
!              ￼def_soil_params.txt -- nearly redundant, can be switched on
!              ￼restart_in.nc -- not strictly required
!
! output file: ￼log_cable.txt
!              ￼out_cable.nc
!              ￼restart_out.nc
!              ￼poolcnpOut.csv -- from CASA-CNP

!==============================================================================


!==============================================================================
! changes since version release on 
! changes made by who on date
!
!==============================================================================


PROGRAM offline_driver
   USE define_dimensions,   ONLY: ms,mp,mvtype,mstype
   USE define_types
   USE io_variables,        ONLY: logn,filename,gswpfile,ncciy,leaps,          &
                                  verbose, fixedCO2,output,check,patchout,     &
                                  patch_type,soilparmnew,                      &
                                  redistrb, wiltParam, satuParam
   USE cable_common_module, ONLY: ktau_gl, kend_gl, knode_gl, cable_user,      &
                                  cable_runtime
   USE input_module,        ONLY: open_met_file,load_parameters,               &
                                  get_met_data,close_met_file
   USE output_module,       ONLY: create_restart,open_output_file,             &
                                  write_output,close_output_file
   USE cbm_module
   
   ! modules related to CASA-CNP
   USE casadimension,       ONLY: icycle 
   USE casavariable,        ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
                                  casa_met, casa_balance
   USE phenvariable,        ONLY: phen_variable
   USE cable_diag_module,   ONLY: cable_stat, cable_diag 

   IMPLICIT NONE
   
   ! CABLE namelist: model configuration, runtime/user switches 
   CHARACTER(LEN=*), PARAMETER :: CABLE_NAMELIST = 'cable.nml'
   
   ! timing variables 
   INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
   
   INTEGER        ::                                                           &
      ktau,       &  ! increment equates to timestep, resets if spinning up
      ktau_tot,   &  ! NO reset when spinning up, total timesteps by model
      kend,       &  ! no. of time steps in run
      ktauday,    &  ! day counter for CASA-CNP
      idoy,       &  ! day of year (1:365) counter for CASA-CNP
      nyear          ! year counter for CASA-CNP

   REAL :: dels                        ! time step size in seconds
   
   ! CABLE variables
   TYPE (met_type)       :: met     ! met input variables
   TYPE (air_type)       :: air     ! air property variables
   TYPE (canopy_type)    :: canopy  ! vegetation variables
   TYPE (radiation_type) :: rad     ! radiation variables
   TYPE (roughness_type) :: rough   ! roughness varibles
   TYPE (balances_type)  :: bal     ! energy and water balance variables
   TYPE (soil_snow_type) :: ssoil   ! soil and snow variables
   
   ! CABLE parameters
   TYPE (soil_parameter_type) :: soil ! soil parameters	
   TYPE (veg_parameter_type)  :: veg  ! vegetation parameters	 
   
   TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
   TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
   
   ! CASA-CNP variables 
   TYPE (casa_biome)     :: casabiome
   TYPE (casa_pool)      :: casapool
   TYPE (casa_flux)      :: casaflux
   TYPE (casa_met)       :: casamet
   TYPE (casa_balance)   :: casabal
   TYPE (phen_variable)  :: phen 
   
   ! declare vars for switches (default .FALSE.) etc declared thru namelist
   ! jhan: YP* these need looking at b4 i move them to switch status 
   LOGICAL, save           :: &
      vegparmnew = .FALSE.,       & ! using new format input file (BP dec 2007)
      spinup = .FALSE.,           & ! model spinup to soil state equilibrium?
      spinConv = .FALSE.,         & ! has spinup converged?
      spincasainput = .FALSE.,    & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                    ! FALSE: READ input to spin CASA-CNP 
      spincasa = .FALSE.,         & ! TRUE: CASA-CNP Will spin mloop times,
                                    ! FALSE: no spin up
      l_casacnp = .FALSE.,        & ! using CASA-CNP with CABLE
      l_laiFeedbk = .FALSE.,      & ! using prognostic LAI
      l_vcmaxFeedbk = .FALSE.       ! using prognostic Vcmax
   
   
   ! jhan: move them to switch status 
   REAL              :: &  
      delsoilM,         & ! allowed variation in soil moisture for spin up
      delsoilT            ! allowed variation in soil temperature for spin up
  
   ! temporary storage for soil moisture/temp. in spin up mode
   REAL, ALLOCATABLE, DIMENSION(:,:)  :: & 
      soilMtemp,                         &   
      soilTtemp      
   
   ! switches etc defined thru namelist (by default cable.nml)
   ! jhan; comments needed for each, tidy, def type nml% ?
   NAMELIST/CABLE/                  &
                  filename,         & ! TYPE, containing input filenames 
                  vegparmnew,       & ! jhan: use new soil param. method
                  soilparmnew,      & ! jhan: use new soil param. method
                  spinup,           & ! spinup model (soil) to steady state 
                  delsoilM,delsoilT,& ! jhan: TOLS: def TYPE for all tols?
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

   cable_runtime%offline = .TRUE.

   ! if cable_user%CASA_dump_read is false then  the code expects just to run 
   ! CABLE. If you want to write CASA_dump.nc in this case then you have 
   ! to re-set %CASA_dump_write=.true. via namelist 
   ! if %CASA_dump_read = .TRUE. then the code expects to run CASA from 
   ! forcing & you cant write CASA_dump.nc in this case 
   IF ( .NOT. cable_user%CASA_dump_read ) THEN 
   
      IF( l_casacnp  .AND. ( icycle == 0 .OR. icycle > 3 ) )                   &
         STOP 'icycle must be 1 to 3 when using casaCNP'
      IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) .AND. ( .NOT. l_casacnp ) )       &
         STOP 'casaCNP required to get prognostic LAI or Vcmax'
      IF( l_vcmaxFeedbk .AND. icycle < 2 )                                     &
         STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
      IF( icycle > 0 .AND. ( .NOT. soilparmnew ) )                             &
         STOP 'casaCNP must use new soil parameters'
   
      ! Open log file:
      OPEN(logn,FILE=filename%log)
    
      ! Check for gswp run
      IF (ncciy /= 0) THEN
         
         PRINT *, 'Looking for global offline run info.'
         
         IF (ncciy < 1986 .OR. ncciy > 1995) THEN
            PRINT *, 'Year ', ncciy, ' outside range of dataset!'
            STOP 'Please check input in namelist file.'
         ELSE
            
            CALL prepareFiles(ncciy)
         
         ENDIF
      
      ENDIF
   
   ENDIF

   ! Open met data and get site information from netcdf file.
   ! This retrieves time step size, number of timesteps, starting date,
   ! latitudes, longitudes, number of sites. 
   CALL open_met_file( dels, kend, spinup )
 
   ! Checks where parameters and initialisations should be loaded from.
   ! If they can be found in either the met file or restart file, they will 
   ! load from there, with the met file taking precedence. Otherwise, they'll
   ! be chosen from a coarse global grid of veg and soil types, based on 
   ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
   CALL load_parameters( met, air, ssoil, veg, bgc,                            &
                         soil, canopy, rough, rad, sum_flux,                   &
                         bal, logn, vegparmnew, casabiome, casapool,           &
                         casaflux, casamet, casabal, phen )

   IF ( .NOT. cable_user%CASA_dump_read ) THEN 
   
     ! Open output file:
     CALL open_output_file( dels, soil, veg, bgc, rough )
   
     ssoil%otss_0 = ssoil%tgg(:,1)
     ssoil%otss = ssoil%tgg(:,1)
     canopy%fes_cor = 0.
     canopy%fhs_cor = 0.
     met%ofsd = 0.1
   
   ENDIF
  
   ! outer loop = spinup loop:
   DO

      ! globally (WRT code) accessible kend through USE cable_common_module
      kend_gl = kend
      knode_gl = 0
      
      ! time step loop over ktau
      DO ktau=kstart, kend 
         
         ! increment total timstep counter
         ktau_tot = ktau_tot + 1
         
         ! globally (WRT code) accessible kend through USE cable_common_module
         ktau_gl = ktau
         
         ! somethings (e.g. CASA-CNP) only need to be done once per day  
         ktauday=int(24.0*3600.0/dels)
         idoy = mod(ktau/ktauday,365)
         IF(idoy==0) idoy=365
         
         ! needed for CASA-CNP
         nyear =INT((kend-kstart+1)/(365*ktauday))
   
         IF ( .NOT. cable_user%CASA_dump_read ) THEN 
  
            canopy%oldcansto=canopy%cansto
   
            ! Get met data and LAI, set time variables.
            ! Rainfall input may be augmented for spinup purposes:
             met%ofsd = met%fsd(:,1) + met%fsd(:,2)
            CALL get_met_data( spinup, spinConv, met, soil,                    &
                               rad, veg, kend, dels ) 
         ENDIF
   
         IF( cable_user%CASA_dump_read ) then 
            
            IF( mod( (ktau-kstart+1), ktauday ) == 0 ) THEN
               CALL read_casa_dump( casamet, casaflux, idoy, kend )
            ENDIF
         
         ELSE 
            
            ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
            IF (l_vcmaxFeedbk) CALL casa_feedback( ktau_gl, veg, casabiome,    &
                                                   casapool, casamet )
   
            IF (l_laiFeedbk) veg%vlai(:) = casamet%glai(:)
   
            ! CALL land surface scheme for this timestep, all grid points:
            CALL cbm( dels, air, bgc, canopy, met,                             &
                      bal, rad, rough, soil, ssoil,                            &
                      sum_flux, veg )
   
            ssoil%smelt = ssoil%smelt*dels
            ssoil%rnof1 = ssoil%rnof1*dels
            ssoil%rnof2 = ssoil%rnof2*dels
            ssoil%runoff = ssoil%runoff*dels
   
         ENDIF
   
   
         !jhan this is insufficient testing. condition for 
         !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
         IF(icycle >0) THEN
            call bgcdriver( ktau_gl, kstart, kend, dels, met,                  &
                            ssoil, canopy, veg, soil, casabiome,               &
                            casapool, casaflux, casamet, casabal,              &
                            phen, spinConv, spinup, ktauday, idoy,             &
                            cable_user%CASA_dump_read,                         &
                            cable_user%CASA_dump_write )
         ENDIF 
   
         ! sumcflux is pulled out of subroutine cbm
         ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
         CALL sumcflux( ktau_gl, kstart, kend, dels, bgc,                      &
                        canopy, soil, ssoil, sum_flux, veg,                    &
                        met, casaflux, l_vcmaxFeedbk )
   
         IF( .NOT. cable_user%CASA_dump_read ) then 
            ! Write time step's output to file if either: we're not spinning up 
            ! or we're spinning up and the spinup has converged:
            IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
               CALL write_output( dels, met, canopy, ssoil,                    &
                                  rad, bal, air, soil, veg )
         ENDIF
   
         IF (cable_user%CONSISTENCY_CHECK) THEN

            IF((.NOT.spinup).OR.(spinup.AND.spinConv)) then 
               ! temporarily first arg = unique UNIT no. for file 
               ! ALSO temporarily, corresponding scripts on request only 
               call cable_diag( 12121, 1, "Sumfluxes", mp, kend, ktau,         &
                                knode_gl, "sumfluxes", canopy%fe + canopy%fh )       
            ENDIF 
         
         ENDIF 
   
       END DO ! END Do loop over timestep ktau



   
      !jhan this is insufficient testing. condition for 
      !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
      ! see if spinup (if conducting one) has converged:
      IF(spinup.AND..NOT.spinConv) THEN
         
         ! Write to screen and log file:
         WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),      &
               ' of data set complete...'
         WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),   &
               ' of data set complete...'
         
         ! IF not 1st run through whole dataset:
         IF( INT( ktau_tot/kend ) > 1 ) THEN 
            
            ! evaluate spinup
            IF( ANY( ABS(ssoil%wb-soilMtemp)>delsoilM).OR.                     &
                ANY(ABS(ssoil%tgg-soilTtemp)>delsoilT) ) THEN
               
               ! No complete convergence yet
               PRINT *, 'ssoil%wb : ', ssoil%wb
               PRINT *, 'soilMtemp: ', soilMtemp
               PRINT *, 'ssoil%tgg: ', ssoil%tgg
               PRINT *, 'soilTtemp: ', soilTtemp
            
            ELSE ! spinup has converged
               
               spinConv = .TRUE.
               ! Write to screen and log file:
               WRITE(*,'(A33)') ' Spinup has converged - final run'
               WRITE(logn,'(A52)')                                             &
                          ' Spinup has converged - final run - writing all data'
               WRITE(logn,'(A37,F7.5,A28)')                                    &
                          ' Criteria: Change in soil moisture < ',             &
                          delsoilM, ' in any layer over whole run'
               WRITE(logn,'(A40,F7.5,A28)' )                                   &
                          '           Change in soil temperature < ',          &
                          delsoilT, ' in any layer over whole run'
            END IF

         ELSE ! allocate variables for storage
         
           ALLOCATE( soilMtemp(mp,ms), soilTtemp(mp,ms) )
         
         END IF
         
         ! store soil moisture and temperature
         soilTtemp = ssoil%tgg
         soilMtemp = REAL(ssoil%wb)

       ELSE

         ! if not spinning up, or spin up has converged, exit:
         EXIT
       
       END IF

   END DO

   IF (icycle > 0) THEN
      
      CALL casa_poolout( ktau, veg, soil, casabiome,                           &
                         casapool, casaflux, casamet, casabal, phen )

      CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
  
   END IF

   ! Write restart file if requested:
   IF(output%restart)                                                          &
      CALL create_restart( logn, dels, soil, veg, ssoil,                       &
                           canopy, rough, rad, bgc, bal )
      
   ! Close met data input file:
   CALL close_met_file
 
   ! Close output file and deallocate main variables:
   CALL close_output_file( bal, air, bgc, canopy, met,                         &
                           rad, rough, soil, ssoil,                            &
                           sum_flux, veg )

   WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck

   ! Close log file
   CLOSE(logn)

END PROGRAM offline_driver


SUBROUTINE prepareFiles(ncciy)
  USE io_variables, ONLY: logn,gswpfile
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncciy

  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy

  CALL renameFiles(logn,gswpfile%rainf,16,ncciy,'rainf')
  CALL renameFiles(logn,gswpfile%snowf,16,ncciy,'snowf')
  CALL renameFiles(logn,gswpfile%LWdown,16,ncciy,'LWdown')
  CALL renameFiles(logn,gswpfile%SWdown,16,ncciy,'SWdown')
  CALL renameFiles(logn,gswpfile%PSurf,16,ncciy,'PSurf')
  CALL renameFiles(logn,gswpfile%Qair,14,ncciy,'Qair')
  CALL renameFiles(logn,gswpfile%Tair,14,ncciy,'Tair')
  CALL renameFiles(logn,gswpfile%wind,15,ncciy,'wind')

END SUBROUTINE prepareFiles


SUBROUTINE renameFiles(logn,inFile,nn,ncciy,inName)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: logn
  INTEGER, INTENT(IN) :: nn
  INTEGER, INTENT(IN) :: ncciy
  CHARACTER(LEN=99), INTENT(INOUT) :: inFile
  CHARACTER(LEN=*),  INTENT(IN)    :: inName
  INTEGER :: idummy

  READ(inFile(nn:nn+3),'(i4)') idummy
  IF (idummy < 1983 .OR. idummy > 1995) THEN
    PRINT *, 'Check position of the year number in input gswp file', inFile
    STOP
  ELSE
    WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
    WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
  ENDIF

END SUBROUTINE renameFiles
