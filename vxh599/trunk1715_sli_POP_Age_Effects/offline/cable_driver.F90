!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from
! http://www.accessimulator.org.au/cable
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
! Purpose: Offline driver for CABLE
!
! Contact: Bernard.Pak@csiro.au
!
! History: Since 1.4b, capability to run global offline (ncciy = YEAR),
!          inclusion of call to CASA-CNP (icycle>0)
!          exclusion of call to cbm (icycle>10)
!          soil_snow_type now ssnow (instead of ssoil)
!
!
! ==============================================================================
! Uses:           cable_def_types_mod
!                 cable_IO_vars_module
!                 cable_common_module
!                 cable_input_module
!                 cable_output_module
!                 cable_cbm_module
!                 casadimension
!                 casavariable
!
! CALLs:       open_met_file
!              load_parameters
!              open_output_file
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
!
! input  file: [SiteName].nc
!              poolcnpIn[SiteName].csv -- for CASA-CNP only
!              gridinfo_CSIRO_1x1.nc
!              def_veg_params.txt
!              def_soil_params.txt -- nearly redundant, can be switched on
!              restart_in.nc -- not strictly required
!
! output file: log_cable.txt
!              out_cable.nc
!              restart_out.nc
!              poolcnpOut.csv -- from CASA-CNP
!==============================================================================

PROGRAM cable_offline_driver
  USE cable_def_types_mod
  USE cable_IO_vars_module, ONLY: logn,gswpfile,ncciy,leaps,                  &
       verbose, fixedCO2,output,check,patchout,    &
       patch_type,soilparmnew,&
       defaultLAI
  USE cable_common_module,  ONLY: ktau_gl, kend_gl, knode_gl, cable_user,     &
       cable_runtime, filename, myhome,            &
       redistrb, wiltParam, satuParam, CurYear,    &
       IS_LEAPYEAR
  USE cable_data_module,    ONLY: driver_type, point2constants
  USE cable_input_module,   ONLY: open_met_file,load_parameters,              &
       get_met_data,close_met_file,                &
       ncid_rain,       &
       ncid_snow,       &
       ncid_lw,         &
       ncid_sw,         &
       ncid_ps,         &
       ncid_qa,         &
       ncid_ta,         &
       ncid_wd
  USE cable_output_module,  ONLY: create_restart,open_output_file,            &
       write_output,close_output_file
  USE cable_write_module,   ONLY: nullify_write
  USE cable_cbm_module
  USE cable_diag_module
  ! modules related to CASA-CNP
  USE casadimension,       ONLY: icycle
  USE casavariable,        ONLY: casafile, casa_biome, casa_pool, casa_flux,  &
       casa_met, casa_balance
  USE phenvariable,        ONLY: phen_variable

  ! modules related to POP
  USE POP_Types,     Only: POP_TYPE
  USE POP_Constants, Only: HEIGHT_BINS, NCOHORT_MAX
  !USE DFLIB

  IMPLICIT NONE

  ! CABLE namelist: model configuration, runtime/user switches
  CHARACTER(LEN=200), PARAMETER :: CABLE_NAMELIST='cable.nml'

  ! timing variables
  INTEGER, PARAMETER ::  kstart = 1   ! start of simulation
  INTEGER, PARAMETER ::  mloop  = 100   ! CASA-CNP PreSpinup loops

  INTEGER        ::                                                           &
       ktau,       &  ! increment equates to timestep, resets if spinning up
       ktau_tot,   &  ! NO reset when spinning up, total timesteps by model
       kend,       &  ! no. of time steps in run
                                !CLN      kstart = 1, &  ! timestep to start at
       koffset = 0, &  ! timestep to start at
       ktauday,    &  ! day counter for CASA-CNP
       idoy,       &  ! day of year (1:365) counter for CASA-CNP
       nyear,      &  ! year counter for CASA-CNP
       casa_it,    &  ! number of calls to CASA-CNP
       YYYY,       &  !
       RYEAR,      &  !
       RRRR,       &  !
       NRRRR,      &  !
       ctime,      &  ! day count for casacnp
       LOY
  REAL :: dels                        ! time step size in seconds

  INTEGER,DIMENSION(:,:),ALLOCATABLE :: GSWP_MID
  CHARACTER     :: dum*9

  ! CABLE variables
  TYPE (met_type)       :: met     ! met input variables
  TYPE (air_type)       :: air     ! air property variables
  TYPE (canopy_type)    :: canopy  ! vegetation variables
  TYPE (radiation_type) :: rad     ! radiation variables
  TYPE (roughness_type) :: rough   ! roughness varibles
  TYPE (balances_type)  :: bal     ! energy and water balance variables
  TYPE (soil_snow_type) :: ssnow   ! soil and snow variables

  ! CABLE parameters
  TYPE (soil_parameter_type) :: soil ! soil parameters
  TYPE (veg_parameter_type)  :: veg  ! vegetation parameters
  TYPE (driver_type)    :: C         ! constants used locally

  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables

  ! CASA-CNP variables
  TYPE (casa_biome)     :: casabiome
  TYPE (casa_pool)      :: casapool
  TYPE (casa_flux)      :: casaflux
  TYPE (casa_met)       :: casamet
  TYPE (casa_balance)   :: casabal
  TYPE (phen_variable)  :: phen
  TYPE ( POP_TYPE )     :: POP
  CHARACTER             :: cyear*4
  CHARACTER             :: ncfile*99

  ! declare vars for switches (default .FALSE.) etc declared thru namelist
  LOGICAL, SAVE           :: &
       vegparmnew     = .FALSE.,       & ! using new format input file (BP dec 2007)
       spinup = .FALSE.,           & ! model spinup to soil state equilibrium?
       spinConv = .FALSE.,         & ! has spinup converged?
       spincasainput = .FALSE.,    & ! TRUE: SAVE input req'd to spin CASA-CNP;
                                ! FALSE: READ input to spin CASA-CNP
       spincasa = .FALSE.,         & ! TRUE: CASA-CNP Will spin mloop times,
                                ! FALSE: no spin up
       l_casacnp = .FALSE.,        & ! using CASA-CNP with CABLE
       l_laiFeedbk = .FALSE.,      & ! using prognostic LAI
       l_vcmaxFeedbk = .FALSE.,    & ! using prognostic Vcmax
       CASAONLY      = .FALSE.,    & ! ONLY Run CASA-CNP
       CALL1 = .TRUE.,             &
       SPINon= .TRUE.


  REAL              :: &
       delsoilM,         & ! allowed variation in soil moisture for spin up
       delsoilT            ! allowed variation in soil temperature for spin up

  ! temporary storage for soil moisture/temp. in spin up mode
  REAL, ALLOCATABLE, DIMENSION(:,:)  :: &
       soilMtemp,                         &
       soilTtemp

  ! timing
  REAL:: etime ! Declare the type of etime(), For receiving user and system time, total time
  REAL, ALLOCATABLE, DIMENSION(:) :: cleaf_max, npp_ann, stemnpp_ann, gpp_ann, &
       leafnpp_ann, gpp_ann_save ! variables for keeping track of maximum annual leaf carbon, annual npp, annual gpp

  ! switches etc defined thru namelist (by default cable.nml)
  NAMELIST/CABLE/                  &
       filename,         & ! TYPE, containing input filenames
       vegparmnew,       & ! jhan: use new soil param. method
       soilparmnew,      & ! jhan: use new soil param. method
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

  REAL :: SUMME
  INTEGER :: i,x
  !INTEGER,dimension(:), ALLOCATABLE :: ALLVEG
  ! END header

  ! Open, read and close the namelist file.
  OPEN( 10, FILE = CABLE_NAMELIST )
  READ( 10, NML=CABLE )   !where NML=CABLE defined above
  CLOSE(10)

  ! Open log file:
  OPEN(logn,FILE=filename%log)

!!!! INISTUFF

  CurYear = CABLE_USER%YearStart

  IF ( icycle .GE. 11 ) THEN
     icycle                     = icycle - 10
     CASAONLY                   = .TRUE.
     CABLE_USER%CASA_DUMP_READ  = .TRUE.
     CABLE_USER%CASA_DUMP_WRITE = .FALSE.
  ELSEIF ( icycle .EQ. 0 ) THEN
     CABLE_USER%CASA_DUMP_READ  = .FALSE.
     spincasa                   = .FALSE.
     CABLE_USER%CALL_POP        = .FALSE.
  ENDIF

  IF ( .NOT. spinup ) THEN
     IF ( spincasa ) THEN
        spincasa = .FALSE.
        WRITE(*,*)   "spinup == .FALSE. -> spincasa set to .F."
        WRITE(logn,*)"spinup == .FALSE. -> spincasa set to .F."
     ENDIF
  ENDIF

  IF ( TRIM(cable_user%MetType) .EQ. 'gswp' ) THEN
     leaps = .FALSE.
     WRITE(*,*)   "gswp data doesn't have leap years!!! leaps -> .FALSE."
     WRITE(logn,*)"gswp data doesn't have leap years!!! leaps -> .FALSE."
  ENDIF

  IF ( TRIM(cable_user%MetType) .EQ. 'gpgs' ) THEN
     leaps = .TRUE.
     cable_user%MetType = 'gswp'
  ENDIF

  cable_runtime%offline = .TRUE.

  ! associate pointers used locally with global definitions
  CALL point2constants( C )

  IF( l_casacnp  .AND. ( icycle == 0 .OR. icycle > 3 ) )                   &
       STOP 'icycle must be 1 to 3 when using casaCNP'
  !IF( ( l_laiFeedbk .OR. l_vcmaxFeedbk ) )       &
  !   STOP 'casaCNP required to get prognostic LAI or Vcmax'
  IF( l_vcmaxFeedbk .AND. icycle < 2 )                                     &
       STOP 'icycle must be 2 to 3 to get prognostic Vcmax'
  IF( icycle > 0 .AND. ( .NOT. soilparmnew ) )                             &
       STOP 'casaCNP must use new soil parameters'

  NRRRR = merge(MAX(CABLE_USER%CASA_NREP,1), 1, CASAONLY)
  ctime = 0

!!!! INISTUFF

  ! Open met data and get site information from netcdf file. (NON-GSWP ONLY!)
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites.
  IF ( TRIM(cable_user%MetType) .NE. "gswp" .AND. &
       TRIM(cable_user%MetType) .NE. "gpgs") THEN

     CALL open_met_file( dels, koffset, kend, spinup, C%TFRZ )
     IF ( koffset .NE. 0 .AND. CABLE_USER%CALL_POP ) THEN
        WRITE(*,*)"When using POP, episode must start at Jan 1st!"
        STOP 991
     ENDIF
  ELSE IF ( NRRRR .GT. 1 ) THEN
     ALLOCATE( GSWP_MID( 8, CABLE_USER%YearStart:CABLE_USER%YearEnd ) )
  ENDIF

  ! outer loop - spinup loop no. ktau_tot :
  ktau_tot = 0
  SPINon   = .TRUE.

  SPINLOOP:DO WHILE ( SPINon )

     NREP: DO RRRR = 1, NRRRR
        YEAR: DO YYYY= CABLE_USER%YearStart,  CABLE_USER%YearEnd
           CurYear = YYYY
           IF ( leaps .and. IS_LEAPYEAR( YYYY ) ) THEN
              LOY = 366
           ELSE
              LOY = 365
           ENDIF
           ! Check for gswp run
           IF ( TRIM(cable_user%MetType) .EQ. 'gswp' ) THEN
              ncciy = CurYear

              IF (ncciy < 1986 .OR. ncciy > 1995) THEN
                 PRINT *, 'Year ', ncciy, ' outside range of dataset!'
                 STOP 'Please check input in namelist file.'
              ELSE IF ( RRRR .EQ. 1 ) THEN
                 PRINT *, 'Looking for global offline run info.'
                 CALL prepareFiles(ncciy)
              ENDIF
              IF ( RRRR .EQ. 1 ) THEN
                 CALL open_met_file( dels, koffset, kend, spinup, C%TFRZ )
                 IF ( NRRRR .GT. 1 ) THEN
                    GSWP_MID(1,YYYY) = ncid_rain
                    GSWP_MID(2,YYYY) = ncid_snow
                    GSWP_MID(3,YYYY) = ncid_lw
                    GSWP_MID(4,YYYY) = ncid_sw
                    GSWP_MID(5,YYYY) = ncid_ps
                    GSWP_MID(6,YYYY) = ncid_qa
                    GSWP_MID(7,YYYY) = ncid_ta
                    GSWP_MID(8,YYYY) = ncid_wd
                 ENDIF
              ELSE
                 PRINT *, 'Use previously loaded met for ',YYYY
                 ncid_rain = GSWP_MID(1,YYYY)
                 ncid_snow = GSWP_MID(2,YYYY)
                 ncid_lw   = GSWP_MID(3,YYYY)
                 ncid_sw   = GSWP_MID(4,YYYY)
                 ncid_ps   = GSWP_MID(5,YYYY)
                 ncid_qa   = GSWP_MID(6,YYYY)
                 ncid_ta   = GSWP_MID(7,YYYY)
                 ncid_wd   = GSWP_MID(8,YYYY)
                 kend      = ktauday * LOY
              ENDIF
           ENDIF

           ! Checks where parameters and initialisations should be loaded from.
           ! If they can be found in either the met file or restart file, they will
           ! load from there, with the met file taking precedence. Otherwise, they'll
           ! be chosen from a coarse global grid of veg and soil types, based on
           ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
           IF ( CALL1 ) THEN
              CALL load_parameters( met, air, ssnow, veg, bgc,           &
                   soil, canopy, rough, rad, sum_flux,                   &
                   bal, logn, vegparmnew, casabiome, casapool,           &
                   casaflux, casamet, casabal, phen, POP, spinup,        &
                   C%EMSOIL, C%TFRZ )

              ! Open output file:
              IF (.not.CASAONLY) THEN
                 IF ( TRIM(filename%out) .EQ. '' ) THEN
                    IF ( CABLE_USER%YEARSTART .GT. 0 ) THEN
                       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART, &
                            CABLE_USER%YEAREND
                       filename%out = TRIM(filename%path)//'/'//&
                            TRIM(cable_user%RunIden)//'_'//&
                            TRIM(dum)//'_cable_out.nc'
                    ELSE
                       filename%out = TRIM(filename%path)//'/'//&
                            TRIM(cable_user%RunIden)//'_cable_out.nc'
                    ENDIF
                 ENDIF
                 IF (RRRR.eq.1) then
                    call nullify_write() ! nullify pointers
                    CALL open_output_file( dels, soil, veg, bgc, rough )
                 endif
              ENDIF

              ssnow%otss_0 = ssnow%tgg(:,1)
              ssnow%otss = ssnow%tgg(:,1)
              canopy%fes_cor = 0.
              canopy%fhs_cor = 0.
              met%ofsd = 0.1

              spinConv = .FALSE. ! initialise spinup convergence variable
              if(.not.spinup)  spinConv=.true.
              IF (.NOT.ALLOCATED(gpp_ann_save)) ALLOCATE(  gpp_ann_save(mp) )
              gpp_ann_save = -999.

              if( icycle>0 .AND. spincasa) then
                 print *, 'EXT spincasacnp enabled with mloop= ', mloop
                 ! CALL read_casa_dump(casafile%dump_cnpspin, casamet, casaflux, kstart, kend)
                 call spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                      casaflux,casamet,casabal,phen,gpp_ann_save)
              endif


              IF (.NOT.ALLOCATED(cleaf_max)) ALLOCATE(  cleaf_max(mp) )
              IF (.NOT.ALLOCATED(npp_ann)) ALLOCATE(  npp_ann(mp) )
              IF (.NOT.ALLOCATED(stemnpp_ann)) ALLOCATE(  stemnpp_ann(mp) )
              IF (.NOT.ALLOCATED(leafnpp_ann)) ALLOCATE(  leafnpp_ann(mp) )
              IF (.NOT.ALLOCATED(gpp_ann)) ALLOCATE(  gpp_ann(mp) )

              cleaf_max = 0.0
              npp_ann = 0.0
              gpp_ann = 0.0
              stemnpp_ann = 0.0
              leafnpp_ann = 0.0


           ENDIF

           ! globally (WRT code) accessible kend through USE cable_common_module
           kend_gl = kend
           knode_gl = 0

           ! somethings (e.g. CASA-CNP) only need to be done once per day
           ktauday=int(24.0*3600.0/dels)

           ! time step loop over ktau
           DO ktau=kstart, kend
              ! increment total timstep counter
              ktau_tot = ktau_tot + 1

              ! globally (WRT code) accessible kend through USE cable_common_module
              ktau_gl = ktau_tot

              idoy = mod(INT(REAL(ktau+koffset)/REAL(ktauday)),LOY)
              IF ( idoy .EQ. 0 ) idoy = LOY

              ! needed for CASA-CNP
              nyear     =INT((kend+koffset)/(LOY*ktauday))
              canopy%oldcansto=canopy%cansto

              ! Get met data and LAI, set time variables.
              ! Rainfall input may be augmented for spinup purposes:
              CALL get_met_data( spinup, spinConv, met, soil,                    &
                   rad, veg, kend, dels, C%TFRZ, ktau+koffset,     &
                   kstart+koffset )

              IF ( TRIM(cable_user%MetType) .NE. 'gswp' ) CurYear = met%year(1)
              met%ofsd = met%fsd(:,1) + met%fsd(:,2)
              ! Zero out lai where there is no vegetation acc. to veg. index
              WHERE ( veg%iveg(:) .GE. 14 ) veg%vlai = 0.



              IF ( .NOT. CASAONLY ) THEN

                 ! Feedback prognostic vcmax and daily LAI from casaCNP to CABLE
                 IF (l_vcmaxFeedbk) CALL casa_feedback( ktau, veg, casabiome,    &
                      casapool, casamet )

                 IF (l_laiFeedbk) veg%vlai(:) = casamet%glai(:)

                 ! CALL land surface scheme for this timestep, all grid points:

                 CALL cbm( ktau, dels, air, bgc, canopy, met,                             &
                      bal, rad, rough, soil, ssnow,                            &
                      sum_flux, veg )


                 ssnow%smelt = ssnow%smelt*dels
                 ssnow%rnof1 = ssnow%rnof1*dels
                 ssnow%rnof2 = ssnow%rnof2*dels
                 ssnow%runoff = ssnow%runoff*dels

              ELSE IF ( mod((ktau-kstart+1+koffset),ktauday)==0 .AND. CABLE_USER%CASA_DUMP_READ ) THEN
                 ! CLN READ FROM FILE INSTEAD !
                 WRITE(CYEAR,FMT="(I4)")CurYear + INT((ktau-kstart+koffset)/(LOY*ktauday))
                 ncfile  = TRIM(casafile%c2cdumppath)//'c2c_'//CYEAR//'_dump.nc'
                 casa_it = NINT( REAL(ktau / ktauday) )
              
                 CALL read_casa_dump( ncfile, casamet, casaflux, casa_it, kend, .FALSE. )       
                         
              ENDIF
              !jhan this is insufficient testing. condition for
              !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
              IF(icycle >0 .OR.  CABLE_USER%CASA_DUMP_WRITE ) THEN


                 call bgcdriver( ktau, kstart, kend, dels, met,                         &
                      ssnow, canopy, veg, soil, casabiome,                   &
                      casapool, casaflux, casamet, casabal,                  &
                      phen, pop, spinConv, spinup, ktauday, idoy, loy,           &
                      CABLE_USER%CASA_DUMP_READ, CABLE_USER%CASA_DUMP_WRITE, &
                      gpp_ann_save )

                 IF(MOD((ktau-kstart+1),ktauday)==0) THEN
                    gpp_ann = gpp_ann + casaflux%cgpp
                    npp_ann = npp_ann + casaflux%cnpp
                    stemnpp_ann = stemnpp_ann + casaflux%cnpp*casaflux%fracCalloc(:,2)*0.7
                    leafnpp_ann = leafnpp_ann + casaflux%cnpp*casaflux%fracCalloc(:,1)
                 ENDIF


                 ! temp variable for max annual cleaf

                 WHERE (cleaf_max.lt.casapool%cplant(:,1))
                    cleaf_max = casapool%cplant(:,1)
                 ENDWHERE
              ENDIF
              ! WRITE CASA OUTPUT
              IF(icycle >0) THEN
                ! IF((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN
                    IF ( MOD ((ktau+koffset),ktauday*CABLE_USER%CASA_OUT_FREQ) == 0  &
                         .OR. ktau .EQ. kend) THEN
                       !ctime = (ktau + ( YYYY - cable_user%YearStart ) * 365 * ktauday)/ktauday
                       ctime = ctime +1
                  
 
                    CALL WRITE_CASA_OUTPUT_NC ( casamet, casapool, casabal, casaflux, &
                            CASAONLY, ctime, ( ktau.EQ.kend .AND. YYYY .EQ.               &
                           cable_user%YearEnd.AND. RRRR .EQ.NRRRR ) )
                    ENDIF
                ! END IF
              ENDIF

              IF ( .NOT. CASAONLY ) THEN

                 ! sumcflux is pulled out of subroutine cbm
                 ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
                 CALL sumcflux( ktau, kstart, kend, dels, bgc,                      &
                      canopy, soil, ssnow, sum_flux, veg,                    &
                      met, casaflux, l_vcmaxFeedbk )

              ENDIF

              ! Write time step's output to file if either: we're not spinning up
              ! or we're spinning up and the spinup has converged:

              IF ( .NOT. CASAONLY ) THEN

                 CALL write_output( dels, ktau_tot, met, canopy, ssnow,                    &
                      rad, bal, air, soil, veg, C%SBOLTZ, &
                      C%EMLEAF, C%EMSOIL )

              ENDIF
              ! dump bitwise reproducible testing data
              IF( cable_user%RUN_DIAG_LEVEL == 'zero') THEN
                 IF (.NOT.CASAONLY) THEN
                   write(*,*) 'before diags'
                    IF((.NOT.spinup).OR.(spinup.AND.spinConv))                         &
                      call cable_diag( 1, "FLUXES", mp, kend, ktau,                   &
                      knode_gl, "FLUXES",                            &
                      canopy%fe + canopy%fh )
                 ENDIF
              ENDIF
           END DO ! END Do loop over timestep ktau

           !jhan this is insufficient testing. condition for
           !spinup=.false. & we want CASA_dump.nc (spinConv=.true.)
           ! see if spinup (if conducting one) has converged:
           !IF(spinup.AND..NOT.spinConv) THEN
           IF(spinup.AND.(.NOT.spinConv).and.(.NOT.CASAONLY) ) THEN

              ! Write to screen and log file:
              WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),      &
                   ' of data set complete...'
              WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(ktau_tot/kend),   &
                   ' of data set complete...'

              ! IF not 1st run through whole dataset:
              !         IF( INT( ktau_tot/kend ) > 1 .AND. RRRR.EQ. NRRRR) THEN
              ! ORI         IF( INT( ktau_tot/kend ) > 1 ) THEN
              IF( MOD( ktau_tot, kend ) .EQ. 0 .AND. ktau_Tot .GT. kend .AND. &
                   YYYY.EQ. CABLE_USER%YearEnd .OR. ( NRRRR .GT. 1 .AND. &
                   RRRR.EQ. NRRRR) ) THEN
                 !         IF( INT( ktau_tot/kend ) > 1 ) THEN

                 ! evaluate spinup
                 IF( ANY( ABS(ssnow%wb-soilMtemp)>delsoilM).OR.                     &
                      ANY( ABS(ssnow%tgg-soilTtemp)>delsoilT) ) THEN

                    ! No complete convergence yet
                    PRINT *, 'ssnow%wb : ', ssnow%wb
                    PRINT *, 'soilMtemp: ', soilMtemp
                    PRINT *, 'ssnow%tgg: ', ssnow%tgg
                    PRINT *, 'soilTtemp: ', soilTtemp

                 ELSE ! spinup has converged

                    spinConv = .TRUE.
                    ! Write to screen and log file:
                    WRITE(*,'(A33)') ' Spinup has converged - final run'
                    WRITE(logn,'(A52)')                                             &
                         ' Spinup has converged - final run - writing all data'
                    WRITE(logn,'(A37,F8.5,A28)')                                    &
                         ' Criteria: Change in soil moisture < ',             &
                         delsoilM, ' in any layer over whole run'
                    WRITE(logn,'(A40,F8.5,A28)' )                                   &
                         '           Change in soil temperature < ',          &
                         delsoilT, ' in any layer over whole run'
                 END IF

              ELSE ! allocate variables for storage

                 IF (.NOT.ALLOCATED(soilMtemp)) ALLOCATE(  soilMtemp(mp,ms) )
                 IF (.NOT.ALLOCATED(soilTtemp)) ALLOCATE(  soilTtemp(mp,ms) )

              END IF

              ! store soil moisture and temperature
              IF ( YYYY.EQ. CABLE_USER%YearEnd ) THEN
                 soilTtemp = ssnow%tgg
                 soilMtemp = REAL(ssnow%wb)
              ENDIF

           ELSE

              ! if not spinning up, or spin up has converged, exit:
              IF ( SpinOn ) THEN
                 PRINT*,"setting SPINON -> FALSE", YYYY, RRRR
                 SPINon = .FALSE.
              END IF

           END IF

           CALL1 = .FALSE.

           !IF ( .NOT. spinup ) THEN
           IF((.NOT.spinup).OR.(spinup.AND.spinConv)) THEN
              IF (icycle > 0) THEN

                 CALL casa_poolout( ktau, veg, soil, casabiome,                           &
                      casapool, casaflux, casamet, casabal, phen )
                 CALL casa_fluxout( nyear, veg, soil, casabal, casamet)
                 CALL write_casa_restart_nc ( casamet, casapool, met, CASAONLY )

              END IF

              !         IF ( .NOT. CASAONLY .AND. CABLE_USER%YEARSTART .GT. 0 ) THEN
              IF ( .NOT. CASAONLY ) THEN
                 ! Write restart file if requested:
                 IF(output%restart)                                                          &
                      CALL create_restart( logn, dels, ktau, soil, veg, ssnow,                 &
                      canopy, rough, rad, bgc, bal, met )

              ELSE
                 ! TESTART OUTPUT FOR CASA
              ENDIF
           ENDIF
           IF ( .NOT. spinup .OR. spinconv ) THEN
              if ( NRRRR .GT. 1 ) THEN
                 RYEAR = YYYY + ( CABLE_USER%YearEnd - CABLE_USER%YearStart + 1 ) &
                      * ( RRRR - 1 )
              ELSE
                 RYEAR = YYYY
              END if
              IF ( cable_user%CALL_POP ) then

                 IF (TRIM(cable_user%POP_out).eq.'epi') THEN
                    CALL POP_IO( pop, casamet, RYEAR, 'WRITE_EPI', &
                      (YYYY.EQ.CABLE_USER%YearEnd .AND. RRRR.EQ.NRRRR) )
                    CALL GLOBFOR_OUT(mp, pop, casapool, veg, rad, cleaf_max, npp_ann, gpp_ann, stemnpp_ann, &
                      leafnpp_ann )

                 ENDIF
                 
              endif
              cleaf_max = 0.0
              npp_ann = 0.0
              gpp_ann_save = gpp_ann
              gpp_ann = 0.0
              stemnpp_ann = 0.0
              leafnpp_ann = 0.0
           ENDIF
           ! Close met data input file:
           IF ( TRIM(cable_user%MetType) .EQ. "gswp" .AND. &
                RRRR .EQ. NRRRR ) THEN
              CALL close_met_file
              IF ( YYYY .EQ. CABLE_USER%YearEnd .AND. &
                   NRRRR .GT. 1 ) DEALLOCATE ( GSWP_MID )
           ENDIF
        ! re-initalise annual flux sums
        casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
        casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0

        END DO YEAR

     END DO NREP

  END DO SPINLOOP

  IF ( SpinConv .AND. .NOT. CASAONLY ) THEN
     ! Close output file and deallocate main variables:
     CALL close_output_file( bal, air, bgc, canopy, met,                         &
          rad, rough, soil, ssnow,                            &
          sum_flux, veg )
  ENDIF

  IF ( cable_user%CALL_POP ) THEN
     IF ( CASAONLY ) THEN
        CALL POP_IO( pop, casamet, RYEAR-1, 'WRITE_INI', .TRUE.)
     ELSE
        CALL POP_IO( pop, casamet, RYEAR-1, 'WRITE_RST', .TRUE.)
     ENDIF
  ENDIF

  IF ( TRIM(cable_user%MetType) .NE. "gswp" ) CALL close_met_file

  WRITE(logn,*) bal%wbal_tot, bal%ebal_tot, bal%ebal_tot_cncheck


  ! Close log file
  CLOSE(logn)
  call cpu_time(etime)
  print *, 'Finished. ', etime, ' seconds needed for ', kend,' hours'

END PROGRAM cable_offline_driver


SUBROUTINE prepareFiles(ncciy)
  USE cable_IO_vars_module, ONLY: logn,gswpfile
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncciy

  WRITE(logn,*) 'CABLE offline global run using gswp forcing for ', ncciy
  PRINT *,      'CABLE offline global run using gswp forcing for ', ncciy

  CALL renameFiles(logn,gswpfile%rainf,ncciy,'rainf')
  CALL renameFiles(logn,gswpfile%snowf,ncciy,'snowf')
  CALL renameFiles(logn,gswpfile%LWdown,ncciy,'LWdown')
  CALL renameFiles(logn,gswpfile%SWdown,ncciy,'SWdown')
  CALL renameFiles(logn,gswpfile%PSurf,ncciy,'PSurf')
  CALL renameFiles(logn,gswpfile%Qair,ncciy,'Qair')
  CALL renameFiles(logn,gswpfile%Tair,ncciy,'Tair')
  CALL renameFiles(logn,gswpfile%wind,ncciy,'wind')

END SUBROUTINE prepareFiles



SUBROUTINE renameFiles(logn,inFile,ncciy,inName)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: logn,ncciy
  INTEGER:: nn
  CHARACTER(LEN=200), INTENT(INOUT) :: inFile
  CHARACTER(LEN=*),  INTENT(IN)    :: inName
  INTEGER :: idummy

  nn = INDEX(inFile,'19')
  READ(inFile(nn:nn+3),'(i4)') idummy
  IF (idummy < 1983 .OR. idummy > 1995) THEN
     PRINT *, 'Check position of the year number in input gswp file', inFile
     STOP
  ELSE
     WRITE(inFile(nn:nn+3),'(i4.4)') ncciy
     WRITE(logn,*) TRIM(inName), ' global data from ', TRIM(inFile)
  ENDIF

END SUBROUTINE renameFiles






