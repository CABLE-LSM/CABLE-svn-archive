#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_imogen.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_imogen()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, IMOGEN_UNIT, NAMELIST_UNIT

  USE datetime_mod, ONLY : SECS_IN_DAY

  USE string_utils_mod, ONLY : to_string

  USE model_time_mod, ONLY : main_run_start, main_run_end, is_spinup,         &
                             timestep_len

  USE dump_mod, ONLY : read_dump

  USE switches, ONLY : l_imogen, l_360

  USE ancil_info, ONLY : land_pts

  USE imogen_constants

  USE imogen_map

  USE imogen_time

  USE imogen_run

  USE imogen_anlg_vals

  USE imogen_clim

  USE imogen_drive_vars

  USE imogen_io_vars

  USE imogen_progs

  USE trifctl, ONLY : cv

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises IMOGEN parameters and arrays
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_clim  ! File to read IMOGEN
                                                 ! climatology from

  INTEGER :: i,j,l,n  ! Index variables

  INTEGER :: error, error_sum  ! Error indicator


!-----------------------------------------------------------------------------


  IF ( .NOT. l_imogen ) RETURN


! Open the IMOGEN namelist file
  OPEN(NAMELIST_UNIT, FILE='imogen.nml',                                      &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_imogen",                                             &
                   "Error opening namelist file imogen.nml " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! There are three namelists to read from this file
  CALL log_info("init_imogen", "Reading IMOGEN_RUN_LIST namelist...")
  READ(NAMELIST_UNIT, nml=imogen_run_list, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_imogen",                                             &
                   "Error reading namelist IMOGEN_RUN_LIST " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_imogen", "Reading IMOGEN_ANLG_VALS_LIST namelist...")
  READ(NAMELIST_UNIT, nml=imogen_anlg_vals_list, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_imogen",                                             &
                   "Error reading namelist IMOGEN_ANLG_VALS_LIST " //         &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_imogen",                                             &
                   "Error closing namelist file imogen.nml " //               &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")



!-----------------------------------------------------------------------------
! Check that the setup of JULES time variables is compatible with IMOGEN and
! set up IMOGEN time variables
!-----------------------------------------------------------------------------
! -> Force 360 day year
  IF ( .NOT. l_360 )                                                          &
    CALL log_fatal("init_imogen", "360 day year must be used with IMOGEN")

! -> Check no spinup has been requested (a full IMOGEN experiment consists of
!    3 JULES runs with the first 2 acting as a really long spinup)
  IF ( is_spinup )                                                            &
    CALL log_fatal("init_imogen", "IMOGEN runs do not require spinup")

! -> Run must start at 00:00 on 1st Jan for some year
!    Since we know there is no spinup, we can check this just by checking
!    the main run start time
  IF ( main_run_start%time /= 0 .OR. main_run_start%day /= 1 .OR.             &
                                     main_run_start%month /= 1 )              &
    CALL log_fatal("init_imogen",                                             &
                   "IMOGEN runs must start at 00:00 on 1st Jan for some year")
! Store the start year in the IMOGEN variables
  YEAR1   = main_run_start%year
  IYSTART = main_run_start%year
  IYEAR   = main_run_start%year

! -> Set the IMOGEN end year variable - although it doesn't make
!    much sense to have an IMOGEN run ending in the middle of
!    a year, it won't mess anything up, so we allow it
  IYEND = main_run_end%year

! -> Set steps per day and check it is not more than the maximum number of
!    timesteps per day allowed by IMOGEN
!    Note that we already know that SECS_IN_DAY MOD timestep_len is 0
!    since this is enforced in init_time
  STEP_DAY = SECS_IN_DAY / timestep_len
  IF(STEP_DAY > NSDMAX)                                                       &
    CALL log_fatal("init_imogen", "Too many timesteps per day")

!-----------------------------------------------------------------------------
! Find corresponding land sites from the imogen grid 'sgind'.
!-----------------------------------------------------------------------------
  CALL GET_IMOGEN_MAP(FILE_POINTS_ORDER)

!-----------------------------------------------------------------------------
! Allocate imogen arrays
!-----------------------------------------------------------------------------
  error_sum = 0
  ALLOCATE(T_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(RAINFALL_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SNOWFALL_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(RH15M_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(UWIND_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(VWIND_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(DTEMP_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(PSTAR_HA_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SW_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(LW_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(F_WET_CLIM(land_pts,mm), STAT=error )
  error_sum = error_sum + error

  ALLOCATE(T_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(CONV_RAIN_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(CONV_SNOW_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(LS_RAIN_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(LS_SNOW_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(QHUM_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(WIND_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(PSTAR_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(SW_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(LW_OUT(land_pts,mm,md,nsdmax), STAT=error )
  error_sum = error_sum + error

  ALLOCATE(LAT(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(LONG(land_pts), STAT=error )
  error_sum = error_sum + error
  ALLOCATE(DCTOT(land_pts), STAT=error )
  error_sum = error_sum + error

! Check for error.
  IF ( error_sum /= 0 ) THEN
    CALL log_fatal("init_imogen", "Error allocating IMOGEN arrays")
  ELSE
! Initialise.
    T_CLIM(:,:)        = 0.0
    RAINFALL_CLIM(:,:) = 0.0
    SNOWFALL_CLIM(:,:) = 0.0
    RH15M_CLIM(:,:)    = 0.0
    UWIND_CLIM(:,:)    = 0.0
    VWIND_CLIM(:,:)    = 0.0
    DTEMP_CLIM(:,:)    = 0.0
    PSTAR_HA_CLIM(:,:) = 0.0
    SW_CLIM(:,:)       = 0.0
    LW_CLIM(:,:)       = 0.0
    F_WET_CLIM(:,:)    = 0.0

    T_OUT(:,:,:,:)         = 0.0
    CONV_RAIN_OUT(:,:,:,:) = 0.0
    CONV_SNOW_OUT(:,:,:,:) = 0.0
    LS_RAIN_OUT(:,:,:,:)   = 0.0
    LS_SNOW_OUT(:,:,:,:)   = 0.0
    QHUM_OUT(:,:,:,:)      = 0.0
    WIND_OUT(:,:,:,:)      = 0.0
    PSTAR_OUT(:,:,:,:)     = 0.0
    SW_OUT(:,:,:,:)        = 0.0
    LW_OUT(:,:,:,:)        = 0.0

    LAT(:)   = 0.0
    LONG(:)  = 0.0
    DCTOT(:) = 0.0
  END IF

! Weather generator is not available at present
  IF ( WGEN )                                                                 &
    CALL log_fatal("init_imogen", "Weather generator not available at present")

!-----------------------------------------------------------------------------
! Check that the configurations proposed are valid.
! This contains the runs that are allowed on the "decision" tree given under
! directory "plots" and called "imogen.jpg"
!-----------------------------------------------------------------------------
  CALL IMOGEN_CHECK(                                                          &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED, OCEAN_FEED,ANLG,ANOM   &
  )

!-----------------------------------------------------------------------------
! This subroutine only allows runs that have been checked.
! Other runs can be made (editing out the "STOP" command),
! but at one's peril!
!-----------------------------------------------------------------------------
  CALL IMOGEN_CONFIRMED_RUN(                                                  &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,                        &
    OCEAN_FEED,.TRUE.,WGEN,ANLG,ANOM                                          &
  )

!-----------------------------------------------------------------------------
! Now get the emissions/CO2 concentrations.
! At present only coded for reading in a file of emission/CO2 concentrat
! imogen projects of the future climate, with carbon cycle feedbacks) or
! in a file of CO2 concentrations for "Hydrology 20th Century" simulatio
!-----------------------------------------------------------------------------
  IF( C_EMISSIONS .AND. INCLUDE_CO2 .AND. ANOM .AND. ANLG ) THEN
    OPEN(IMOGEN_UNIT, FILE=FILE_SCEN_EMITS,                                   &
                      STATUS='old', POSITION='rewind', ACTION='read')

    DO N = 1,NYR_EMISS
      READ(IMOGEN_UNIT,*) YR_EMISS(N),C_EMISS(N)
    ENDDO

    CLOSE(IMOGEN_UNIT)
  ENDIF

!-----------------------------------------------------------------------
! Read in monthly control climate data.
!-----------------------------------------------------------------------
  DO J = 1,MM
    FILE_CLIM = TRIM(DIR_CLIM) // DRIVE_MONTH(J)

    OPEN(IMOGEN_UNIT, FILE=FILE_CLIM,                                         &
                      STATUS='old', POSITION='rewind', ACTION='read')

    READ(IMOGEN_UNIT,*) LONGMIN_CLIM,LATMIN_CLIM,LONGMAX_CLIM, LATMAX_CLIM

    DO I = 1,n_imogen_land
      IF (SGINDINV(I) > 0) THEN
        L = SGINDINV(I)
        READ(IMOGEN_UNIT,*) LONG(L), LAT(L), T_CLIM(L,J), RH15M_CLIM(L,J),    &
                            UWIND_CLIM(L,J), VWIND_CLIM(L,J), LW_CLIM(L,J),   &
                            SW_CLIM(L,J), DTEMP_CLIM(L,J),                    &
                            RAINFALL_CLIM(L,J), SNOWFALL_CLIM(L,J),           &
                            PSTAR_HA_CLIM(L,J), F_WET_CLIM(L,J)
      ELSE
        READ(IMOGEN_UNIT,*)
      ENDIF
    ENDDO

    CLOSE(IMOGEN_UNIT)
  ENDDO

!-----------------------------------------------------------------------------
! Set up the initial conditions
!-----------------------------------------------------------------------------
  CO2_PPMV = CO2_INIT_PPMV

  IF ( INCLUDE_CO2 ) CO2_CHANGE_PPMV = 0.0

  IF ( ANLG ) THEN
    DTEMP_O(:) = 0.0

    IF ( INCLUDE_CO2 .AND. OCEAN_FEED .AND. C_EMISSIONS ) FA_OCEAN(:)=0.0
  ENDIF

! Initiate seeding values for subdaily rainfall
  SEED_RAIN(1) = 9465
  SEED_RAIN(2) = 1484
  SEED_RAIN(3) = 3358
  SEED_RAIN(4) = 8350
! Initiate seeding values for the weather generator
  IF ( WGEN ) THEN
    SEED_WG(1) = 5810
    SEED_WG(2) = 5575
    SEED_WG(3) = 5817
    SEED_WG(4) = 9119
  ENDIF

  cv(:) = 0.0

! Override the IMOGEN variables from given dump file if we have been asked to
  IF ( initialise_from_dump )                                                 &
    CALL read_dump(dump_file, (/ 'co2_ppmv       ', 'co2_change_ppmv',        &
                                 'dtemp_o        ', 'fa_ocean       ',        &
                                 'seed_rain      ', 'cv             ' /))

  RETURN

END SUBROUTINE init_imogen
#endif
