!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: I/O module for GSWP3 for CABLE offline version
!
! Contact: c.carouge@unsw.edu.au
!
! History: Developed by Claire Carouge
!
! ==============================================================================
!
! MODULEs used:
!
!==============================================================================

MODULE CABLE_GSWP3_MOD

  use netcdf
  USE CABLE_COMMON_MODULE, ONLY: HANDLE_ERR
  USE CABLE_IO_VARS_MODULE, ONLY: logn
  USE CABLE_METFORCING_MOD

  IMPLICIT NONE


  ! Some local parameters

  ! These will allow to access variable indexes by name instead of number.
  INTEGER, PRIVATE, PARAMETER :: &
       rain     =  1, &
       snow     =  2, &
       lwdn     =  3, &
       swdn     =  4, &
       pres     =  5, &
       rhum     =  6, &
       tair     =  7, &
       wind     =  8, &
       time     =  9

  ! Prefix name for files. Must be in same order as the variable indexes above.
  CHARACTER(len=6), DIMENSION(9), PARAMETER, PRIVATE :: &
       PREF = (/ "Rainf ", "Snowf ", "LWdown", "SWdown", "PSurf ", "Qair  ",   &
                 "Tair  ", "Wind  ", "time  " /)

CONTAINS

  SUBROUTINE READ_SETTINGS( METFORC )

    USE CABLE_COMMON_MODULE, ONLY: GET_UNIT
    USE cable_IO_vars_module, ONLY: logn

    TYPE (FORC_MET_TYPE) :: METFORC

    INTEGER              :: iu
    CHARACTER(len=200)   :: BasePath, LandMaskFile

    NAMELIST /GSWP3/ BasePath, LandMaskFile


    CALL GET_UNIT(iu)
    OPEN (iu,FILE="cable.nml",STATUS='OLD',ACTION='READ')
    READ (iu,NML=GSWP3)
    CLOSE(iu)

    METFORC%BasePath     = BasePath
    METFORC%LandMaskFile = LandMaskFile

    ! Print settings

    WRITE(*   ,*)"========================================= GSWP3 ============"
    WRITE(*   ,*)"GSWP3 settings chosen:"
    WRITE(*   ,*)" BasePath: ",TRIM(METFORC%BasePath)
    WRITE(*   ,*)" LandMask: ",TRIM(METFORC%LandMaskFile)
    WRITE(logn,*)"========================================= GSWP3 ============"
    WRITE(logn,*)"GSWP3 settings chosen:"
    WRITE(logn,*)" BasePath: ",TRIM(METFORC%BasePath)
    WRITE(logn,*)" LandMask: ",TRIM(METFORC%LandMaskFile)

  END SUBROUTINE READ_SETTINGS

  SUBROUTINE GSWP3_INIT(METFORC)
    ! Initialise met forcings. Reads time invariate. Populate names, etc in
    ! METFORC structure.
    !
    ! Arguments:
    ! METFORC: structure to keep info on met. forcings

    USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, LEAP_DAY,                      &
                                   HANDLE_ERR, GET_UNIT

    USE cable_IO_vars_module, ONLY: latitude, longitude, nmetpatches,          &
                                    mask, metGrid, sdoy, smoy, syear, shod,    &
                                    xdimsize, ydimsize, lat_all, lon_all,      &
                                    land_x, land_y

    USE cable_def_types_mod,  ONLY: mland

    IMPLICIT NONE

    TYPE (FORC_MET_TYPE), INTENT(INOUT) :: METFORC ! Type defined in cable_metforcing_mod

    TYPE(file_var_names)               :: f_varnames  ! Type defined in cable_metforcing_mod
    INTEGER                            :: FID
    INTEGER                            :: cnt, x, y ! loop indexes
    INTEGER, DIMENSION(2)              :: shape_all


    ! SET PARAMETERS, NAMES AND INPUTS
    !---------------------------------------------------------------------------
    ! Set number of variables and allocate array
    METFORC%NMET = 9
    IF (.NOT.ALLOCATED(METFORC%MET)) ALLOCATE(METFORC%MET( METFORC%NMET ))

    ! Set Leap-years according to dataset
    METFORC%LeapYears = .TRUE.

    ! Set variable names in files  !!!! Does not work with 1 structure for all !!!!
    METFORC%MET(rain)%VAR_NAME = "Rainf"
    METFORC%MET(snow)%VAR_NAME = "Snowf"
    METFORC%MET(lwdn)%VAR_NAME = "LWdown"
    METFORC%MET(swdn)%VAR_NAME = "SWdown"
    METFORC%MET(pres)%VAR_NAME = "PSurf"
    METFORC%MET(rhum)%VAR_NAME = "Qair"
    METFORC%MET(tair)%VAR_NAME = "Tair"
    METFORC%MET(wind)%VAR_NAME = "Wind"

    ! Set names of dimensions
    f_varnames%lat_dim_name = "lat"
    f_varnames%lat_var_name = "lat"
    f_varnames%lon_dim_name = "lon"
    f_varnames%lon_var_name = "lon"
    f_varnames%mask_var_name = "mask"

    ! Read GSWP3 settings
    CALL READ_SETTINGS(METFORC)

    WRITE(*   ,*)"========================================= GSWP3 ============"
    WRITE(logn,*)"========================================= GSWP3 ============"

    ! GET DIMENSIONS
    !---------------------------------------------------------------------------
    ! Now read landmask file, get LAT, LON etc. from there
    FID = OPEN_NETCDFFILE( METFORC%LandMaskFile )

    ! Read dimensions. Store dimensions in global arrays lat_all and lon_all.
    CALL READ_DIMENSIONS( f_varnames, FID, METFORC%LandMaskFile )

    ! Save # of longitudes and latitudes for allocation of arrays
    shape_all = SHAPE( lat_all )
    METFORC%xdimsize = shape_all(1)
    METFORC%ydimsize = shape_all(2)

    ! GET MASK
    !---------------------------------------------------------------------------
    CALL READ_LANDMASK( METFORC, f_varnames , FID, METFORC%LandMaskFile )

    ! Save # of land points
    METFORC%mland = COUNT(METFORC%landmask)

    ! Allocate variable arrrays. Time array allocated later because
    ! no time in landmask file.
    DO x = 1, METFORC%NMET-1
       ALLOCATE( METFORC%MET(x)%VAL(METFORC%mland) )
    END DO


    ! Set global CABLE variables.
    metGrid     = "mask"
    mland       = METFORC%mland
    nmetpatches = 1

    ! MAP TO LAND ONLY ARRAYS
    !---------------------------------------------------------------------------
    CALL MAP2LANDONLY(METFORC%landmask)


!    ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
!    shod        = 0.
!    sdoy        = 1
!    smoy        = 1
!    syear       = METFORC%CYEAR

!    IF ( ERR ) THEN
!       WRITE(logn,*)"Invalid settings in GSWP3_INIT"
!       STOP "Invalid settings in GSWP3_INIT"
!    ENDIF
    ! Close mask file
    CALL CLOSE_NETCDFFILE( FID, 'Closing file: '//TRIM(METFORC%LandMaskFile) )

  END SUBROUTINE GSWP3_INIT

  SUBROUTINE GSWP3_GET_FILENAME( METFORC, syear, par )
    ! Construct complete filenames for each variable. Assumes 1 var / file.

    TYPE(FORC_MET_TYPE), INTENT(INOUT) :: METFORC
    INTEGER, INTENT(IN)                :: par ! parameter ID
    INTEGER, INTENT(IN)                :: syear

    CHARACTER(len=4)                   :: year

    write(year,*) syear

    METFORC%MET(par)%MetFile = TRIM(METFORC%BasePath)
    METFORC%MET(par)%MetFile = TRIM(METFORC%MET(par)%MetFile)//PREF(par)//"/"
    METFORC%MET(par)%MetFile = TRIM(METFORC%MET(par)%MetFile)//"GSWP3.BC."//PREF(par)//".3hrMap."//year//".nc"

  END SUBROUTINE GSWP3_GET_FILENAME

  SUBROUTINE OPEN_GSWP3_MET( METFORC, syear )
    ! Open met files.

    TYPE(FORC_MET_TYPE), INTENT(INOUT) :: METFORC
    INTEGER,             INTENT(IN   ) :: syear

    INTEGER :: par

    DO par = 1, METFORC%NMET-1

       CALL GSWP3_GET_FILENAME( METFORC, syear, par )

       ! OPEN MET FILES and access variables' IDs
       WRITE(*   ,*) 'Opening met data file: ', METFORC%MET(par)%MetFile
       WRITE(logn,*) 'Opening met data file: ', METFORC%MET(par)%MetFile

       METFORC%MET(par)%F_ID = OPEN_NETCDFFILE( METFORC%MET(par)%MetFile )

       METFORC%MET(par)%V_ID = GET_VARID( METFORC%MET(par) )

    ENDDO

    ! Special case for time. Read from first file
    METFORC%MET(time)%MetFile  = METFORC%MET(1)%MetFile
    METFORC%MET(time)%VAR_NAME = PREF(time)
    METFORC%MET(time)%F_ID     = METFORC%MET(1)%F_ID

    METFORC%MET(time)%V_ID     = GET_VARID( METFORC%MET(time) )

  END SUBROUTINE OPEN_GSWP3_MET

  SUBROUTINE GET_TIME( METFORC )

     TYPE(FORC_MET_TYPE), INTENT(INOUT) :: METFORC
     INTEGER                            :: ntime

     ! Number of time steps in met file
     ntime = GET_TIME_DIMENSION( METFORC%MET(time) )

     ! Allocate time array
     IF (.NOT.ALLOCATED(METFORC%MET(time)%VAL)) ALLOCATE(METFORC%MET(time)%VAL(ntime))

     ! Read in time
     CALL GET_TIME_VAR( METFORC, time )

     ! Read in time coordinate attribute to check if local or GMT time.
     METFORC%time_coord = GET_TIME_COORD(METFORC%MET(time))

  END SUBROUTINE GET_TIME

  SUBROUTINE ASSIGN2MET(met, METVAR, NMET, rad)
    USE cable_IO_vars_module,    ONLY: landpt
    USE cable_def_types_mod,     ONLY: MET_TYPE, mland, RADIATION_TYPE
    USE cable_radiation_module,  ONLY: sinbet
    USE cable_checks_module,     ONLY: rh_sh

    ! strangely, met% is not save over Wait all in MPI...!
    ! Arguments:
    INTEGER,                             INTENT(IN   ) :: NMET
    TYPE(VAR_MET_TYPE), DIMENSION(NMET), INTENT(IN   ) :: METVAR
    TYPE(RADIATION_TYPE),                INTENT(IN   ) :: rad
    TYPE(MET_TYPE),                      INTENT(INOUT) :: met

    ! Local variables
    INTEGER :: i, is, ie   ! Dummy variables

    DO i = 1, mland
      is = landpt(i)%cstart
      ie = landpt(i)%cend
      met%precip    (is:ie)   = METVAR(rain)%VAL(i) * METVAR(rain)%convert &
                              + METVAR(snow)%VAL(i) * METVAR(snow)%convert
      met%precip_sn (is:ie)   = METVAR(snow)%VAL(i) * METVAR(snow)%convert
      met%fld       (is:ie)   = METVAR(lwdn)%VAL(i) * METVAR(lwdn)%convert
      met%fsd       (is:ie,1) = METVAR(swdn)%VAL(i) * METVAR(swdn)%convert &
                              * 0.5
      met%fsd       (is:ie,2) = METVAR(swdn)%VAL(i) * METVAR(swdn)%convert &
                              * 0.5
      met%tk        (is:ie)   = METVAR(tair)%VAL(i) + METVAR(tair)%convert
      met%ua        (is:ie)   = METVAR(wind)%VAL(i) * METVAR(wind)%convert
      met%pmb       (is:ie)   = METVAR(pres)%VAL(i) * METVAR(pres)%convert
     ! Set cosine of zenith angle (provided by GCM when online):
      met%coszen    (is:ie)   = sinbet(met%doy(is:ie), rad%latitude(is:ie),    &
                                       met%hod(is:ie))

      ! compute qv
      IF ( METVAR(rhum)%convert > 0. ) THEN ! We have specific humidity
           met%qv   (is:ie)   = METVAR(rhum)%VAL(i) ! WATCH "rhum'" is actually specific humidity (contrast with other data sets)
        ELSE ! Convert relative humidity to specific humidity
           CALL rh_sh ( METVAR(rhum)%VAL(i), met%tk(is), met%pmb(is), met%qv(is) )
           met%qv   (is:ie)   = met%qv(is)
        ENDIF
     END DO

     ! initialise within canopy air temp
     met%tvair     = met%tk
     met%tvrad     = met%tk

  END SUBROUTINE ASSIGN2MET

  SUBROUTINE GSWP3_GET_MET(dels, METFORC, met, rad, ktau, kend, tfrz, ncciy )
    !==============================================================================
    !
    ! Name: GSWP3_GET_MET
    !
    ! Purpose: Reading met. forcing from GSWP3 dataset.
    !
    !
    !
    !
    !
    !
    !
    !
    !
    ! CALLed from: cable_<offline>_driver
    !
    ! CALLs:
    !
    ! Input file: [GSWP3_landmask].nc
    !             [GSWP3_Snowf].nc
    !             [GSWP3_LWdown].nc
    !             [GSWP3_SWdown].nc
    !             [GSWP3_PSurf].nc
    !             [GSWP3_Qair].nc
    !             [GSWP3_Tair].nc
    !             [GSWP3_wind].nc
    !             [GSWP3_Rainf].nc
    !
    !==============================================================================

    USE cable_def_types_mod,    ONLY: MET_TYPE, radiation_type
    USE cable_IO_vars_module,   ONLY: timevar, timeunits, LANDPT, latitude,    &
                                      syear, smoy, sdoy, shod
    USE cable_common_module,    ONLY: YMDHMS2DOYSOD
    USE cable_checks_module,    ONLY: rh_sh
    USE cable_radiation_module, ONLY: sinbet

    IMPLICIT NONE

    INTEGER, INTENT(IN   ) :: ktau
    INTEGER, INTENT(IN   ) :: ncciy
    INTEGER, INTENT(OUT  ) :: kend
    REAL,    INTENT(IN   ) :: tfrz

    TYPE(FORC_MET_TYPE),  INTENT(INOUT) :: METFORC
    TYPE(MET_TYPE),       INTENT(INOUT) :: met
    TYPE(RADIATION_TYPE), INTENT(IN)    :: rad
    REAL,                 INTENT(OUT)   :: dels  ! Time-step for CABLE

    INTEGER   :: i, dY, dM, dD, is, ie, par
    INTEGER   :: sdoytmp, sod
    REAL      :: dt, etime
    CHARACTER :: LMFILE*200

    LOGICAL,        SAVE :: CALL1 = .TRUE.


    ! time-step relevant settings

    METFORC%ktau  = ktau

    IF ( CALL1 ) THEN
      ! Open met files
      CALL OPEN_GSWP3_MET( METFORC, syear )

      ! Time variables.
      ! Get time variable from met file
      CALL GET_TIME( METFORC )

      METFORC%MET(time)%VAL = (METFORC%MET(time)%VAL - (1))*3600.0 + 1.5*3600.0  !convert hours to seconds

      ! Store time in CABLE global variable
      kend = SIZE( METFORC%MET(time)%VAL )
      ALLOCATE( timevar( kend ))

      ! For GSWP3, need to convert to seconds and from start of current year.
      timevar = ( METFORC%MET(time)%VAL - METFORC%MET(time)%VAL(1) )           &   ! remove offset
                * 3600.0 + 1.5 * 3600.0                                            ! time-unit conversion
      !Hack the GSWP3 time units to make from start of year
      write(*,*) 'writing timeunits'
      write (timeunits, "('seconds since ',I4.4,'-01-01 00:00:00')") ncciy
      write(*,*) 'wrote time units'

      ! Set start time variables
      ! Assume run starts on jan 1st every year
      syear=ncciy
      smoy=1
      sdoytmp=1
      shod=0

      ! Convert to Day of year
      CALL YMDHMS2DOYSOD(syear, smoy, sdoytmp, INT(shod), 0, 0, sdoy, sod)

      ! Set end time variables
!      sim_days = INT((timevar(kend)-timevar(1) + sod)/(3600.0 * 24.0))
!      CALL ENDTIME(sim_days)

      ! Calculate time-step as the met. forcing time-step
      dels = REAL(timevar(2) - timevar(1))

      ! Unit conversions. Check units of all variables to see if they need
      ! some conversion.
      METFORC%MET(rain)%convert = CONVERT_RAIN_UNIT( METFORC%MET(rain), dels )
      METFORC%MET(snow)%convert = CONVERT_RAIN_UNIT( METFORC%MET(snow), dels )
      METFORC%MET(pres)%convert = CONVERT_PRES_UNIT( METFORC%MET(pres) )
      METFORC%MET(rhum)%convert = CONVERT_QAIR_UNIT( METFORC%MET(rhum) )
      METFORC%MET(tair)%convert = CONVERT_TAIR_UNIT( METFORC%MET(tair), tfrz )

      ! Other variables don't need conversion
      METFORC%MET(lwdn)%convert = 1.0
      METFORC%MET(swdn)%convert = 1.0
      METFORC%MET(wind)%convert = 1.0

    END IF

   ! Now get the Met data for this time step
   ! Set time.
   CALL SET_TIMEVARS(met, ktau, METFORC%time_coord, dels)

   ! Read in variables
   DO par = 1, METFORC%NMET-1  ! Don't read in last var. because it is time.
     CALL GET_VAR( METFORC, par, ktau )
   END DO

   ! assign to cable variables
   CALL ASSIGN2MET( met, METFORC%MET, METFORC%NMET, rad )

    ! Finally closing files when done
    IF (ktau .EQ. kend ) THEN
       DO par=1, METFORC%NMET-1
          CALL CLOSE_NETCDFFILE(METFORC%MET(par)%F_ID, 'Closing file: '//TRIM(METFORC%MET(par)%MetFile))
       END DO
    END IF

    ! CALL1 is over...
    CALL1 = .FALSE.

  END SUBROUTINE GSWP3_GET_MET

END MODULE CABLE_GSWP3_MOD
