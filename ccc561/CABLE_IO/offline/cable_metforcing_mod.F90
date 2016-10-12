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
! Purpose: Generic subroutines for met. forcing reading in CABLE offline version
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

MODULE CABLE_METFORCING_MOD
  USE netcdf
  USE cable_common_module, only : HANDLE_ERR, GET_UNIT

  IMPLICIT NONE

  !ccc Define met field variable type.
  TYPE VAR_MET_TYPE
     REAL, DIMENSION(:), ALLOCATABLE :: VAL           ! Values
     CHARACTER(len=12)               :: VAR_NAME      ! Variable name
     CHARACTER(len=200)              :: MetFile       ! Complete path to file
     INTEGER                         :: F_ID, V_ID    ! File and variable IDs
     REAL                            :: convert       ! Conversion value from file unit to CABLE unit
  END TYPE VAR_MET_TYPE

  TYPE FORC_MET_TYPE
     INTEGER  :: mland, NMET, xdimsize, ydimsize, tdimsize
     INTEGER  :: CYEAR, MetStart, MetEnd, CTSTEP, DT, ktau
     REAL,   DIMENSION(:)  ,ALLOCATABLE :: AVG_LWDN, CO2VALS
     LOGICAL  :: DirectRead, LeapYears
     LOGICAL,DIMENSION(:,:),ALLOCATABLE :: LandMask
     CHARACTER(len=15) :: Run,Forcing,RCP, CO2, NDEP,RCPdir
     CHARACTER(len=200):: BasePath, MetPath, LandMaskFile
     CHARACTER(len=3)  :: time_coord ! To know if local or GMT time.
     TYPE(VAR_MET_TYPE), DIMENSION(:), ALLOCATABLE :: MET
  END TYPE FORC_MET_TYPE

  TYPE file_var_names
     CHARACTER(len=9) :: lat_dim_name
     CHARACTER(len=9) :: lat_var_name
     CHARACTER(len=9) :: lon_dim_name
     CHARACTER(len=9) :: lon_var_name
     CHARACTER(len=8) :: mask_var_name
  END TYPE file_var_names


CONTAINS
   FUNCTION OPEN_NETCDFFILE( Filename, message ) RESULT( FID )
! Opens a netcdf file and returns its file ID.
!
! Arguments:
! Filename: full path of the file to open
! message: Optional error message to write if unsuccessful to open file
!
! Result:
! FID: file ID number, integer
    USE cable_IO_vars_module, only : logn

    ! Arguments
    CHARACTER(len=200), INTENT(IN)           :: FILENAME
    CHARACTER(len=256), INTENT(IN), OPTIONAL :: message

    ! Local variables
    INTEGER              :: STATUS
    INTEGER              :: FID
    CHARACTER(len=256)   :: mess

    mess = 'Opening file: '
    IF ( PRESENT(message) ) mess = message


    WRITE(*   ,*) mess,TRIM(FILENAME)
    WRITE(logn,*) mess,TRIM(FILENAME)

    STATUS = NF90_OPEN(TRIM(FILENAME), NF90_NOWRITE, FID)
    CALL HANDLE_ERR(STATUS, TRIM(mess)//TRIM(FILENAME))

  END FUNCTION OPEN_NETCDFFILE

!**************

  SUBROUTINE READ_DIMENSIONS(f_varnames, FID, FileName)
    ! Read latitudes and longitudes from FID file. FID is file ID from
    ! OPEN_NETCDFFILE() call. Save values in lat_all and lon_all.
    !
    ! Arguments:
    ! f_varnames: the dimension and variable names in the land mask file
    ! FID: file ID
    ! FileName: optional filename for error messages.

    USE cable_IO_vars_module, ONLY: lat_all, lon_all

    ! Arguments
    INTEGER, INTENT(IN)                      :: FID
    TYPE(file_var_names), INTENT(IN)         :: f_varnames
    CHARACTER(len=200), INTENT(IN), OPTIONAL :: FileName

    ! Local variables
    INTEGER                         :: STATUS, latID, lonID, ydimsize, xdimsize
    INTEGER                         :: x, y
    REAL,DIMENSION(:)  ,ALLOCATABLE :: lats, lons
    CHARACTER(len=200)              :: mess
    CHARACTER(len=9)                :: lat_dim, lon_dim, lat_var, lon_var

    mess=""
    IF (PRESENT(FileName)) mess = FileName

    ! Store the dimension_names in local variables to shorten var. names.
    lat_dim = f_varnames%lat_dim_name
    lat_var = f_varnames%lat_var_name
    lon_dim = f_varnames%lon_dim_name
    lon_var = f_varnames%lon_var_name

    ! lat dimension
    STATUS = NF90_INQ_DIMID(FID,TRIM(lat_dim),latID)
    STATUS = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
    CALL HANDLE_ERR(STATUS, "Inquiring "//TRIM(lat_dim)//" "//TRIM(mess))

    ALLOCATE( lats( ydimsize ) )
    STATUS = NF90_INQ_VARID(FID,TRIM(lat_var),latID)
    CALL HANDLE_ERR(STATUS, "Inquiring "//TRIM(lat_var)//" "//TRIM(mess))
    STATUS = NF90_GET_VAR(FID,latID,lats)
    CALL HANDLE_ERR(STATUS, "Reading "//TRIM(lat_var)//" "//TRIM(mess))

    ! lon dimension
    STATUS = NF90_INQ_DIMID(FID,TRIM(lon_dim),lonID)
    STATUS = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
    CALL HANDLE_ERR(STATUS, "Inquiring "//TRIM(lon_dim)//" "//TRIM(mess))

    ALLOCATE( lons( xdimsize ) )
    STATUS = NF90_INQ_VARID(FID,TRIM(lon_var),lonID)
    CALL HANDLE_ERR(STATUS, "Inquiring "//TRIM(lon_var)//" "//TRIM(mess))
    STATUS = NF90_GET_VAR(FID,lonID,lons)
    CALL HANDLE_ERR(STATUS, "Reading "//TRIM(lon_var)//" "//TRIM(mess))

    ! Save latitudes in lat_all global array
    ALLOCATE( lat_all(xdimsize, ydimsize) )
    ALLOCATE( lon_all(xdimsize, ydimsize) )
    DO x = 1, xdimsize
       lat_all(x,:) = lats
    END DO
    DO y = 1, ydimsize
       lon_all(:,y) = lons
    END DO

    ! Deallocate local variables
    DEALLOCATE(lats, lons)

  END SUBROUTINE READ_DIMENSIONS

!**************

  SUBROUTINE READ_LANDMASK( METFORC, f_varnames , FID, FileName )
    ! Read in the land/sea mask and store values in METFORC%landmask and
    ! the global variable mask.
    ! Assumes that mask on land points has values > 0. Water points have
    ! values < 0.
    !
    ! Arguments:
    ! METFORC: the met. forcing structure
    ! f_varnames: the dimension and variable names in the land mask file
    ! FID: file ID
    ! FileName: optional filename for error messages.

    USE cable_IO_vars_module, ONLY : mask

    ! Arguments
    TYPE(forc_met_type), INTENT(INOUT)       :: METFORC
    INTEGER, INTENT(IN)                      :: FID
    TYPE(file_var_names), INTENT(IN)         :: f_varnames
    CHARACTER(len=200), INTENT(IN), OPTIONAL :: FileName

    ! Local variables
    INTEGER                                  :: STATUS, landID
    REAL, DIMENSION(:,:), ALLOCATABLE        :: landmask
    CHARACTER(len=200)                       :: mess

    ! Optional arguments
    mess=" "
    IF (PRESENT( FileName )) mess = FileName

    ! Allocate arrays
    ALLOCATE( METFORC%landmask( METFORC%xdimsize, METFORC%ydimsize) )
    ALLOCATE( mask( METFORC%xdimsize, METFORC%ydimsize) )
    ALLOCATE( landmask( METFORC%xdimsize, METFORC%ydimsize) )

    ! get mask
    STATUS = NF90_INQ_VARID(FID,f_varnames%mask_var_name,landID)
    CALL HANDLE_ERR(STATUS, "Inquiring "//TRIM(f_varnames%mask_var_name)//" "//TRIM(mess))
    STATUS = NF90_GET_VAR(FID,landID,landmask)
    CALL HANDLE_ERR(STATUS, "Reading "//TRIM(f_varnames%mask_var_name)//" "//TRIM(mess))

    ! Populate global arrays
    WHERE ( landmask .GT. 0 )
       METFORC%landmask = .TRUE.
       mask           = 1
    ELSEWHERE
       METFORC%landmask = .FALSE.
       mask           = 0
    END WHERE

    ! Deallocate local variables
    DEALLOCATE( landmask )

  END SUBROUTINE READ_LANDMASK

!**************
  subroutine MAP2LANDONLY(landmask)

    USE cable_IO_vars_module, ONLY: latitude, longitude, land_x, land_y,       &
                                    lon_all, lat_all, xdimsize, ydimsize

    USE cable_def_types_mod, ONLY: mland

    ! Arguments
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: landmask

    ! Local variables
    INTEGER :: cnt, y, x

    ! Allocate global variables
    ALLOCATE( latitude( mland ), longitude( mland ) )
    ALLOCATE( land_y  ( mland ), land_x   ( mland ) )

    ! Map all to Landgrid arrays
    cnt = 1
    DO y = 1, ydimsize
       DO x = 1, xdimsize
          IF ( .NOT. landmask(x,y) ) CYCLE
          WRITE(6,FMT='(A15,I6,2(1X,F8.2),2(1x,I3))')"i, lo,la, x,y",cnt,      &
                     lon_all(x,y),lat_all(x,y),x, y

          land_x   (cnt) = x
          land_y   (cnt) = y
          longitude(cnt) = lon_all(x,y)
          latitude (cnt) = lat_all(x,y)
          cnt = cnt + 1
       END DO
    END DO

  end subroutine
!**************

   SUBROUTINE CLOSE_NETCDFFILE( FID, message )
! Closes a netcdf file
!
! Arguments:
! FID: file ID to close
! message: Optional error message to write if unsuccessful to close file
!

    USE cable_IO_vars_module, only : logn

    ! Arguments
    INTEGER, INTENT(IN)                      :: FID
    CHARACTER(len=256), INTENT(IN), OPTIONAL :: message

    ! Local variables
    INTEGER              :: STATUS
    CHARACTER(len=256)   :: mess

    mess = 'Closing file: '
    IF ( PRESENT(message) ) mess = message

    WRITE(*   ,*) mess
    WRITE(logn,*) mess

    STATUS = NF90_CLOSE(FID)
    CALL HANDLE_ERR(STATUS, mess)

  END SUBROUTINE CLOSE_NETCDFFILE

!**************

  FUNCTION GET_VARID( MET_STRUCT ) RESULT( VID )
    ! Get the variables' IDs from met. forcing file.
    USE cable_common_module, only : HANDLE_ERR

    TYPE(VAR_MET_TYPE), INTENT(IN) :: MET_STRUCT
    INTEGER                        :: VID

    INTEGER :: STATUS

    STATUS = NF90_INQ_VARID(MET_STRUCT%F_ID,TRIM(MET_STRUCT%VAR_NAME), VID)
    CALL HANDLE_ERR(STATUS, "Inquiring forcing var "//TRIM(MET_STRUCT%VAR_NAME)// &
         " in "//MET_STRUCT%MetFile )

  END FUNCTION GET_VARID

!**************

  FUNCTION GET_TIME_DIMENSION( METVAR ) RESULT( NTIME )

    TYPE(VAR_MET_TYPE), INTENT(IN) :: METVAR
    INTEGER                        :: NTIME
    INTEGER                        :: STATUS
    CHARACTER(len=256)             :: mess

    ! Error message
    mess = "Inquiring "//TRIM(METVAR%VAR_NAME)//" in file "//TRIM(METVAR%MetFile)

    ! Get dimension
    STATUS = NF90_INQUIRE_DIMENSION(METVAR%F_ID,METVAR%V_ID,len=NTIME)

    ! Handle error gracefully
    CALL HANDLE_ERR(STATUS, TRIM(mess))

  END FUNCTION GET_TIME_DIMENSION

!**************

  SUBROUTINE GET_VAR( METFORC, par, tstart )

    USE cable_IO_vars_module, ONLY : land_x, land_y

    ! Read the forcing values. Note we can't do a simple read as we only want to
    ! keep the land points! Also only read in 1 time step.
    TYPE(FORC_MET_TYPE), INTENT(INOUT) :: METFORC ! Info structure on met forcing
    INTEGER, INTENT(IN)                :: par     ! Variable to read in
    INTEGER, INTENT(IN)                :: tstart  ! Time step to read in

    INTEGER                            :: STATUS
    INTEGER                            :: xds, yds, k
    REAL, DIMENSION(:,:), ALLOCATABLE  :: tmparr
    REAL                               :: tmp

    xds = METFORC%xdimsize
    yds = METFORC%ydimsize

    IF ( METFORC%DirectRead ) THEN

      DO k = 1, METFORC%mland

        STATUS = NF90_GET_VAR( METFORC%MET(par)%F_ID, METFORC%MET(par)%V_ID, tmp, &
                               start=(/land_x(k),land_y(k),tstart/) )
        CALL HANDLE_ERR(STATUS, "Reading direct from "//METFORC%MET(par)%MetFile )
        METFORC%MET(par)%VAL(k) = tmp
      END DO

    ELSE
      ALLOCATE( tmparr( xds, yds ))

      STATUS = NF90_GET_VAR( METFORC%MET(par)%F_ID, METFORC%MET(par)%V_ID, tmparr, &
                             start=(/1,1,tstart/),count=(/xds,yds,1/) )
      CALL HANDLE_ERR(STATUS, "Reading from "//METFORC%MET(par)%MetFile )

      DO k = 1, METFORC%mland
        METFORC%MET(par)%VAL(k) = tmparr( land_x(k), land_y(k) )
      END DO

      DEALLOCATE( tmparr )
    ENDIF

  END SUBROUTINE GET_VAR

!**************

  SUBROUTINE GET_TIME_VAR( METFORC, par )
    TYPE(FORC_MET_TYPE), INTENT(INOUT) :: METFORC
    INTEGER, INTENT(IN)                :: par

    INTEGER                           :: STATUS

    STATUS = NF90_GET_VAR( METFORC%MET(par)%F_ID, METFORC%MET(par)%V_ID,       &
                           METFORC%MET(par)%VAL)
    CALL HANDLE_ERR(STATUS, "Reading from "//METFORC%MET(par)%MetFile )

  END SUBROUTINE GET_TIME_VAR

!**************

function GET_TIME_COORD(METVAR) RESULT ( time_coord )

  USE cable_def_types_mod, ONLY: mland
  ! Arguments:
  TYPE(VAR_MET_TYPE), INTENT(INOUT) :: METVAR

  ! Local variables
  INTEGER                           :: STATUS
  CHARACTER(len=3)                  :: time_coord

  ! Get coordinate field:
  STATUS = NF90_GET_ATT(METVAR%F_ID,METVAR%V_ID,'coordinate',time_coord)

  ! If error getting coordinate field (i.e. it doesn't exist):
  IF(STATUS /= NF90_NOERR) THEN
     ! Assume default time coordinate:
     IF(mland==1) THEN ! If single site, this is local time
        time_coord = 'LOC' ! 12am is 12am local time, at site/gridcell
     ELSE ! If multiple/global/regional, use GMT
        time_coord = 'GMT' ! 12am is GMT time, local time set by longitude
     END IF
  ELSE IF((STATUS==NF90_NOERR.AND.time_coord=='LOC'.AND.mland>1)) THEN
     ! Else if local time is selected for regional simulation, abort:
     CALL abort('"time" variable must be GMT for multiple site simulation!' &
          //' Check "coordinate" field in time variable.' &
          //' (SUBROUTINE GET_TIME_COORD)')
  ELSE IF(time_coord/='LOC'.AND.time_coord/='GMT') THEN
     CALL abort('Meaningless time coordinate in met data file!' &
          // ' (SUBROUTINE GET_TIME_COORD)')
  END IF

end function GET_TIME_COORD

!**************

  ! subroutine ENDTIME(sim_days)
  !   USE cable_IO_vars_module, ONLY: timevar, syear, smoy, sdoy, shod,          &
  !                                            eyear, emoy, edoy, ehod
  !
  !   !Arguments:
  !   INTEGER, INTENT(IN) :: sim_days ! Total number of days for the simulation.
  !
  !   !Local variable:
  !   INTEGER :: days_in_year
  !
  !   edoy  = sdoy + sim_days
  !   eyear = syear
  !   ehod  = shod
  !
  !   ! If running more than 1 year.
  !   DO WHILE (edoy > days_in_year)
  !     ! Check total number of days in the year for each year.
  !     IF ( IS_LEAPYEAR( eyear )) THEN
  !       days_in_year = 366
  !     ELSE
  !       days_in_year = 365
  !     END IF
  !
  !     edoy = edoy - days_in_year
  !     eyear = eyear + 1
  !   END DO
  !
  ! end subroutine ENDTIME

!**************

  subroutine SET_TIMEVARS(met, ktau, time_coord, dels)

    USE cable_IO_vars_module, ONLY: landpt, shod, syear, sdoy, smoy, longitude
    USE cable_common_module,  ONLY: DOYSOD2YMDHMS, IS_LEAPYEAR
    USE cable_def_types_mod,  ONLY: MET_TYPE, mland

    ! Arguments:
    INTEGER,          INTENT(IN   ) :: ktau  ! Current Time-step
    REAL,             INTENT(IN   ) :: dels  ! Time-step
    CHARACTER(len=3), INTENT(IN   ) :: time_coord ! Local or GMT time
    TYPE(MET_TYPE),   INTENT(INOUT) :: met   ! met forcing for CABLE

    ! Local variables:
    INTEGER :: i, fcell, dd ! Dummy variables
    INTEGER :: sod  ! seconds of day

    DO i=1,mland ! over all land points/grid cells
      ! First set timing variables:
      ! All timing details below are initially written to the first patch
      ! of each gridcell, then dumped to all patches for the gridcell.
      fcell = landpt(i)%cstart

      IF(ktau==1) THEN ! initialise...
        SELECT CASE(time_coord)
        CASE('LOC')! i.e. use local time by default
            ! hour-of-day = starting hod
            met%hod(fcell) = shod
            met%doy(fcell) = sdoy
            met%moy(fcell) = smoy
            met%year(fcell) = syear
        CASE('GMT')! use GMT
            ! hour-of-day = starting hod + offset from GMT time:
            met%hod(fcell) = shod + (longitude(i)/180.0)*12.0
            ! Note above that all met%* vars have dim mp,
            ! while longitude and latitude have dimension mland.
            met%doy(fcell) = sdoy
            met%moy(fcell) = smoy
            met%year(fcell) = syear
        CASE DEFAULT
            CALL abort('Unknown time coordinate! ' &
                //' (SUBROUTINE SET_TIMEVARS)')
        END SELECT
      ELSE
        ! increment hour-of-day by time step size:
        met%hod(fcell) = met%hod(fcell) + dels/3600.0
      END IF
!
      IF(met%hod(fcell)<0.0) THEN ! may be -ve since longitude
        ! has range [-180,180]
        ! Reduce day-of-year by one and ammend hour-of-day:
        met%doy(fcell) = met%doy(fcell) - 1
        met%hod(fcell) = met%hod(fcell) + 24.0
        sod = met%hod(fcell) * 3600.

        ! Update year.
        ! If doy is 0 then it's the last day of previous year
        IF ( met%doy(fcell) == 0. ) THEN
          met%year = met%year - 1

          IF ( IS_LEAPYEAR( syear-1 ) ) THEN
            met%doy(fcell) = 366
          ELSE
            met%doy(fcell) = 365
          END IF
        END IF


      ELSE IF(met%hod(fcell)>=24.0) THEN
        ! increment or GMT adj has shifted day
        ! Adjust day-of-year and hour-of-day:
        met%doy(fcell) = met%doy(fcell) + 1
        met%hod(fcell) = met%hod(fcell) - 24.0

        ! First day of next year case
        IF (( IS_LEAPYEAR(met%year(fcell)) .AND. met%doy(fcell) == 367) .OR. &
            ( .NOT. IS_LEAPYEAR(met%year(fcell)) .AND. met%doy(fcell) == 366 ))THEN
            met%doy(fcell) = 1
            met%year(fcell) = met%year(fcell) + 1
        END IF

      END IF ! if increment has pushed hod to a different day

      ! Update month
      CALL DOYSOD2YMDHMS( met%year(fcell), INT(met%doy(fcell)), sod,              &
                          met%moy(fcell), dd)

      ! Now copy these values to all veg/soil patches in the current grid cell:
      met%hod(fcell:landpt(i)%cend) = met%hod(fcell)
      met%doy(fcell:landpt(i)%cend) = met%doy(fcell)
      met%moy(fcell:landpt(i)%cend) = met%moy(fcell)
      met%year(fcell:landpt(i)%cend) = met%year(fcell)
    ENDDO

  end subroutine SET_TIMEVARS

!*********

  FUNCTION CONVERT_RAIN_UNIT( VAR_MET, dels ) RESULT( convert )
    TYPE(VAR_MET_TYPE), INTENT(IN) :: VAR_MET
    REAL,               INTENT(IN) :: dels  ! Time-step

    INTEGER                        :: STATUS
    CHARACTER(len=20)              :: unit_in
    REAL                           :: convert

    ! Get Rainf units:
    STATUS = NF90_GET_ATT(VAR_MET%F_ID,VAR_MET%V_ID,'units',unit_in)
    CALL HANDLE_ERR(STATUS, "Reading unit in "//TRIM( VAR_MET%MetFile ))

    IF(unit_in(1:8 ) == 'kg/m^2/s'   .OR.                                      &
       unit_in(1:10) == 'kgm^-2s^-1' .OR.                                      &
       unit_in(1:4 ) == 'mm/s'       .OR.                                      &
       unit_in(1:6 ) == 'mms^-1'     .OR.                                      &
       unit_in(1:7 ) == 'kg/m^2s'    .OR.                                      &
       unit_in(1:10) == 'kg m-2 s-1'     ) THEN
       ! Change from mm/s to mm/time step:
       convert = dels
    ELSE IF(                                                                   &
       unit_in(1:4 ) == 'mm/h'       .OR.                                      &
       unit_in(1:6 ) == 'mmh^-1'         ) THEN
       ! Change from mm/h to mm/time step:
       convert = dels/3600.0
    ELSE
       WRITE(*,*) unit_in
       CALL abort('Unknown units for Rainf'// &
         ' in '//TRIM(VAR_MET%MetFile)//' (FUNCTION convert_rain_unit)')
    END IF
  END FUNCTION CONVERT_RAIN_UNIT

!**************

  FUNCTION CONVERT_TAIR_UNIT( VAR_MET, tfrz ) RESULT( convert )
    TYPE(VAR_MET_TYPE), INTENT(IN) :: VAR_MET
    REAL,               INTENT(IN) :: tfrz

    INTEGER                        :: STATUS
    CHARACTER(len=20)              :: unit_in
    REAL                           :: convert

    ! Get Tair units:
    STATUS = NF90_GET_ATT(VAR_MET%F_ID,VAR_MET%V_ID,'units',unit_in)
    CALL HANDLE_ERR(STATUS, "Reading unit in "//TRIM( VAR_MET%MetFile ))

    IF(unit_in(1:1) == 'C' .OR.                                                &
       unit_in(1:1) == 'c'     ) THEN
       ! Change from Celsius to Kelvin step:
       convert = tfrz
    ELSE IF(                                                                   &
       unit_in(1:1) == 'K' .OR.                                                &
       unit_in(1:1) == 'k'     ) THEN
       ! Correct unit
       convert = 0.
    ELSE
       WRITE(*,*) unit_in
       CALL abort('Unknown units for Tair'// &
         ' in '//TRIM(VAR_MET%MetFile)//' (FUNCTION convert_tair_unit)')
    END IF
  END FUNCTION CONVERT_TAIR_UNIT

!**************

  FUNCTION CONVERT_QAIR_UNIT( VAR_MET ) RESULT( convert )
    TYPE(VAR_MET_TYPE), INTENT(IN) :: VAR_MET

    INTEGER                        :: STATUS
    CHARACTER(len=20)              :: unit_in
    REAL                           :: convert

    ! Get Qair units:
    STATUS = NF90_GET_ATT(VAR_MET%F_ID,VAR_MET%V_ID,'units',unit_in)
    CALL HANDLE_ERR(STATUS, "Reading unit in "//TRIM( VAR_MET%MetFile ))

    IF(unit_in(1:1) == '%' .OR.                                                &
       unit_in(1:1) == '-'     ) THEN
       ! Change from relative humidity to specific humidity:
       convert = -999.0
    ELSE IF(                                                                   &
       unit_in(1:3) == 'g/g'       .OR.                                        &
       unit_in(1:5) == 'kg/kg'     .OR.                                        &
       unit_in(1:3) == 'G/G'       .OR.                                        &
       unit_in(1:5) == 'KG/KG'     .OR.                                        &
       unit_in(1:7) == 'kg kg-1'       ) THEN
       ! Units are correct:
       convert = 1.0
    ELSE
       WRITE(*,*) unit_in
       CALL abort('Unknown units for Qair'// &
         ' in '//TRIM(VAR_MET%MetFile)//' (FUNCTION convert_qair_unit)')
    END IF
  END FUNCTION CONVERT_QAIR_UNIT

!**************

  FUNCTION CONVERT_PRES_UNIT( VAR_MET ) RESULT( convert )
    TYPE(VAR_MET_TYPE), INTENT(IN) :: VAR_MET

    INTEGER                        :: STATUS
    CHARACTER(len=20)              :: unit_in
    REAL                           :: convert

    ! Get pressure units:
    STATUS = NF90_GET_ATT(VAR_MET%F_ID,VAR_MET%V_ID,'units',unit_in)
    CALL HANDLE_ERR(STATUS, "Reading unit in "//TRIM( VAR_MET%MetFile ))

    IF(unit_in(1:2) == 'Pa' .OR.                                               &
       unit_in(1:2) == 'pa' .OR.                                               &
       unit_in(1:2) == 'PA'     ) THEN
       ! Change from pa to mbar:
       convert = 0.01
    ELSE IF(                                                                   &
       unit_in(1:2) == 'KP' .OR.                                               &
       unit_in(1:2) == 'kP' .OR.                                               &
       unit_in(1:2) == 'Kp' .OR.                                               &
       unit_in(1:2) == 'kp'     ) THEN
       ! Change from kPa to mbar:
       convert = 10.0
    ELSE IF(                                                                   &
       unit_in(1:2) == 'MB' .OR.                                               &
       unit_in(1:2) == 'mB' .OR.                                               &
       unit_in(1:2) == 'Mb' .OR.                                               &
       unit_in(1:2) == 'mb'     ) THEN
       ! Units are correct:
       convert = 1.0
    ELSE
       WRITE(*,*) unit_in
       CALL abort('Unknown units for pressure'// &
         ' in '//TRIM(VAR_MET%MetFile)//' (FUNCTION convert_pres_unit)')
    END IF
  END FUNCTION CONVERT_PRES_UNIT

END MODULE CABLE_METFORCING_MOD
