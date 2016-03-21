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
! Purpose: tracking climate variables for use in phenology and potential pft modules
!
! Called from: cable_driver
!
! History: Vanessa Haverd Jan 2015

! ==============================================================================
MODULE cable_climate_mod

 Use cable_def_types_mod, ONLY: met_type, climate_type, canopy_type, mp, &
      r_2, alloc_cbm_var
 USE TypeDef,              ONLY: i4b, dp
 USE cable_IO_vars_module, ONLY: patch
 USE CABLE_COMMON_MODULE, ONLY: CurYear, filename, cable_user, HANDLE_ERR

CONTAINS
! ==============================================================================


SUBROUTINE cable_climate(ktau,kstart,kend,ktauday,idoy,LOY,met,climate, canopy)


  IMPLICIT NONE

  INTEGER,      INTENT(IN) :: ktau ! integration step number
  INTEGER,      INTENT(IN) :: kstart ! starting value of ktau
  INTEGER,      INTENT(IN) :: kend ! total # timesteps in run

  INTEGER,      INTENT(IN)                  :: idoy ,LOY ! day of year (1-365) , Length oy
  INTEGER,      INTENT(IN)                  :: ktauday
  TYPE (met_type), INTENT(IN)       :: met  ! met input variables
  TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
  TYPE (canopy_type), INTENT(IN) :: canopy ! vegetation variables
  INTEGER :: d, y, k
  INTEGER, PARAMETER:: COLDEST_DAY_NHEMISPHERE = 355
  INTEGER, PARAMETER:: COLDEST_DAY_SHEMISPHERE = 172
  real,      dimension(mp)  :: mtemp_last, mmoist_last
  integer :: startyear
  integer :: MonthDays(12)
  integer::  DaysInMonth, nmonth, tmp
  logical :: IsLastDay ! last day of month?

  climate%doy = idoy
  ! accumulate daily temperature
  IF(MOD(ktau,ktauday)==1) THEN
     climate%dtemp = met%tk - 273.15
     climate%dmoist = canopy%fwsoil
  ELSE
     climate%dtemp = climate%dtemp + met%tk - 273.15
     climate%dmoist = climate%dmoist + canopy%fwsoil
  ENDIF

  IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
     climate%dtemp = climate%dtemp/FLOAT(ktauday)
     climate%dmoist = climate%dmoist/FLOAT(ktauday)
  ENDIF



  IF(MOD((ktau-kstart+1),ktauday)==0) THEN  ! end of day
     ! get month and check if end of month
     IsLastDay = .FALSE.
     MonthDays = (/31,28,31,30,31,30,31,31,30,31,30,31/)
     IF (LOY==366) MonthDays(2) = MonthDays(2) + 1
     nmonth = 1
     tmp = MonthDays(1)
     DO WHILE(tmp.LT.idoy)
        tmp = tmp +  MonthDays(nmonth+1)
        nmonth = nmonth + 1
     ENDDO

     if (idoy == sum(MonthDays(1:nmonth))) IsLastDay = .TRUE.

     ! On first day of year ...
     IF (idoy==1) THEN

        ! ... reset annual GDD5 counter
        climate%agdd5=0.0

     ENDIF

     WHERE ((patch%latitude>=0.0 .and. idoy==COLDEST_DAY_NHEMISPHERE).OR. &
          (patch%latitude<0.0 .and. idoy==COLDEST_DAY_SHEMISPHERE) )

        ! In midwinter, reset GDD counter for summergreen phenology
        climate%gdd5=0.0
     END WHERE

     ! Update GDD counters and chill day count

     climate%gdd5 = climate%gdd5 + max(0.0,climate%dtemp-5.0)
     climate%agdd5= climate%agdd5 + max(0.0,climate%dtemp-5.0)
     WHERE (climate%dtemp<5.0 .and. climate%chilldays<=365)
        climate%chilldays = climate%chilldays + 1
     ENDWHERE

     ! Save yesterday's mean temperature for the last month
     mtemp_last=climate%mtemp

     ! Update daily temperatures, and mean overall temperature, for last 31 days

     climate%mtemp=climate%dtemp
     climate%mmoist = climate%dmoist
     DO d=1,30
        climate%dtemp_31(:,d)=climate%dtemp_31(:,d+1)
        climate%mtemp= climate%mtemp + climate%dtemp_31(:,d)
        climate%dmoist_31(:,d)=climate%dmoist_31(:,d+1)
        climate%mmoist = climate%mmoist + climate%dmoist_31(:,d)
     ENDDO
     climate%dtemp_31(:,31)=climate%dtemp
     climate%dmoist_31(:,31)=climate%dmoist
     climate%mtemp = climate%mtemp/31.0
     climate%mmoist = climate%mmoist/31.0

     ! Reset GDD and chill day counter if mean monthly temperature falls below base
     ! temperature

     WHERE(mtemp_last>=5.0 .and. climate%mtemp<5.0)
        climate%gdd5=0.0
        climate%chilldays=0
     ENDWHERE

     ! On last day of month ...

     if (IsLastDay) THEN

        ! Update mean temperature for the last 12 months
        ! atemp_mean_new = atemp_mean_old * (11/12) + mtemp * (1/12)

        climate%atemp_mean=climate%atemp_mean*(11./12.)+climate%mtemp*(1./12.)

        ! Record minimum and maximum monthly temperatures

        if (nmonth==1) THEN
           climate%mtemp_min=climate%mtemp;
           climate%mtemp_max=climate%mtemp;

        else
           where (climate%mtemp<climate%mtemp_min) &
                climate%mtemp_min=climate%mtemp
           where (climate%mtemp>climate%mtemp_max) &
                climate%mtemp_max=climate%mtemp
        ENDIF  ! first month of year

        ! On 31 December update records of minimum monthly temperatures for the last
        ! 20 years and find minimum monthly temperature for the last 20 years

        if (nmonth==12) THEN
           climate%nyears = climate%nyears +1


           startyear=20-min(19,climate%nyears-1)
           climate%mtemp_min20=0.0
           climate%mtemp_max20=0.0

           if (startyear<20) then
              DO y=startyear,19
                 climate%mtemp_min_20(:,y)=climate%mtemp_min_20(:,y+1)
                 climate%mtemp_min20=climate%mtemp_min20+ &
                      climate%mtemp_min_20(:,y)
                 climate%mtemp_max_20(:,y)=climate%mtemp_max_20(:,y+1)
                 climate%mtemp_max20 =climate%mtemp_max20 + &
                      climate%mtemp_max_20(:,y)
              ENDDO

              climate%mtemp_min20=climate%mtemp_min20/real(21-startyear)
              climate%mtemp_max20=climate%mtemp_max20/real(21-startyear)
           else
              climate%mtemp_min20 = climate%mtemp_min
              climate%mtemp_max20 = climate%mtemp_max
           endif

           climate%mtemp_min_20(:,20)=climate%mtemp_min
           climate%mtemp_max_20(:,20)=climate%mtemp_max

        ENDIF  ! last month of year

!!$k=2
!!$write(559,92), climate%nyears, nmonth, idoy, climate%chilldays(k), &
!!$ climate%gdd5(k), climate%agdd5(k), &
!!$climate%mtemp(k), climate%mtemp_min(k), climate%mtemp_max(k), &
!!$climate%mtemp_min20(k), climate%mtemp_max20(k),climate%atemp_mean(k), &
!!$climate%mtemp_min_20(k,:), climate%mtemp_max_20(k,:), climate%mmoist(k)
!!$
!!$92    format(4(i6,',',2x),100(e18.6,',',2x))

     ENDIF     ! last day of month


  ENDIF ! end of day

END SUBROUTINE cable_climate

! ==============================================================================

SUBROUTINE climate_init ( climate,np )
IMPLICIT NONE

TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables
INTEGER, INTENT(IN) :: np
INTEGER :: d

CALL alloc_cbm_var(climate,np)

if (cable_user%climate_fromzero) then

   DO d=1,31
      !climate%dtemp_31(:,d)= climate%dtemp
     ! climate%dmoist_31(:,d)= climate%dmoist
      climate%dtemp_31(:,d)= 0
      climate%dmoist_31(:,d)= 0

   ENDDO
   climate%atemp_mean=0



   climate%nyears = 0
   climate%chilldays = 0
   climate%mtemp = 0
   climate%mmoist = 0
   climate%mtemp_min = 0
   climate%mtemp_max=0
   climate%mtemp_min20 =0
   climate%mtemp_max20=0
   climate%AGDD5=0
   climate%GDD5=0
   climate%GDD0=0
   climate%alpha_PT=0
   climate%aevap_PT=0
   climate%aevap=0
   climate%mtemp_min_20=0
   climate%mtemp_max_20=0


else
   CALL READ_CLIMATE_RESTART_NC (climate)

endif

END SUBROUTINE climate_init

! ==============================================================================

SUBROUTINE WRITE_CLIMATE_RESTART_NC ( climate )

  USE netcdf


  IMPLICIT NONE

  TYPE (climate_type), INTENT(IN)       :: climate  ! climate variables

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, i
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! 0 dim arrays
  CHARACTER(len=20),DIMENSION(2) :: A0
  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(13) :: A1
 ! 1 dim arrays (integer) (npt )
  CHARACTER(len=20),DIMENSION(1) :: AI1
  ! 2 dim arrays (npt,20)
  CHARACTER(len=20),DIMENSION(2) :: A2
  ! 2 dim arrays (npt,31)
  CHARACTER(len=20),DIMENSION(2) :: A3


  INTEGER*4 ::  VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1)), &
       VID2(SIZE(A2)), VID3(SIZE(A3))

  mp4=int(mp,fmp4)
  A0(1) = 'nyears'
  A0(2) = 'year'

  A1(1) = 'latitude'
  A1(2) = 'longitude'
  A1(3) = 'dtemp'
  A1(4) = 'mtemp'
  A1(5) = 'mtemp_min'
  A1(6) = 'mtemp_max'
  A1(7) = 'mtemp_min20'
  A1(8) = 'mtemp_max20'
  A1(9) = 'atemp_mean'
  A1(10) = 'AGDD5'
  A1(11) = 'GDD5'
  A1(12) = 'GDD0'
  A1(13) = 'mmoist'

  AI1(1) = 'chilldays'

  A2(1) = 'mtemp_min_20'
  A2(2) = 'mtemp_max_20'

  A3(1) = 'dtemp_31'
  A3(2) = 'dmoist_31'


  ! Get File-Name
  WRITE(CYEAR, FMT='(I4)') CurYear + 1
  fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
       '_climate_rst.nc'
  ! Create NetCDF file:
  STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  ! Put the file in define mode:
  STATUS = NF90_redef(FILE_ID)

  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )

  ! Define dimensions:
  ! Land (number of points)
  STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  ! number of years (stored for 20 y running means0
  STATUS = NF90_def_dim(FILE_ID, 'nyear' , 20 , nyear_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  ! number of days (stored for 31 day monthly means)
  STATUS = NF90_def_dim(FILE_ID, 'nday' , 31 , nday_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  DO i = 1, SIZE(A0)
     STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_INT ,VID0(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO


  DO i = 1, SIZE(A1)
     STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(AI1)
     STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A2)
     STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,nyear_ID/),VID2(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A3)
     STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,nday_ID/),VID3(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO



  ! End define mode:
  STATUS = NF90_enddef(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  ! PUT nyears and current year
  STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), climate%nyears )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), CurYear )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT LAT / LON
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), patch%latitude )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), patch%longitude )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT VARS
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), climate%dtemp )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), climate%mtemp )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), climate%mtemp_min )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), climate%mtemp_max )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), climate%mtemp_min20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), climate%mtemp_max20  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), climate%atemp_mean  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), climate%AGDD5  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), climate%GDD5  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), climate%GDD0 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), climate%mmoist )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), climate%chilldays )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), climate%mtemp_min_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), climate%mtemp_max_20 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), climate%dtemp_31 )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! Close NetCDF file:
  STATUS = NF90_close(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


END SUBROUTINE WRITE_CLIMATE_RESTART_NC
! ==============================================================================

SUBROUTINE READ_CLIMATE_RESTART_NC ( climate )

  USE netcdf


  IMPLICIT NONE

  TYPE (climate_type), INTENT(INOUT)       :: climate  ! climate variables

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, nyear_ID, nday_ID, dID, i, land_dim
  CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! 0 dim arrays
  CHARACTER(len=20),DIMENSION(2) :: A0
  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(13) :: A1
 ! 1 dim arrays (integer) (npt )
  CHARACTER(len=20),DIMENSION(1) :: AI1
  ! 2 dim arrays (npt,20)
  CHARACTER(len=20),DIMENSION(2) :: A2
  ! 2 dim arrays (npt,31)
  CHARACTER(len=20),DIMENSION(2) :: A3

  REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
  REAL(r_2)                         :: TMP2(mp,20),TMP3(mp,31)
  INTEGER*4 :: TMPI(mp), TMPI0
  LOGICAL            ::  EXISTFILE

  mp4=int(mp,fmp4)
  A0(1) = 'nyears'
  A0(2) = 'year'

  A1(1) = 'latitude'
  A1(2) = 'longitude'
  A1(3) = 'dtemp'
  A1(4) = 'mtemp'
  A1(5) = 'mtemp_min'
  A1(6) = 'mtemp_max'
  A1(7) = 'mtemp_min20'
  A1(8) = 'mtemp_max20'
  A1(9) = 'atemp_mean'
  A1(10) = 'AGDD5'
  A1(11) = 'GDD5'
  A1(12) = 'GDD0'
  A1(13) = 'mmoist'

  AI1(1) = 'chilldays'

  A2(1) = 'mtemp_min_20'
  A2(2) = 'mtemp_max_20'

  A3(1) = 'dtemp_31'
  A3(2) = 'dmoist_31'


  ! Get File-Name
  WRITE(CYEAR, FMT='(I4)') CurYear + 1
  fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
       '_climate_rst.nc'

  INQUIRE( FILE=TRIM( fname ), EXIST=EXISTFILE )

        IF ( .NOT.EXISTFILE) write(*,*) fname, ' does not exist!!'

  ! Open NetCDF file:
  STATUS = NF90_OPEN(fname, NF90_NOWRITE, FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


  ! dimensions:
  ! Land (number of points)
  STATUS = NF90_INQ_DIMID(FILE_ID, 'land'   , dID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  ! number of years (stored for 20 y running means
  STATUS = NF90_INQ_DIMID(FILE_ID, 'nyear'   , dID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  ! number of days (stored for 31 d monthly means
  STATUS = NF90_INQ_DIMID(FILE_ID, 'nday'   , dID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  IF ( land_dim .NE. SIZE(patch%latitude)) THEN
     WRITE(*,*) "Dimension misfit, ", fname
     WRITE(*,*) "land_dim", land_dim
     STOP
  ENDIF

  ! LAT & LON
  STATUS = NF90_INQ_VARID( FILE_ID, A1(1), dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_GET_VAR( FILE_ID, dID, LAT )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


  STATUS = NF90_INQ_VARID( FILE_ID, A1(2), dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_GET_VAR( FILE_ID, dID, LON )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

! READ scalar fields
  DO i = 1, SIZE(A0)
     STATUS = NF90_INQ_VARID( FILE_ID, A0(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI0 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A0(i)))
     CASE ('nyears'      ) ; climate%nyears      = TMPI0
     END SELECT
  END DO


! READ 1-dimensional real fields
  DO i = 3, SIZE(A1)
     STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A1(i)))
     CASE ('mtemp'      ) ; climate%mtemp       = TMP
     CASE ('mtemp_min'   ) ; climate%mtemp_min   = TMP
     CASE ('mtemp_max'  ) ; climate%mtemp_max  = TMP
     CASE ('mtemp_min20' ) ; climate%mtemp_min20 = TMP
     CASE ('mtemp_max20'  ) ; climate%mtemp_max20  = TMP
     CASE ('atemp_mean'  ) ; climate%atemp_mean  = TMP
     CASE ('AGDD5'  ) ; climate%AGDD5  = TMP
     CASE ('GDD5'  ) ; climate%GDD5  = TMP
     CASE ('GDD0'  ) ; climate%GDD0  = TMP
     END SELECT
  END DO

! READ 1-dimensional integer fields
  DO i = 1, SIZE(AI1)

     write(*,*)  TRIM(AI1(i))
     STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMPI )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(AI1(i)))
     CASE ('chilldays'      ) ; climate%chilldays      = TMPI
     END SELECT
  END DO


 ! READ 2-dimensional fields (nyear)
  DO i = 1, SIZE(A2)
     write(*,*)  TRIM(A2(i))
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('mtemp_min_20' ) ; climate%mtemp_min_20 = TMP2
     CASE ('mtemp_max_20' ) ; climate%mtemp_max_20 = TMP2
     END SELECT
  END DO

 ! READ 2-dimensional fields (nday)
  DO i = 1, SIZE(A3)
     STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A3(i)))
     CASE ('dtemp_31' ) ; climate%dtemp_31 = TMP3
     END SELECT
  END DO


  ! Close NetCDF file:
  STATUS = NF90_close(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


END SUBROUTINE  READ_CLIMATE_RESTART_NC


END MODULE cable_climate_mod
