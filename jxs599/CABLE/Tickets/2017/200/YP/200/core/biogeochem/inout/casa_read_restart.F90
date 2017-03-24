SUBROUTINE READ_CASA_RESTART_NC (  casamet, casapool, casaflux,phen )

  USE CASAVARIABLE
  USE phenvariable
  USE CABLE_COMMON_MODULE
  USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, r_2, mp
  USE netcdf

  IMPLICIT NONE

  !INTEGER, INTENT(in)    :: YEAR
  TYPE (casa_met) , INTENT(inout) :: casamet
  TYPE (casa_pool), INTENT(inout) :: casapool
  TYPE (casa_flux), INTENT(inout) :: casaflux
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS, i
  INTEGER*4   :: FILE_ID, dID, land_dim, mp_dim, ml_dim, ms_dim
  CHARACTER :: FRST_IN*99, CYEAR*4, CDATE*12, RSTDATE*12, FNAME*99

  ! ! 1 dim arrays (npt )
  ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
  !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
  ! ! 2 dim arrays (npt,mplant)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
  ! ! 2 dim arrays (npt,mlitter)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
  ! ! 2 dim arrays (npt,msoil)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)
  REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
  REAL(r_2)                         :: TMP2(mp,mplant),TMP3(mp,mlitter),TMP4(mp,msoil)

  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(12) :: A1
  CHARACTER(len=20),DIMENSION(2) :: AI1
  ! 2 dim arrays (npt,mplant)
  CHARACTER(len=20),DIMENSION(3) :: A2
  ! 2 dim arrays (npt,mlitter)
  CHARACTER(len=20),DIMENSION(3) :: A3
  ! 2 dim arrays (npt,msoil)
  CHARACTER(len=20),DIMENSION(3) :: A4
  INTEGER :: VID1(SIZE(A1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))
  LOGICAL            ::  EXISTFILE, EXISTFILE1
  mp4=int(mp,fmp4)
  A1(1) = 'latitude'
  A1(2) = 'longitude'
  A1(3) = 'glai'
  A1(4) = 'clabile'
  A1(5) = 'psoillab'
  A1(6) = 'psoilsorb'
  A1(7) = 'psoilocc'
  A1(8) = 'frac_sapwood'
  A1(9) = 'sapwood_area'
  A1(10) = 'phen'
  A1(11) = 'aphen'
  A1(12) = 'nsoilmin'

  AI1(1) = 'phase'
  AI1(2) = 'doyphase3'

  A2(1) = 'cplant'
  A2(2) = 'nplant'
  A2(3) = 'pplant'
  A3(1) = 'clitter'
  A3(2) = 'nlitter'
  A3(3) = 'plitter'
  A4(1) = 'csoil'
  A4(2) = 'nsoil'
  A4(3) = 'psoil'

!fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
!       '_casa_rst.nc'
  fname =  TRIM(casafile%cnpipool)
  INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
  IF (EXISTFILE) THEN
     STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     PRINT *, 'initial pool from restart file: ', fname
  ELSE
     write(*,*) 'CASA restart file:', TRIM(fname), ' does not exist'
     fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
          '_casa_rst.nc'   
     INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE1 )
     IF (EXISTFILE1) THEN
        STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
        IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
        PRINT *, 'initial pool from restart file: ', fname
     ELSE
        write(*,*) 'CASA restart file:', TRIM(fname), ' does not exist either'
        write(*,*) 'Set cable_user%CASA_fromZero to true to initialise without restart file.'
        write(*,*) 'Otherwise set casafile%cnpipool to netcdf restart file name in cable.nml'
        stop
     ENDIF
  ENDIF

  ! TIME
  STATUS = NF90_GET_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", RSTDATE )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$
  WRITE(CYEAR, FMT="(I4)") CurYear
  CDATE = '01/01/'//CYEAR
  ! compare current year with restart year (only for non-site type met data)
  IF ( CDATE .NE. RSTDATE .and. &
      TRIM(cable_user%MetType).NE.'' .and. TRIM(cable_user%MetType).NE.'site' ) THEN
     WRITE(*,*)"Restart Date in rst file doesn't match start date of Run!"
     WRITE(*,*)"File: "//RSTDATE//' Run: '//CDATE
    ! STOP
  ENDIF

  ! DIMS
  STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'mplant', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mp_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'mlitter', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ml_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'msoil', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ms_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  IF ( land_dim .NE. SIZE(casamet%lon) .OR. mp_dim .NE. mplant .OR. &
       ml_dim   .NE. mlitter             .OR. ms_dim .NE. msoil ) THEN
     WRITE(*,*)"Dimension misfit!"
     WRITE(*,*)"Restart file      Run"
     WRITE(*,*)"# points  ",land_dim,"     ",SIZE(casamet%lon)
     WRITE(*,*)"# mplant  ",mp_dim,"     ",mplant
     WRITE(*,*)"# mlitter ",ml_dim,"     ",mlitter
     WRITE(*,*)"# msoil   ",ms_dim,"     ",msoil
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

  ! CHECK FOR VALID LONS

  ! READ 1-dimensional fields
  DO i = 3, SIZE(A1)
     STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A1(i)))
     CASE ('glai'      ) ; casamet%glai       = TMP
     CASE ('clabile'   ) ; casapool%clabile   = TMP
     CASE ('frac_sapwood' ) ; casaflux%frac_sapwood  = TMP
     CASE ( 'sapwood_area' ) ; casaflux%sapwood_area  = TMP
     CASE ( 'phen' ) ; phen%phen  = TMP
     CASE ( 'aphen' ) ; phen%aphen  = TMP
     CASE ( 'nsoilmin' ) ; casapool%Nsoilmin  = TMP
     END SELECT
  END DO
IF (icycle==3) then
  DO i = 3, SIZE(A1)
     STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A1(i)))
     CASE ('psoillab'  ) ; casapool%psoillab  = TMP
     CASE ('psoilsorb' ) ; casapool%psoilsorb = TMP
     CASE ('psoilocc'  ) ; casapool%psoilocc  = TMP
     END SELECT
  END DO
ENDIF

  DO i = 1, SIZE(AI1)
     STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(AI1(i)))
     CASE ( 'phase' ) ; phen%phase  = TMP
     CASE ( 'doyphase3' ) ; phen%doyphase(:,3)  = TMP
     END SELECT
  END DO

  ! READ 2-dimensional fields (mplant)
  DO i = 1, SIZE(A2)
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('cplant' ) ; casapool%cplant = TMP2
     CASE ('nplant' ) ; casapool%nplant = TMP2
     END SELECT
  END DO


IF (icycle==3) then
   DO i = 1, SIZE(A2)
     STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A2(i)))
     CASE ('pplant' ) ; casapool%pplant = TMP2
     END SELECT
  END DO
ENDIF

  ! READ 2-dimensional fields (mlitter)
  DO i = 1, SIZE(A3)
     STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A3(i)))
     CASE ('clitter' ) ; casapool%clitter = TMP3
     CASE ('nlitter' ) ; casapool%nlitter = TMP3
     END SELECT
  END DO

IF (icycle==3) then

  DO i = 1, SIZE(A3)
     STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A3(i)))
     CASE ('plitter' ) ; casapool%plitter = TMP3
     END SELECT
  END DO


ENDIF

  ! READ 2-dimensional fields (msoil)
  DO i = 1, SIZE(A4)
     STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     SELECT CASE ( TRIM(A4(i)))
     CASE ('csoil' ) ; casapool%csoil = TMP4
     CASE ('nsoil' ) ; casapool%nsoil = TMP4
     END SELECT
  END DO
IF (icycle==3) then
 DO i = 1, SIZE(A4)
     STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
     STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

     SELECT CASE ( TRIM(A4(i)))
     CASE ('psoil' ) ; casapool%psoil = TMP4
     END SELECT
  END DO
ENDIF

  STATUS = NF90_CLOSE( FILE_ID )

END SUBROUTINE READ_CASA_RESTART_NC

