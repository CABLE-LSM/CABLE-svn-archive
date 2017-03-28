!Ticket200:YP's version of write
SUBROUTINE WRITE_CASA_RESTART_NC ( casamet, casapool, casaflux, phen, CASAONLY )

  USE CASAVARIABLE, ONLY : casa_met, casa_pool, casa_flux, icycle, mplant, mlitter, msoil
  USE CABLE_COMMON_MODULE
  USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, mp
  USE phenvariable
  USE casavariable
  USE netcdf

  IMPLICIT NONE

 
  TYPE (casa_met),  INTENT(IN) :: casamet
  TYPE (casa_pool),  INTENT(IN) :: casapool
  TYPE (casa_flux),           INTENT(IN) :: casaflux
  TYPE (phen_variable),       INTENT(IN) :: phen

  INTEGER*4 :: mp4
  INTEGER*4, parameter   :: pmp4 =0
  INTEGER, parameter   :: fmp4 = kind(pmp4)
  INTEGER*4   :: STATUS
  INTEGER*4   :: FILE_ID, land_ID, plnt_ID, litt_ID, soil_ID, i
  LOGICAL   :: CASAONLY
  CHARACTER :: CYEAR*4, FNAME*99,dum*50

  ! ! 1 dim arrays (npt )
  ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
  !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
  ! ! 2 dim arrays (npt,mplant)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
  ! ! 2 dim arrays (npt,mlitter)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
  ! ! 2 dim arrays (npt,msoil)
  ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)

  ! 1 dim arrays (npt )
  CHARACTER(len=20),DIMENSION(12) :: A1
  CHARACTER(len=20),DIMENSION(2) :: AI1
  ! 2 dim arrays (npt,mplant)
  CHARACTER(len=20),DIMENSION(3) :: A2
  ! 2 dim arrays (npt,mlitter)
  CHARACTER(len=20),DIMENSION(3) :: A3
  ! 2 dim arrays (npt,msoil)
  CHARACTER(len=20),DIMENSION(3) :: A4
  INTEGER*4 :: VID1(SIZE(A1)), VIDI1(SIZE(AI1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))

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

  ! Get File-Name
  WRITE(CYEAR, FMT='(I4)') CurYear + 1

  !fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
  !     '_'//CYEAR//'_casa_rst.nc'
  fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
       '_casa_rst.nc'
  ! Create NetCDF file:
  STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
write(*,*) 'writing casa restart', fname 
  ! Put the file in define mode:
  STATUS = NF90_redef(FILE_ID)

  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )
  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle", icycle  )
  IF ( CASAONLY ) THEN
     dum = 'CASA-ONLY run'
  ELSE
     dum = 'CABLE-CASA coupled run'
  ENDIF
  STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

  ! Define dimensions:
  ! Land (number of points)
  STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  DO i = 1, SIZE(A1)
     STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(AI1)
     STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A2)
     STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID/),VID2(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A3)
     STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID/),VID3(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  DO i = 1, SIZE(A4)
     STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID/),VID4(i))
     IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  END DO

  ! End define mode:
  STATUS = NF90_enddef(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  ! PUT LAT / LON
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), casamet%lat )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), casamet%lon )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  ! PUT VARS
  STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), casamet%glai )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), casapool%clabile )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


  STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), casaflux%frac_sapwood )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), casaflux%sapwood_area )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), phen%phen )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), phen%aphen )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casapool%Nsoilmin )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), phen%phase )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), phen%doyphase(:,3) )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter  )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
  
  STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil )
  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

  IF (icycle ==3) then
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), casapool%psoillab )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), casapool%psoilsorb )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), casapool%psoilocc )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     

     STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant  )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
     
     STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter )
     IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
      
  ENDIF
  ! Close NetCDF file:
  STATUS = NF90_close(FILE_ID)
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

END SUBROUTINE WRITE_CASA_RESTART_NC


