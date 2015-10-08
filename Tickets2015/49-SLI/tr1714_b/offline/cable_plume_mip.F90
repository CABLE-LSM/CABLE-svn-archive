MODULE CABLE_PLUME_MIP

  USE netcdf 
  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, LEAP_DAY, &
       HANDLE_ERR, GET_UNIT

  USE cable_IO_vars_module, ONLY: logn, land_x, land_y
  
  IMPLICIT NONE

  ! The whole PLUME information

  TYPE PLUME_MET_TYPE 
     REAL, DIMENSION(:), ALLOCATABLE :: VAL
  END TYPE PLUME_MET_TYPE

  TYPE PLUME_MIP_TYPE
     INTEGER  :: mland, NMET, xdimsize, ydimsize, tdimsize
     INTEGER  :: CYEAR, MetStart, MetEnd, CTSTEP, DT, ktau
     REAL,   DIMENSION(:)  ,ALLOCATABLE :: AVG_LWDN, CO2VALS
     LOGICAL,DIMENSION(:,:),ALLOCATABLE :: LandMask 
     CHARACTER(len=15) :: Run,Forcing,RCP, CO2, NDEP,RCPdir
     CHARACTER(len=200):: BasePath, MetPath, LandMaskFile
     INTEGER,           DIMENSION(9) :: F_ID, V_ID
     CHARACTER(len=12) ,DIMENSION(9) :: VAR_NAME
     CHARACTER(len=200),DIMENSION(9) :: MetFile
     TYPE(PLUME_MET_TYPE), DIMENSION(11) :: MET
  END TYPE PLUME_MIP_TYPE
  
  TYPE (PLUME_MIP_TYPE):: PLUME

  ! Some local parameters

  INTEGER, PRIVATE, PARAMETER :: & 
       prec     =  1, &
       snow     =  2, &
       lwdn     =  3, &
       swdn     =  4, &
       pres     =  5, &
       rhum     =  6, &
       tmax     =  7, &
       tmin     =  8, &
       wind     =  9, &
       prevTmax = 10, &
       nextTmin = 11

  INTEGER, PRIVATE :: STATUS

  REAL, PRIVATE, PARAMETER :: SecDay = 86400.

  CHARACTER(len=6),DIMENSION(9),PARAMETER,PRIVATE  :: &
       PREF = (/ "pr    ","prsn  ","rlds  ","rsds  ","ps    ","hurs  ","tasmax","tasmin","wind  " /)
      
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PLUME routines 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PLUME_MIP_INIT( PLUME )

  USE cable_IO_vars_module, ONLY: latitude, longitude, nmetpatches, &
       mask, metGrid, sdoy, smoy, syear, shod, xdimsize, ydimsize
  
  USE cable_def_types_mod,  ONLY: mland

  IMPLICIT NONE

  TYPE (PLUME_MIP_TYPE):: PLUME

  INTEGER              :: STATUS, iu
  INTEGER              :: FID, latID, lonID, landID, timID, tdimsize
  INTEGER              :: cnt, x, y
  LOGICAL              :: ERR = .FALSE.
  CHARACTER(len=15)    :: Run, Forcing, RCP, CO2, NDEP
  CHARACTER(len=200)   :: BasePath, LandMaskFile, LMFILE
  REAL                 :: DT
  REAL,DIMENSION(:)  ,ALLOCATABLE :: plume_lats, plume_lons
  INTEGER,DIMENSION(:,:),ALLOCATABLE :: landmask

  NAMELIST /PLUMENML/ BasePath, LandMaskFile, Run, Forcing, RCP, CO2, NDEP, DT

  ! Read PLUME settings

  PLUME%NMET = 9
  IF ( TRIM(PLUME%Forcing) .EQ. "watch") THEN
     PLUME%VAR_NAME(prec)  = "Rainf"
     PLUME%VAR_NAME(snow)  = "Snowf"
     PLUME%VAR_NAME(lwdn)  = "LWdown"
     PLUME%VAR_NAME(swdn)  = "SWdown"
     PLUME%VAR_NAME(pres)  = "PSurf"
     PLUME%VAR_NAME(rhum)  = "Qmean"
     PLUME%VAR_NAME(tmax)  = "Tmax"
     PLUME%VAR_NAME(tmin)  = "Tmin"
     PLUME%VAR_NAME(wind)  = "Wind"
  ELSE
     PLUME%VAR_NAME(prec)  = "prAdjust"
     PLUME%VAR_NAME(snow)  = "prsnAdjust"
     PLUME%VAR_NAME(lwdn)  = "rldsAdjust"
     PLUME%VAR_NAME(swdn)  = "rsdsAdjust"
     PLUME%VAR_NAME(pres)  = "psAdjust"
     PLUME%VAR_NAME(rhum)  = "hurs"
     PLUME%VAR_NAME(tmax)  = "tasmaxAdjust"
     PLUME%VAR_NAME(tmin)  = "tasminAdjust"
     PLUME%VAR_NAME(wind)  = "windAdjust"
  END IF

  CALL GET_UNIT(iu)
  OPEN (iu,FILE="plume.nml")
  READ (iu,NML=PLUMENML)
  CLOSE(iu)

  PLUME%BasePath     = BasePath
  PLUME%LandMaskFile = LandMaskFile
  PLUME%Run          = Run    
  PLUME%Forcing      = Forcing
  PLUME%RCP          = RCP    
  PLUME%CO2          = CO2    
  PLUME%NDEP         = NDEP   
  PLUME%DT           = DT * 3600.  ! in seconds
  ! Print settings
  
  WRITE(*   ,*)"========================================= PLUME ============"
  WRITE(*   ,*)"PLUME-MIP settings chosen:"
  WRITE(*   ,*)" BasePath: ",TRIM(PLUME%BasePath)
  WRITE(*   ,*)" LandMask: ",TRIM(PLUME%LandMaskFile)
  WRITE(*   ,*)" Run     : ",PLUME%Run    
  WRITE(*   ,*)" Forcing : ",PLUME%Forcing
  WRITE(*   ,*)" RCP     : ",PLUME%RCP    
  WRITE(*   ,*)" CO2     : ",PLUME%CO2    
  WRITE(*   ,*)" NDEP    : ",PLUME%NDEP   
  WRITE(*   ,*)" DT      : ",PLUME%DT   
  WRITE(logn,*)"========================================= PLUME ============"
  WRITE(logn,*)"PLUME-MIP settings chosen:"
  WRITE(logn,*)" BasePath: ",TRIM(PLUME%BasePath)
  WRITE(logn,*)" LandMask: ",TRIM(PLUME%LandMaskFile)
  WRITE(logn,*)" Run     : ",PLUME%Run    
  WRITE(logn,*)" Forcing : ",PLUME%Forcing
  WRITE(logn,*)" RCP     : ",PLUME%RCP    
  WRITE(logn,*)" CO2     : ",PLUME%CO2    
  WRITE(logn,*)" NDEP    : ",PLUME%NDEP   
  WRITE(logn,*)" DT      : ",PLUME%DT   

PRINT*,"DOUBLECHECK FOR EACH CASE IN FILE_SWITCH WHETHER ANNUAL OR OTHER!"

  ! check for valid Run Identifier

  SELECT CASE (TRIM(PLUME%Run))
  CASE( "spinup","1850_1900","2006_2099" ) ; CONTINUE
  CASE( "1901_2001" )
     IF ( TRIM(PLUME%Forcing) .NE. "watch" ) THEN
        WRITE(*   ,*)"Run 1901_2001 must use 'watch' forcing!"
        WRITE(logn,*)"Run 1901_2001 must use 'watch' forcing!"
        ERR = .TRUE.
     ENDIF
  CASE( "1901_2005" )
     IF ( TRIM(PLUME%Forcing) .EQ. "watch" ) THEN
        WRITE(*   ,*)"Run 1901_2005 cannot use 'watch' forcing!"
        WRITE(logn,*)"Run 1901_2005 cannot use 'watch' forcing!"
        ERR = .TRUE.
     ENDIF
  CASE  default 
     WRITE(*   ,*)"Wrong PLUME%Run: ",PLUME%Run
     WRITE(*   ,*)"Use: Spinup, 1850_1900, 1901_2001, 1901_2005 or 2006_2099!"
     WRITE(logn,*)"Wrong PLUME%Run: ",PLUME%Run
     WRITE(logn,*)"Use: Spinup, 1850_1900, 1901_2001, 1901_2005 or 2006_2099!"
     ERR = .TRUE.
  END SELECT

  ! check for valid Forcing

  SELECT CASE (TRIM(PLUME%Forcing))
  CASE ("hadgem2-es", "ipsl-cm5a-lr", "miroc-esm-chem", &
        "gfdl-esm2m", "noresm1-m","ccsm4", "watch" ) ; CONTINUE
  CASE default
     WRITE(*   ,*)"Wrong PLUME%Forcing: ",PLUME%Forcing
     WRITE(*   ,*)"Please choose any of"
     WRITE(*   ,*)" hadgem2-es"
     WRITE(*   ,*)" ipsl-cm5a-lr"
     WRITE(*   ,*)" miroc-esm-chem"
     WRITE(*   ,*)" gfdl-esm2m"
     WRITE(*   ,*)" noresm1-m"
     WRITE(*   ,*)" ccsm4"
     WRITE(*   ,*)" watch" 
     WRITE(logn,*)"Wrong PLUME%Forcing: ",PLUME%Forcing
     WRITE(logn,*)"Please choose any of"
     WRITE(logn,*)" hadgem2-es"
     WRITE(logn,*)" ipsl-cm5a-lr"
     WRITE(logn,*)" miroc-esm-chem"
     WRITE(logn,*)" gfdl-esm2m"
     WRITE(logn,*)" noresm1-m"
     WRITE(logn,*)" ccsm4"
     WRITE(logn,*)" watch" 
     ERR = .TRUE.
  END SELECT

  ! check for valid RCP

  SELECT CASE (TRIM(PLUME%RCP))
  CASE ("hist","2.6","4.5","6.0","8.5") 
     
  CASE default
     WRITE(*   ,*)"Wrong PLUME%RCP: ",PLUME%RCP
     WRITE(*   ,*)"Please choose any of"
     WRITE(*   ,*)" hist, 2.6, 4.5, 6.0, 8.5"
     WRITE(logn,*)"Wrong PLUME%RCP: ",PLUME%RCP
     WRITE(logn,*)"Please choose any of"
     WRITE(logn,*)" hist, 2.6, 4.5, 6.0, 8.5"
     ERR = .TRUE.
  END SELECT

  IF ( ERR ) THEN
     WRITE(logn,*)"Invalid settings in PLUME_INIT"
     STOP "Invalid settings in PLUME_INIT"
  ENDIF

  ! Determine paths to met-files re settings
  IF ( TRIM(PLUME%Forcing) .EQ. "watch" ) THEN
     PLUME%MetPath = TRIM(PLUME%BasePath)//"/WATCH/"   
  ELSE
     PLUME%MetPath = TRIM(PLUME%BasePath)//"/GCM/"//TRIM(PLUME%Forcing)//"/"
  ENDIF
   
  IF (TRIM(PLUME%Run) .EQ. "spinup" .OR. TRIM(PLUME%Run) .EQ. "1850_1900") THEN
     PLUME%MetPath = TRIM(PLUME%MetPath)//"spinup_data"//"/hist/"
  ELSE
     SELECT CASE(TRIM(PLUME%RCP))
     CASE ( "hist" ); PLUME%RCPdir = "hist"   
     CASE ( "2.6"  ); PLUME%RCPdir = "rcp2p6"
     CASE ( "4.5"  ); PLUME%RCPdir = "rcp4p5"
     CASE ( "6.0"  ); PLUME%RCPdir = "rcp6p0"
     CASE ( "8.5"  ); PLUME%RCPdir = "rcp8p5"  
     END SELECT
     PLUME%MetPath = TRIM(PLUME%MetPath)//TRIM(PLUME%RCPdir)//"/"
  ENDIF
  
  WRITE(*   ,*)"========================================= PLUME ============"
  WRITE(logn,*)"========================================= PLUME ============"
  
  ! Now read landmask file
  ! Landmask file into init! Get LAt, LON etc. from there
!  LMFILE = TRIM(PLUME%BasePath)//"/plumber_landmask_1x1.nc"
!  LMFILE = TRIM(PLUME%BasePath)//"/plumber_landmask_.5x.5_F.nc"
  LMFILE = TRIM(PLUME%BasePath)//"/"//TRIM(PLUME%LandMaskFile)
  WRITE(*   ,*) 'Opening PLUME landmask file: ',TRIM(LMFILE)
  WRITE(logn,*) 'Opening PLUME landmask file: ',TRIM(LMFILE) 
  
  STATUS = NF90_OPEN(TRIM(LMFILE), NF90_NOWRITE, FID)
  CALL HANDLE_ERR(STATUS, "Opening PLUME Land-mask file"//TRIM(LMFILE))

  ! lat dimension
  STATUS = NF90_INQ_DIMID(FID,'latitude',latID)
  STATUS = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
  CALL HANDLE_ERR(STATUS, "Inquiring 'lat'"//TRIM(LMFILE))
  PLUME%ydimsize = ydimsize

  ALLOCATE( plume_lats ( ydimsize ) )
  STATUS = NF90_INQ_VARID(FID,'latitude',latID)
  CALL HANDLE_ERR(STATUS, "Inquiring 'latitudes'"//TRIM(LMFILE))
  STATUS = NF90_GET_VAR(FID,latID,plume_lats)
  CALL HANDLE_ERR(STATUS, "Reading 'latitudes'"//TRIM(LMFILE))
  
  ! lon dimension
  STATUS = NF90_INQ_DIMID(FID,'longitude',lonID)
  STATUS = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
  CALL HANDLE_ERR(STATUS, "Inquiring 'lon'"//TRIM(LMFILE))
  PLUME%xdimsize = xdimsize

  ALLOCATE( plume_lons ( xdimsize ) )
  STATUS = NF90_INQ_VARID(FID,'longitude',lonID)
  CALL HANDLE_ERR(STATUS, "Inquiring 'longitudes'"//TRIM(LMFILE))
  STATUS = NF90_GET_VAR(FID,lonID,plume_lons)
  CALL HANDLE_ERR(STATUS, "Reading 'longitudes'"//TRIM(LMFILE))
  
  ! Allocate PLUME arrays
  ALLOCATE( PLUME%landmask ( xdimsize, ydimsize) )
  ALLOCATE( landmask ( xdimsize, ydimsize) )
  
  ! get mask
  STATUS = NF90_INQ_VARID(FID,'land',landID)
  CALL HANDLE_ERR(STATUS, "Inquiring 'land' "//TRIM(LMFILE))
  STATUS = NF90_GET_VAR(FID,landID,landmask)
  CALL HANDLE_ERR(STATUS, "Reading 'land' "//TRIM(LMFILE))

  ! Number of land cells
  WHERE ( landmask .GT. 0 ) 
     PLUME%landmask = .TRUE.
     mask           = 1
  ELSEWHERE
     PLUME%landmask = .FALSE.
     mask           = 0
  END WHERE
  PLUME%mland = COUNT(PLUME%landmask)

  ! Allocate CABLE arrays
  ALLOCATE( latitude(PLUME%mland), longitude(PLUME%mland) )
  ALLOCATE( land_y  (PLUME%mland), land_x   (PLUME%mland) )
  DO x = 1, PLUME%NMET
     ALLOCATE( PLUME%MET(x)%VAL(PLUME%mland) )   
  END DO
  ! Extra fields for Tmin and Tmax for Weathergenerator
  ALLOCATE( PLUME%MET(prevTmax)%VAL(PLUME%mland) )
  ALLOCATE( PLUME%MET(nextTmin)%VAL(PLUME%mland) )

  ! Map all to Landgrid arrays
  cnt = 1
  DO y = 1, ydimsize
     DO x = 1, xdimsize
        IF ( .NOT. PLUME%landmask(x,y) ) CYCLE
PRINT*," lo,la, x,y   ",plume_lons(x),plume_lats(y),x, y

        land_x   (cnt) = x
        land_y   (cnt) = y
        longitude(cnt) = plume_lons(x)
        latitude (cnt) = plume_lats(y)
        cnt = cnt + 1
     END DO
  END DO

  ! set global cable variables
  metGrid     = "mask"
  ALLOCATE( mask(xdimsize, ydimsize) )
  mask        = landmask
  mland       = PLUME%mland
  nmetpatches = 1
  ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
  shod        = 0.
  sdoy        = 1
  smoy        = 1
  syear       = PLUME%CYEAR
  ! and a PLUME var
  ALLOCATE( PLUME%AVG_LWDN(mland) )
  
  DEALLOCATE ( landmask, plume_lats, plume_lons )

  STATUS = NF90_CLOSE(FID)
  CALL HANDLE_ERR(STATUS, "Closing mask-file"//TRIM(LMFILE))

END SUBROUTINE PLUME_MIP_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PLUME_GET_FILENAME ( PLUME, cyear, par, FN )

  IMPLICIT NONE
  
  TYPE(PLUME_MIP_TYPE) :: PLUME
  INTEGER,              INTENT(IN)  :: cyear, par
  CHARACTER(LEN=200),   INTENT(OUT) :: FN
  INTEGER                           :: i, idx
  CHARACTER            :: cy*4, sfy*12, mp*200,fc*15,rcp*15 
  
  INTEGER,         DIMENSION(21),PARAMETER :: &
       syear = (/ 1901, 1911, 1921, 1931, 1941, 1951, 1961, 1971, 1981, 1991, 2001, &
       2006, 2011, 2021, 2031, 2041, 2051, 2061, 2071, 2081, 2091 /) , &
       eyear = (/ 1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2005, &
       2010, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2099 /) 
  
  FN = "                                                                    "

  WRITE(cy,FMT='(I4)')cyear
  mp = PLUME%MetPath
  fc = PLUME%Forcing
  rcp= PLUME%RCP
     
  IF ( TRIM(fc) .EQ. "watch" ) THEN
     ! WATCH data comes in annual files
     IF ( TRIM(PLUME%Run) .EQ. "Spinup" ) THEN
        SELECT CASE ( par )
        CASE(prec) ; FN = TRIM(mp)//"/1901_1930/Rainf_daily_WFD_TYX_format_detrended_" //cy//".nc"
        CASE(snow) ; FN = TRIM(mp)//"/1901_1930/Snowf_daily_WFD_TYX_format_detrended_" //cy//".nc"
        CASE(lwdn) ; FN = TRIM(mp)//"/1901_1930/LWdown_daily_WFD_TYX_format_detrended_"//cy//".nc"
        CASE(swdn) ; FN = TRIM(mp)//"/1901_1930/SWdown_daily_WFD_TYX_format_detrended_"//cy//".nc"
        CASE(pres) ; FN = TRIM(mp)//"/1901_1930/PSurf_daily_WFD_TYX_format_detrended_" //cy//".nc"
        CASE(rhum) ; FN = TRIM(mp)//"/1901_1930/Qmean_WFD_TYX_format_detrended_"       //cy//".nc"
        CASE(tmax,PrevTmax) ; FN = TRIM(mp)//"/1901_1930/Tmax_WFD_TYX_format_detrended_"        //cy//".nc"
        CASE(tmin,NextTmin) ; FN = TRIM(mp)//"/1901_1930/Tmin_WFD_TYX_format_detrended_"        //cy//".nc"
        CASE(wind) ; FN = TRIM(mp)//"/1901_1930/Wind_daily_WFD_TYX_format_"            //cy//".nc"
        END SELECT
     ELSE IF ( TRIM(PLUME%Run) .EQ. "1850_1900" ) THEN
        STOP "Not yet implemented! PLUME: GET_FILE_NAMES"
     ELSE
        SELECT CASE ( par )
        CASE(prec) ; FN = TRIM(mp)//"/Rainf_daily_WFD/Rainf_daily_WFD_TYX_format_"  //cy//".nc"
        CASE(snow) ; FN = TRIM(mp)//"/Snowf_daily_WFD/Snowf_daily_WFD_TYX_format_"  //cy//".nc"
        CASE(lwdn) ; FN = TRIM(mp)//"/LWdown_daily_WFD/LWdown_daily_WFD_TYX_format_"//cy//".nc"
        CASE(swdn) ; FN = TRIM(mp)//"/SWdown_daily_WFD/SWdown_daily_WFD_TYX_format_"//cy//".nc"
        CASE(pres) ; FN = TRIM(mp)//"/PSurf_daily_WFD/PSurf_daily_WFD_TYX_format_"  //cy//".nc"
        CASE(rhum) ; FN = TRIM(mp)//"/Qmean_daily_WFD/Qmean_daily_WFD_TYX_format_"  //cy//".nc"
        CASE(tmax,PrevTmax) ; FN = TRIM(mp)//"/Tmax_daily_WFD/Tmax_daily_WFD_TYX_format_"    //cy//".nc"
        CASE(tmin,NextTmin) ; FN = TRIM(mp)//"/Tmin_daily_WFD/Tmin_daily_WFD_TYX_format_"    //cy//".nc"
        CASE(wind) ; FN = TRIM(mp)//"/Wind_daily_WFD/Wind_daily_WFD_TYX_format_"    //cy//".nc"
        END SELECT
     ENDIF
  ELSE
     ! find proper file for current time
        ! hist spinup only
     IF ( TRIM(PLUME%Run) .EQ. "Spinup" ) THEN
        FN = TRIM(mp)//"/"//TRIM(PREF(par))
        IF ( par .NE. rhum ) FN = TRIM(FN)//"Adjust"
        FN = TRIM(FN)//"_"//TRIM(fc)//"_"//TRIM(rcp)//"_"
        IF ( par .NE. wind ) FN = TRIM(FN)//"detrended_"
        FN = TRIM(FN)//cy
        
     ELSE IF ( TRIM(PLUME%Run) .EQ. "1850_1900" ) THEN
        STOP "Not yet implemented! PLUME: GET_FILE_NAMES"
     ELSE 
        ! real runs
        IF ( cyear .LT. syear(1) .OR. cyear .GT. 2099) THEN
           WRITE(*   ,*)"Wrong year (Must be 1901 <= y <= 2099) ",CYEAR
           WRITE(logn,*)"Wrong year (Must be 1901 <= y <= 2099) ",CYEAR
           STOP "Error in PLUME_GET_FILENAME"
        ENDIF
        idx = SIZE(syear)
        DO i = 1, SIZE(syear)
           IF ( syear(i) .GT. cyear ) THEN
              idx = i - 1
              EXIT
           ENDIF
        END DO
        WRITE(sfy,FMT='(I4,A1,I4,A3)')syear(idx),"-",eyear(idx),".nc"
        PLUME%MetStart = syear(idx)
        PLUME%MetEnd   = eyear(idx)
        FN = TRIM(mp)//"/"//TRIM(PREF(par))//"_"
        IF ( .NOT. par .EQ. rhum ) FN = TRIM(FN)//"bced_1960_1999_"
        FN = TRIM(FN)//TRIM(fc)//"_"//TRIM(PLUME%RCPdir)//"_"//sfy
     END IF
  END IF

END SUBROUTINE PLUME_GET_FILENAME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION FILE_SWITCH ( PLUME, action )

  IMPLICIT NONE

  TYPE ( PLUME_MIP_TYPE ), INTENT(IN) :: PLUME
  CHARACTER(len=5)  , INTENT(IN) :: action
  LOGICAL :: FILE_SWITCH 
  
  FILE_SWITCH = .FALSE.

  IF ( TRIM(PLUME%Forcing) .EQ. "watch" ) THEN
     FILE_SWITCH = .TRUE.
     RETURN
  ENDIF

  IF ( INDEX( action,"OPEN") .GE. 1 ) THEN
     !Open files 
     IF (TRIM(PLUME%Run) .EQ. "spinup" ) THEN
        IF ( PLUME%CYEAR .EQ. 1900 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "1850_1900" ) THEN
!WRONG like spinuup
        IF ( PLUME%CYEAR .EQ. 1850 .OR. PLUME%CYEAR .EQ. 1880 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "1901_2005" ) THEN
        IF ( MOD(PLUME%CYEAR-1,10) .EQ. 0 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "2006_2099" ) THEN
        IF ( MOD(PLUME%CYEAR-1,10) .EQ. 0 .OR. PLUME%CYEAR .EQ. 2006 ) FILE_SWITCH = .TRUE. 
      
     ENDIF
  ELSE IF ( INDEX(action, "CLOSE") .GE. 1 ) THEN
     ! Closing files
     IF (TRIM(PLUME%Run) .EQ. "spinup" ) THEN
        IF ( PLUME%CYEAR .EQ. 1930 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "1850_1900" ) THEN
        IF ( PLUME%CYEAR .EQ. 1879 .OR. PLUME%CYEAR .EQ. 1900 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "1901_2005" ) THEN
        IF ( MOD(PLUME%CYEAR  ,10) .EQ. 0 ) FILE_SWITCH = .TRUE.       

     ELSE IF ( TRIM(PLUME%Run) .EQ. "2006_2099" ) THEN
        IF ( MOD(PLUME%CYEAR  ,10) .EQ. 0 .OR. PLUME%CYEAR .EQ. 2006 ) FILE_SWITCH = .TRUE. 
      
     ENDIF
  ELSE 
     STOP "Wrong action in PLUME FILE_SWITCH! <OPEN|CLOSE> "
  ENDIF

END FUNCTION FILE_SWITCH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GET_PLUME_CO2( PLUME, CO2air )

  IMPLICIT NONE
  
  TYPE(PLUME_MIP_TYPE) :: PLUME
  REAL, INTENT(OUT)    :: CO2air

  INTEGER              :: i, iu, f, IOS = 0
  CHARACTER            :: CO2FILE*200

  IF ( TRIM(PLUME%Run) .EQ. "Spinup" ) THEN
     ! fixed 1850 value 
     CO2air = 284.72501
  ELSE IF ( TRIM(PLUME%CO2) .EQ. "static1990" ) THEN
     ! fixed 1990 value 
     CO2air = 353.85501
  ELSE IF ( TRIM(PLUME%CO2) .EQ. "static2085" ) THEN
     ! fixed 2085 value 
     STOP " static 2085 not yet implemented! cable_plume_mip.F90"
  
  ELSE 
     
     ! In 2006 new file must be opened -> old array discarded

     IF ( ALLOCATED( PLUME%CO2VALS ) .AND. PLUME%CYEAR .EQ. 2006 ) &
          DEALLOCATE ( PLUME%CO2VALS )      

     ! Read from file if none availale

     IF ( .NOT. ALLOCATED( PLUME%CO2VALS) ) THEN
        IF ( PLUME%CYEAR .LE. 2005 ) THEN
           ALLOCATE( PLUME%CO2VALS( 1850:2005 ) )
           CO2FILE = TRIM(PLUME%BasePath)//"/CO2/co2_1850_2005_hist.dat"
        ELSE
           ALLOCATE( PLUME%CO2VALS( 2006:2100 ) )
           CO2FILE = TRIM(PLUME%BasePath)//"/CO2/co2_2006_2100_"
           SELECT CASE(TRIM(PLUME%RCP))
           CASE ( "2.6"  ); CO2FILE = TRIM(CO2FILE)//"rcp2p6.dat"
           CASE ( "4.5"  ); CO2FILE = TRIM(CO2FILE)//"rcp4p5.dat"
           CASE ( "6.0"  ); CO2FILE = TRIM(CO2FILE)//"rcp6p0.dat"
           CASE ( "8.5"  ); CO2FILE = TRIM(CO2FILE)//"rcp8p5.dat"  
           END SELECT
        ENDIF
        
        ! open CO2 file and read

        CALL GET_UNIT(iu)
        OPEN (iu, FILE=TRIM(CO2FILE), STATUS="OLD", ACTION="READ")
        DO WHILE( IOS .EQ. 0 )
           READ(iu, FMT=*, IOSTAT=IOS)f,PLUME%CO2VALS(f)
        END DO
        CLOSE(iu)
        
     ENDIF

     CO2air = PLUME%CO2VALS( PLUME%CYEAR ) 
  ENDIF

END SUBROUTINE GET_PLUME_CO2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE OPEN_PLUME_MET( PLUME )

  USE cable_IO_vars_module, ONLY: timeunits

  IMPLICIT NONE

  TYPE( PLUME_MIP_TYPE )   :: PLUME
  INTEGER             :: i, iy, yy, tID
  LOGICAL, SAVE       :: CALL1 = .TRUE.

  ! find time 
! ipsl_cm5a_lr/rcp4p5/hurs_ipsl-cm5a-lr_rcp4p5_2006-2010.nc                                                          
! ipsl_cm5a_lr/rcp4p5/pr_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2010.nc                                             
! Wind_daily_WFD_TYX_format_1901.nc


  DO i = 1, PLUME%NMET
     CALL PLUME_GET_FILENAME( PLUME, PLUME%CYEAR, i, PLUME%MetFile(i) )

     ! OPEN NEW MET FILES and access variables
     WRITE(*   ,*) 'Opening met data file: ', PLUME%MetFile(i)
     WRITE(logn,*) 'Opening met data file: ', PLUME%MetFile(i)
     STATUS = NF90_OPEN(TRIM(PLUME%MetFile(i)), NF90_NOWRITE, PLUME%F_ID(i))
     CALL HANDLE_ERR(STATUS, "Opening PLUME file "//PLUME%MetFile(i) )
     STATUS = NF90_INQ_VARID(PLUME%F_ID(i),TRIM(PLUME%VAR_NAME(i)), PLUME%V_ID(i))
     CALL HANDLE_ERR(STATUS, "Inquiring PLUME var in "//PLUME%MetFile(i) )
  END DO

  ! Set internal counter
  PLUME%CTSTEP = 1
  ! For first call there might be an offset
  IF ( PLUME%CYEAR .GT. PLUME%MetStart ) THEN 
     DO yy = PLUME%MetStart, PLUME%CYEAR - 1
        PLUME%CTSTEP = PLUME%CTSTEP + 365 + LEAP_DAY( yy )
     END DO
  END IF

  ! Get unit of time
  IF ( CALL1 ) THEN
     STATUS = NF90_INQ_VARID(PLUME%F_ID(prec),'time',tID)
     CALL HANDLE_ERR(STATUS, "Inquiring PLUME time  in "//PLUME%MetFile(i) )
     STATUS = NF90_GET_ATT(PLUME%F_ID(prec),tID,'units',timeunits)
     CALL HANDLE_ERR(STATUS, "Inquiring PLUME timeunit in "//PLUME%MetFile(i) )
  ENDIF

  CALL1 = .FALSE.
  
END SUBROUTINE OPEN_PLUME_MET
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PLUME_GET_DAILY_MET( PLUME, TminFlag, islast )

  IMPLICIT NONE

  TYPE(PLUME_MIP_TYPE) :: PLUME
  LOGICAL, INTENT(IN)  :: TminFlag, islast
  REAL    :: tmparr(720,360)
  INTEGER :: t, i, ii, k, realk
  INTEGER :: fid, vid, tid
  INTEGER :: xds, yds, tds
  LOGICAL, SAVE :: CALL1 = .TRUE. 
  CHARACTER(LEN=200) :: filename


  ! Set previous day's Tmax to current (last step's) Tmax
  IF ( .NOT. CALL1 ) THEN
     PLUME%MET(prevTmax)%VAL(:) = PLUME%MET(  Tmax  )%VAL(:) 
     PLUME%MET(  Tmin  )%VAL(:) = PLUME%MET(NextTmin)%VAL(:) 
  ENDIF

  xds = PLUME%xdimsize
  yds = PLUME%ydimsize

  DO i= 1, PLUME%NMET
     
     IF ( i .EQ. Tmin )THEN 
        ii = nextTmin
        
        t  = PLUME%CTSTEP + 1
        
        ! On first occasion read t=1 as well  
        IF ( CALL1 ) THEN
           STATUS = NF90_GET_VAR(PLUME%F_ID(i), PLUME%V_ID(i), tmparr, &
                start=(/1,1,t-1/),count=(/xds,yds,1/) )
           CALL HANDLE_ERR(STATUS, "Reading from "//PLUME%MetFile(i) )

           DO k = 1, PLUME%mland
              PLUME%MET(i)%VAL(k) = tmparr( land_x(k), land_y(k) )
           END DO
        END IF
           
     ELSE
        ii = i
        t  = PLUME%CTSTEP
     ENDIF

     IF ( i .EQ. Tmin .AND. TminFlag ) THEN
        IF ( .NOT. islast ) THEN
           ! Open next file for Quick access or use same if very last t-step.
           t  = 1
           CALL PLUME_GET_FILENAME( PLUME, PLUME%CYEAR+1, Tmin, filename )
           STATUS = NF90_OPEN(TRIM(filename), NF90_NOWRITE, fid)
           CALL HANDLE_ERR(STATUS, "Opening PLUME file "//filename )
           
           STATUS = NF90_INQ_VARID(fid,TRIM(PLUME%VAR_NAME(i)), vid)
           CALL HANDLE_ERR(STATUS, "Inquiring PLUME var "//filename )
           
           STATUS = NF90_GET_VAR(fid, vid, tmparr, &
                start=(/1,1,t/),count=(/xds,yds,1/) )
           CALL HANDLE_ERR(STATUS, "Reading from "//filename )
           
           DO k = 1, PLUME%mland
              PLUME%MET(ii)%VAL(k) = tmparr( land_x(k), land_y(k) )
           END DO
           
           STATUS = NF90_CLOSE(fid)
           CALL HANDLE_ERR(STATUS, "Closing PLUME file "//filename)
           
        ELSE
           CONTINUE ! The values simply remain the same
        ENDIF
     ELSE

        ! STANDARD READ 
        ! variables from open files
        STATUS = NF90_GET_VAR(PLUME%F_ID(i), PLUME%V_ID(i), tmparr, &
             start=(/1,1,t/),count=(/xds,yds,1/) )
        CALL HANDLE_ERR(STATUS, "Reading from "//PLUME%MetFile(i) )
        DO k = 1, PLUME%mland
           PLUME%MET(ii)%VAL(k) = tmparr( land_x(k), land_y(k) )
        END DO
     ENDIF

     IF ( (i .EQ. Tmax ) .AND. CALL1 ) THEN
        ii = prevTmax

        IF ( PLUME%CYEAR .GT. 1900 ) THEN
           ! on
           CALL PLUME_GET_FILENAME( PLUME, PLUME%CYEAR-1, i, filename )
           STATUS = NF90_OPEN(TRIM(filename), NF90_NOWRITE, fid)
           CALL HANDLE_ERR(STATUS, "Opening PLUME file "//filename )
           
           STATUS = NF90_INQ_DIMID(fid,'time',tid)
           STATUS = NF90_INQUIRE_DIMENSION(fid,tid,len=tds)
           CALL HANDLE_ERR(STATUS, "Inquiring 'time'"//TRIM(filename))
           
           STATUS = NF90_INQ_VARID(fid,TRIM(PLUME%VAR_NAME(i)), vid)
           CALL HANDLE_ERR(STATUS, "Inquiring PLUME var "//filename )
           
           STATUS = NF90_GET_VAR(fid, vid, tmparr, &
                start=(/1,1,tds/),count=(/xds,yds,1/) )
           CALL HANDLE_ERR(STATUS, "Reading from "//filename )
           
           DO k = 1, PLUME%mland
              PLUME%MET(ii)%VAL(k) = tmparr( land_x(k), land_y(k) )
           END DO
           
           STATUS = NF90_CLOSE(fid)
           CALL HANDLE_ERR(STATUS, "Closing PLUME file "//filename)
           
        ELSE
           PLUME%MET(ii)%VAL(:) = PLUME%MET(i) %VAL(:)
        ENDIF
     END IF
        
        
        
  END DO

  ! Convert pressure Pa -> hPa

  PLUME%MET(pres)%VAL(:) = PLUME%MET(pres)%VAL(:) / 100.
  
  ! update internal counter

  PLUME%CTSTEP = PLUME%CTSTEP + 1

  ! Wrap up
  IF ( CALL1 )CALL1 = .FALSE.

END SUBROUTINE PLUME_GET_DAILY_MET

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PLUME_MIP_GET_MET(PLUME, MET, CurYear, ktau, kend, islast )
!==============================================================================
!
! Name: open_met_file
!
! Purpose: 
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
! CALLs: abort
!        nc_abort
!        date_and_time
!
! Input file: plume_landmask.nc
!             [GSWP_Snowf].nc
!             [GSWP_LWdown].nc
!             [GSWP_SWdown].nc
!             [GSWP_PSurf].nc
!             [GSWP_Qair].nc
!             [GSWP_Tair].nc
!             [GSWP_wind].nc
!             [GSWP_Rainf].nc
!
!==============================================================================

  USE cable_def_types_mod,   ONLY: MET_TYPE
  USE cable_IO_vars_module,  ONLY: LANDPT, latitude
  USE cable_common_module,   ONLY: DOYSOD2YMDHMS
  USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
                                   WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
  USE cable_checks_module,   ONLY: rh_sh

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: CurYear, ktau, kend
  LOGICAL, INTENT(IN)  :: islast

  TYPE(PLUME_MIP_TYPE) :: PLUME
  TYPE(MET_TYPE)       :: MET

  LOGICAL   :: newday                  
  INTEGER   :: i, dY, dM, dD
  REAL      :: dt, CO2air
  CHARACTER :: LMFILE*200

  TYPE(WEATHER_GENERATOR_TYPE), SAVE :: WG
  LOGICAL,                      SAVE :: CALL1 = .TRUE.

  dt = PLUME%DT

  ! On first step read and check PLUME settings and read land-mask 

  IF ( CALL1 ) CALL WGEN_INIT( WG, PLUME%mland, latitude, dt )

  ! time-step relevant settings

  PLUME%CYEAR = CurYear
  PLUME%ktau  = ktau
  
  !!!!  this only works with CANBERRA cable_driver, as ktau    !!!!
  !!!!  restarts on Jan 1                                      !!!!

  met%hod (landpt(1)%cstart) = REAL(MOD( (ktau-1) * NINT(dt), INT(SecDay)) ) / 3600.
  met%doy (landpt(:)%cstart) = INT(REAL(ktau-1) * dt / SecDay ) + 1
  met%year(landpt(:)%cstart) = Curyear

  CALL DOYSOD2YMDHMS(CurYear, INT(met%doy(1)), INT(met%hod(1)) * 3600, dY, dM, dD )

  met%moy (landpt(:)%cstart) = dM

  newday = ( met%hod(landpt(1)%cstart).EQ. 0 ) 

  ! Beginning-of-year accounting
  IF ( ktau .EQ. 1 ) THEN 

     ! Update CO2
     CALL GET_PLUME_CO2( PLUME, CO2air )

     met%ca(:) = CO2air / 1.e+6

     ! Open/close Met-files if necessary
     IF (FILE_SWITCH( PLUME, 'OPEN ' ) .OR. CALL1)  CALL OPEN_PLUME_MET( PLUME )

  ENDIF

  ! Now get the Met data for this day

  IF ( newday ) THEN

     CALL PLUME_GET_DAILY_MET( PLUME, (ktau.EQ.kend-((SecDay/dt)-1) .AND. FILE_SWITCH( PLUME, 'CLOSE' )), islast )

     ! Air pressure assumed to be constant over day

     met%pmb = PLUME%MET(pres)%VAL !CLN interpolation??

     WG%WindDay        = PLUME%MET(  wind  )%VAL
     WG%TempMinDay     = PLUME%MET(  Tmin  )%VAL - 273.15
     WG%TempMaxDay     = PLUME%MET(  Tmax  )%VAL - 273.15   
     WG%TempMinDayNext = PLUME%MET(NextTmin)%VAL - 273.15   
     WG%TempMaxDayPrev = PLUME%MET(PrevTmax)%VAL - 273.15   
     WG%SolarMJDay     = PLUME%MET(  swdn  )%VAL * 1.e-6 * SecDay ! ->[MJ/m2/d]
     WG%PrecipDay      = PLUME%MET(  prec  )%VAL * SecDay / 1000. ! ->[m/d]
     WG%SnowDay        = PLUME%MET(  snow  )%VAL * SecDay / 1000. ! ->[m/d]

     CALL WGEN_DAILY_CONSTANTS( WG, PLUME%mland, INT(met%doy(1))+1 )
     
     ! To get the diurnal cycle for lwdn get whole day and scale with 
     ! LWDN from file later 

     PLUME%AVG_LWDN(:) = 0.
     DO i = 1, NINT(SecDay/dt) 
        CALL WGEN_SUBDIURNAL_MET( WG, PLUME%mland, i-1 )
        PLUME%AVG_LWDN = PLUME%AVG_LWDN + WG%PhiLD
     END DO
     PLUME%AVG_LWDN = PLUME%AVG_LWDN / (SecDay/dt)
        
  END IF
  
  ! Decision has been made, that first tstep of the day is at 0:01 am

  CALL WGEN_SUBDIURNAL_MET( WG, PLUME%mland, NINT(met%hod(1)*3600./dt) )

  ! assign to cable variables

  met%precip    = WG%Precip 
  met%precip_sn = WG%Snow
  met%fld       = WG%PhiLD * PLUME%MET(lwdn )%VAL / PLUME%AVG_LWDN
  met%fsd(:,1)  = WG%PhiSD * 0.5
  met%fsd(:,2)  = WG%PhiSD * 0.5
  met%tk        = WG%Temp + 273.15
  met%ua        = WG%Wind
  met%coszen    = WG%coszen
  met%tvrad     = met%tk !!!! Exclamation mark!
  DO i = 1, PLUME%mland
     ! compute qv 
     CALL rh_sh ( PLUME%MET(rhum)%VAL(i), met%tk(i), met%pmb(i), met%qv(i) )
  END DO

  IF ( ktau.EQ.1 ) THEN
     WRITE(*,*)"Pmax,Tmin, max,Ntmin", PLUME%MET(PrevTmax)%VAL(19), &
          PLUME%MET(Tmin)%VAL(19),PLUME%MET(Tmax)%VAL(19),PLUME%MET(NextTmin)%VAL(19) 
     WRITE(*,*)"#  qv       Precip   snow   LWDin PhiLD  rPhiLD  PhiSD Temp    Wind    coszen"
    
     
  ENDIF
  
  IF ( ktau .lt. 9 ) &
  WRITE(*,FMT='(I2,3(X,F8.5),X,3(F5.1,x),F7.3,2(X,F5.1),3(X,ES12.4))')&
       MOD((ktau-1)*NINT(dt/3600.),24), met%qv(19), &
       WG%Precip(19), met%precip_sn(19), PLUME%MET(lwdn)%VAL(19), &
       met%fld(19), WG%PhiLD(19) , WG%PhiSD(19) , met%tk(19) , WG%Wind(19), &
       WG%coszen(19)

  ! Finally closing files when done

  IF ((ktau .EQ. kend .AND. FILE_SWITCH( PLUME, 'CLOSE' )) .OR. islast) THEN
     DO i=1, PLUME%NMET
        STATUS = NF90_CLOSE(PLUME%F_ID(i))
        CALL HANDLE_ERR(STATUS, "Closing PLUME file"//PLUME%MetFile(i))
     END DO
  END IF
  
  ! CALL1 is over...
  CALL1 = .FALSE.

END SUBROUTINE PLUME_MIP_GET_MET

END MODULE CABLE_PLUME_MIP
