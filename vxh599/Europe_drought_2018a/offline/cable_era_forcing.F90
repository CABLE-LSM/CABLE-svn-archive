MODULE CABLE_ERA

  USE netcdf                         ! Access to netcdf routines
  USE CABLE_COMMON_MODULE, ONLY: &   ! Selected cable_common.f90 routines:
      HANDLE_ERR,  &                 ! Print error status info returned by netcdf file operations
      GET_UNIT                       ! Finds an unused unit number for file opens

  USE cable_IO_vars_module, ONLY: &  ! Selected cable_iovars.F90 variables:
      logn,            &             ! Log file unit number
      land_x, land_y,  &             ! Col (x) & row (y) indices of each land point in land mask (dimension mland)
      exists                         ! Only for exists%Snowf, which we will set to .FALSE. because there is no snow 
                                     ! in CRU-NCEP. Setting this ensures snow will be determined in CABLE from temperature.

  IMPLICIT NONE

  ! Define a type for ERA information, and the subtype METVALS

  TYPE ERA_MET_TYPE 
    REAL, DIMENSION(:), ALLOCATABLE :: METVALS  ! Define a spatial vector of meteorology for one timestep
  END TYPE ERA_MET_TYPE

  TYPE ERA_TYPE
    INTEGER  :: mland                    ! Number of land cells
    INTEGER  :: NMET                     ! Number of met variable types (rain, lwdn etc)
    INTEGER  :: xdimsize, ydimsize       ! Landmask grid size dimensions (x=cols, y=rows)
    INTEGER  :: tdimsize                 ! Time dimension of metfiles (met data timesteps per annual file)
    INTEGER  :: CYEAR                    ! Current run year, same as CurYear, not necessarily the same as MetYear
    INTEGER  :: MetStart                 ! First year of met
    INTEGER  :: MetEnd                   ! Last year of met
    INTEGER  :: CTSTEP                   ! Current met data timestep (1 to tdimsize, i.e. 365 for CRU-NCEP annual daily files)
    INTEGER  :: DTsecs                   ! Model timestep in seconds, converted from namelist value in hours
    INTEGER  :: ktau                     ! Current model timestep, reset at the start of a new year of met
    INTEGER  :: F_ID                     ! NetCDF object id's for files (NetCDF bookkeeping stuff) 
	INTEGER, DIMENSION(8) :: V_ID  ! NetCDF object id's for variables (NetCDF bookkeeping stuff) 
    REAL, DIMENSION(:), ALLOCATABLE :: CO2VALS  ! Global annual CO2 values (dim is the number of years of data, or 1 if time-invariant)
    REAL, DIMENSION(:), ALLOCATABLE :: SeqYears ! Sequence of MetYears in the S0 runs
    !REAL, DIMENSION(:), ALLOCATABLE :: MetYear  ! MetYear used for runs
    INTEGER :: MetYear  ! MetYear used for runs
    LOGICAL  :: DirectRead     ! Flag to do with reading small numbers of points efficiently. Set true for small numbers of points
    LOGICAL  :: LeapYears      ! Flag for whether leaps years occur, required by CABLE. (If False -> no Feb 29th) 
    LOGICAL,DIMENSION(:,:),ALLOCATABLE :: LandMask ! Logical landmask, true for land, false for non-land
  !
    CHARACTER(len=30)  :: Run            ! Where run type is      : "S0_TRENDY", "S1_TRENDY", "S2_TRENDY"
    CHARACTER(len=15)  :: CO2            ! CO2 takes value        : "static1860", "1860_1900", "1901_2015" 
    CHARACTER(len=15)  :: Ndep           ! Ndep takes value        : "static1860", "1860_1900", "1901_2015"
    CHARACTER(len=15)  :: Forcing        ! Met Forcing takes value: "spinup",        "spinup", "1901_2015" 
  !         
    CHARACTER(len=200) :: BasePath       ! Full path for the location of data used for CRU runs "/x/y"
    CHARACTER(len=200) :: MetPath        ! Full path for the location of the met files "/x/y"
    CHARACTER(len=200) :: LandMaskFile   ! Land mask filename, without path
    CHARACTER(len=30) ,DIMENSION(8) :: VAR_NAME  ! Netcdf variable 'Name' for each type of met (dim=# of met vars). Note: Name, not 'Title'
    CHARACTER(len=200) :: MetFile   ! Met file names incl metpath, constructed in ERA_GET_FILENAME (dim=# of met vars)
    TYPE(ERA_MET_TYPE), DIMENSION(8) :: MET  ! Met data vectors (METVALS) for one timestep, dim=# of met vars + 2 for prev Tmax and next Tmin  
    REAL,  DIMENSION(:), ALLOCATABLE :: NdepVALS
    INTEGER :: NdepF_ID, NdepV_ID
    INTEGER  :: Ndep_CTSTEP   ! counter for Ndep in input file   
  END TYPE ERA_TYPE

  TYPE (ERA_TYPE):: ERA  ! Define the variable ERA, of type ERA_TYPE

! Define local parameter names representing the position of each met var within variable MET. 
! prevTmax and nextTmin are special cases of Tmax and Tmin that do not count as extra met variables per se.  
  INTEGER, PRIVATE, PARAMETER :: &
       rainf   =  1, &
       snowf   =  2, &
       lwdown  =  3, &
       swdown  =  4, &
       psurf   =  5, &
       qair    =  6, &
       tair    =  7, &
       wind    =  8


! Error status of various operations (mostly netcdf-related). Typically 0 means ok, > 0 means unexpected condition.
  INTEGER, PRIVATE :: ErrStatus

  REAL, PRIVATE, PARAMETER :: SecDay = 86400. ! Number of seconds in a day
  REAL, PRIVATE, PARAMETER :: Pi = 3.14159265 ! copied from cable_weathergenerator.f90
  
! Filename prefix expected in the names of met files. Used by CRU_GET_FILENAME to construct met file names.
!  CHARACTER(len=6), DIMENSION(9), PARAMETER, PRIVATE :: &
!     !  PREF = (/ "rain  ", "lwdown", "swdown", "press ", "qair  ", "tmax  ", "tmin  ", "uwind ", "vwind " /)
!
!         PREF = (/ "pre   ", "dlwrf ", "dswrf ", "pres  ", "spfh  ", "tmax  ", "tmin  ", "ugrd  ", "vgrd  " /)

CONTAINS

!**************************************************************************************************

  SUBROUTINE ERA_INIT( ERA )

! Initialise the contents of the CRU defined type collection, from the CRU namelist file
! and by obtaining dimensions from the landmask
! JK: same for ERA

!**************************************************************************************************

    USE cable_IO_vars_module, ONLY: &
         latitude, longitude, & ! (R) Lat and long of landcells only (?)
         nmetpatches,         & ! (I) Size of patch dimension in met file, if it exists
         mask,                & ! (I) Land/sea mask (1,0)
         metGrid,             & ! (C4) Either 'land' or 'mask' for whether the data are packed or not (?)
         sdoy, smoy, syear,   & ! (I) Start time day of year, month, year
         shod,                & ! (R) Start time hour of day
         xdimsize, ydimsize,  & ! (I) Size of grid dimensions
         lat_all, lon_all       ! (R) Grids with the lat or lon of each cell (i.e. repetition along rows/cols), for CABLE. 

    USE cable_def_types_mod,  ONLY: mland  ! (I) Number of land cells

    IMPLICIT NONE

    TYPE (ERA_TYPE):: ERA

    INTEGER              :: ErrStatus  ! Error status returned by nc routines (zero=ok, non-zero=error) 
    INTEGER              :: nmlunit    ! Unit number for reading namelist file
    INTEGER              :: FID        ! NetCDF id for the landmask file
    INTEGER              :: latID, lonID, timID  ! NetCDF ids for dimensions in the landmask file
    INTEGER              :: landID     ! NetCDF id for the landmask variable in the landmask file
    INTEGER              :: tdimsize   ! Time dimension in the met file, not used here apparently (??)
    INTEGER              :: landcnt    ! Manually incremented counter for the number of land cells
    INTEGER              :: xcol, yrow ! Column and row position in the data file grids
    INTEGER              :: imetvar    ! loop counter through met variables

    ! Temporary local names for ERA% variables as they are read from the namelist file. 
    ! Note that ERA%CO2 and ERA%Forcing are assigned based on the value of Run, not read as options from the namelist file.
    LOGICAL              :: DirectRead = .FALSE.
    CHARACTER(len=30)    :: Run 
    CHARACTER(len=200)   :: BasePath
    CHARACTER(len=200)   :: MetPath
    CHARACTER(len=200)   :: LandMaskFile
    REAL                 :: DThrs   ! CABLE timestep (hrs), converted immediately to integer seconds for ERA%DTsecs
    REAL,DIMENSION(:)     ,ALLOCATABLE :: ERA_lats, ERA_lons  ! Lat/long values for each grid rows/cols from landmask.
    INTEGER,DIMENSION(:,:),ALLOCATABLE :: landmask

    ! Flag for errors
    LOGICAL              :: ERR = .FALSE.

    NAMELIST /ERANML/ BasePath, MetPath, LandMaskFile, Run, DThrs, DirectRead

    ! Read ERA namelist settings
    CALL GET_UNIT(nmlunit)  ! CABLE routine finds spare unit number
    OPEN (nmlunit,FILE="era.nml",STATUS='OLD',ACTION='READ')
    READ (nmlunit,NML=ERANML)
    CLOSE(nmlunit)

    ! Assign namelist settings to corresponding ERA defined-type elements
    ERA%BasePath     = BasePath
    ERA%MetPath      = MetPath
    ERA%LandMaskFile = LandMaskFile
    ERA%Run          = Run
    ERA%DTsecs       = int(DThrs * 3600.)  ! in seconds
    ERA%DirectRead   = DirectRead

    ! Assign Forcing and CO2 labels based only on the value of ERA%Run
    SELECT CASE (TRIM(ERA%Run))
    CASE( "S0" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "static1979"
       ERA%Ndep    = "static1979"
       WRITE(*   ,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1979'"
       WRITE(logn,*)"Run = 'spinup': Therefore Forcing = 'spinup', CO2 = 'static1979'"
    CASE( "S2.1" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "1979_2018"
       ERA%Ndep    = "1979_2018"
       WRITE(*   ,*)"Run = 'S2.1': Therefore Forcing = 'S2.1', CO2 = '1979_2018'"
       WRITE(logn,*)"Run = 'S2.1': Therefore Forcing = 'S2.1', CO2 = '1979_2018'"
    CASE( "S2.2.0" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "1979_2018"
       ERA%Ndep    = "1979_2018"
       WRITE(*   ,*)"Run = 'S2.2.0': Therefore Forcing = 'S2.2.0', CO2 = '1979_2018'"
       WRITE(logn,*)"Run = 'S2.2.0': Therefore Forcing = 'S2.2.0', CO2 = '1979_2018'"
    CASE( "S2.3.0" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "1979_2018"
       ERA%Ndep    = "1979_2018"
       WRITE(*   ,*)"Run = 'S2.3.0': Therefore Forcing = 'S2.3.0', CO2 = '1979_2018'"
       WRITE(logn,*)"Run = 'S2.3.0': Therefore Forcing = 'S2.3.0', CO2 = '1979_2018'"
    CASE( "S2.2.1-4" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "1979_2018"
       ERA%Ndep    = "1979_2018"
       WRITE(*   ,*)"Run = 'S2.2.1-4': Therefore Forcing = 'S2.2.1-4', CO2 = '1979_2018'"
       WRITE(logn,*)"Run = 'S2.2.1-4': Therefore Forcing = 'S2.2.1-4', CO2 = '1979_2018'"
    CASE( "S2.3.1-4" ) 
       ERA%Forcing = "spinup"
       ERA%CO2     = "1979_2018"
       ERA%Ndep    = "1979_2018"
       WRITE(*   ,*)"Run = 'S2.3.1-4': Therefore Forcing = 'S2.3.1-4', CO2 = '1979_2018'"
       WRITE(logn,*)"Run = 'S2.3.1-4': Therefore Forcing = 'S2.3.1-4', CO2 = '1979_2018'"
    CASE  default 
       WRITE(*   ,*)"Wrong ERA%Run: ",ERA%Run
       WRITE(*   ,*)"Use one of: S0, S2.1, S2.2.0, S2.3.0, S2.2.1-4, or S2.3.1-4!"
       WRITE(logn,*)"Wrong ERA%Run: ",ERA%Run
       WRITE(logn,*)"Use one of: S0, S2.1, S2.2.0, S2.3.0, S2.2.1-4, or S2.3.1-4!"
       ERR = .TRUE.
    END SELECT
   
    ! Print settings
    WRITE(*   ,*)"========================================= ERA ============"
    WRITE(*   ,*)"ERA settings chosen:"
    WRITE(*   ,*)" BasePath: ",TRIM(ERA%BasePath)
    WRITE(*   ,*)" LandMask: ",TRIM(ERA%LandMaskFile)
    WRITE(*   ,*)" Run                : ",TRIM(ERA%Run)
    WRITE(*   ,*)" Forcing (assigned) : ",TRIM(ERA%Forcing)
    WRITE(*   ,*)" CO2     (assigned) : ",TRIM(ERA%CO2)
    WRITE(*   ,*)" Ndep     (assigned) : ",TRIM(ERA%Ndep)
    WRITE(*   ,*)" DT(secs): ",ERA%DTsecs
    WRITE(logn,*)"========================================= ERA ============"
    WRITE(logn,*)"ERA settings chosen:"
    WRITE(logn,*)" BasePath: ",TRIM(ERA%BasePath)
    WRITE(logn,*)" LandMask: ",TRIM(ERA%LandMaskFile)
    WRITE(logn,*)" Run                : ",TRIM(ERA%Run)
    WRITE(logn,*)" Forcing (assigned) : ",TRIM(ERA%Forcing)
    WRITE(logn,*)" CO2     (assigned) : ",TRIM(ERA%CO2)
    WRITE(logn,*)" Ndep     (assigned) : ",TRIM(ERA%Ndep)
    WRITE(logn,*)" DT(secs): ",ERA%DTsecs

    ! Error trap for bad namelist. 
    IF ( ERR ) THEN
       WRITE(logn,*)"Invalid settings in ERA_INIT"
       STOP "Invalid settings in ERA_INIT"
    ENDIF

    ! If this is a S0_TRENDY run look for met data in the spinup directory instead.
    !IF (TRIM(CRU%Run) .EQ. "S0_TRENDY") THEN
    !   CRU%MetPath = TRIM(CRU%MetPath)//"/spinup_data"
    !ENDIF

    ! Set variable names to their NetCDF 'Names' (i.e. not their 'Titles')
    ERA%NMET = 8
    ERA%VAR_NAME(rainf)   = "Rainf"
    ERA%VAR_NAME(snowf)   = "Snowf"
    ERA%VAR_NAME(lwdown)  = "LWdown"
    ERA%VAR_NAME(swdown)  = "SWdown"
    ERA%VAR_NAME(psurf)   = "PSurf"
    ERA%VAR_NAME(qair)    = "Qair"
    ERA%VAR_NAME(tair)    = "Tair"
    ERA%VAR_NAME(wind)    = "Wind"

	

    WRITE(*   ,*)"========================================= ERA ============"
    WRITE(logn,*)"========================================= ERA ============"

    ! Now read landmask file
    ! Landmask file into init! Get LAt, LON etc. from there
    ! LMFILE = TRIM(CRU%LandMaskFile)
    WRITE(*   ,*) 'Opening ERA landmask file: ',TRIM(LandMaskFile)
    WRITE(logn,*) 'Opening ERA landmask file: ',TRIM(LandMaskFile)

    ! Open the land mask file
    ErrStatus = NF90_OPEN(TRIM(LandMaskFile), NF90_NOWRITE, FID)
    CALL HANDLE_ERR(ErrStatus, "Opening ERA Land-mask file"//TRIM(LandMaskFile))

    ! Latitude: Get the dimension ID, find the size of the dimension, assign it to ERA structure.
    ErrStatus = NF90_INQ_DIMID(FID,'latitude',latID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lat'"//TRIM(LandMaskFile))
    ERA%ydimsize = ydimsize

    ! Collect the latitudes into ERA_lats
    ALLOCATE( ERA_lats ( ydimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'latitude',latID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'latitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,latID,ERA_lats)
    CALL HANDLE_ERR(ErrStatus, "Reading 'latitudes'"//TRIM(LandMaskFile))

    ! Longitude: Get the dimension ID, find the size of the dimension, assign it to ERA structure.
    ErrStatus = NF90_INQ_DIMID(FID,'longitude',lonID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lon'"//TRIM(LandMaskFile))
    ERA%xdimsize = xdimsize

    ! Collect the longitudes into ERA_lons
    ALLOCATE( ERA_lons ( xdimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'longitude',lonID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'longitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,lonID,ERA_lons)
    CALL HANDLE_ERR(ErrStatus, "Reading 'longitudes'"//TRIM(LandMaskFile))

    ! Allocate the landmask arrays for... 
    ALLOCATE( ERA%landmask ( xdimsize, ydimsize) )  ! Passing out to other ERA routines (logical)
    ALLOCATE( landmask ( xdimsize, ydimsize) )      ! Local use in this routine (integer)
    ALLOCATE ( mask( xdimsize, ydimsize) )          ! Use by CABLE

    ! Check that the land mask variable is called "land" in the land mask file,
    ! and read it into local variable landmask
    ErrStatus = NF90_INQ_VARID(FID,'lsm',landID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lsm' "//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,landID,landmask)
    CALL HANDLE_ERR(ErrStatus, "Reading 'lsm' "//TRIM(LandMaskFile))

    ! Convert the integer landmask into the logical ERA%landmask
    WHERE ( landmask .GT. 0 )
      ERA%landmask = .TRUE.
      mask           = 1
    ELSEWHERE
      ERA%landmask = .FALSE.
      mask           = 0
    END WHERE

    ! Count the number of land cells -> mland
    ERA%mland = COUNT(ERA%landmask)

    ! Allocate CABLE land-only vectors for lat/long and row/col values/indices.
    ALLOCATE( latitude(ERA%mland), longitude(ERA%mland) )
    ALLOCATE( land_y  (ERA%mland), land_x   (ERA%mland) )
   
    ! Allocate vectors for each of the different met quantities, including extra 
    ! prev/next temperatures for the Cesarracio temperature calculations in the
    ! weather generator.
    DO imetvar = 1, ERA%NMET
      ALLOCATE( ERA%MET(imetvar)%METVALS(ERA%mland) )
    END DO
    !ALLOCATE( CRU%MET(prevTmax)%METVALS(CRU%mland) )
    !ALLOCATE( CRU%MET(nextTmin)%METVALS(CRU%mland) )
    ! allocate array for Nitrogen deposition input data
    ALLOCATE( ERA%NdepVALS(ERA%mland) )

    ! Copy the col/row and lat/long positions of each land cell into the corresponding
    ! land only CABLE vectors. Q: We know mland at this point. Why not use landcnt to confirm
    ! the correct value of mland? 
    landcnt = 1
    DO yrow = 1, ydimsize
       DO xcol = 1, xdimsize
          IF ( .NOT. ERA%landmask(xcol,yrow) ) CYCLE   ! Go to next iteration if not a land cell

!          WRITE(6,FMT='(A15,I5,2(1X,F8.2),2(1x,I3))')"i, lo,la, xcol,yrow",landcnt,CRU_lons(xcol),CRU_lats(yrow),xcol, yrow

          land_x(landcnt)    = xcol
          land_y(landcnt)    = yrow
          longitude(landcnt) = ERA_lons(xcol)
          latitude(landcnt)  = ERA_lats(yrow)
          landcnt = landcnt + 1
       END DO
    END DO

    ! Set global CABLE variables
    metGrid     = "mask"
    ALLOCATE( mask(xdimsize, ydimsize) )
    mask        = landmask
    mland       = ERA%mland
    nmetpatches = 1
    ALLOCATE( lat_all(xdimsize, ydimsize), lon_all(xdimsize, ydimsize) )
    DO xcol = 1, xdimsize
       lat_all(xcol,:) = ERA_lats
    END DO
    DO yrow = 1, ydimsize
       lon_all(:,yrow) = ERA_lons
    END DO

    ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
    shod        = 0.
    sdoy        = 1
    smoy        = 1
    syear       = ERA%CYEAR
    
    ! Used to rescale the diurnal cycle from Swinbank calculation to match CRU-NCEP provided value.
    !ALLOCATE( CRU%AVG_LWDN(mland) )

    DEALLOCATE ( landmask, ERA_lats, ERA_lons )

    ErrStatus = NF90_CLOSE(FID)
    CALL HANDLE_ERR(ErrStatus, "Closing mask-file"//TRIM(LandMaskFile))

  END SUBROUTINE ERA_INIT

  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE GET_ERA_CO2( ERA, CO2air )

! Get CO2 values for use with a ERA run. Assign a static 1979 value if specified otherwise
! on the first call read all the annual values from a file into the ERA%CO2VALS array. On the first
! and subsequent   

  IMPLICIT NONE
  
  TYPE(ERA_TYPE) :: ERA           ! All the info needed for CRU met runs
  REAL, INTENT(OUT)    :: CO2air  ! A single annual value of CO2air in ppm for the current year.

  INTEGER              :: i, iunit, iyear, IOS = 0
  CHARACTER            :: CO2FILE*200
  LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine 

! For S0_TRENDY, use only static 1979 CO2 value and return immediately
  IF ( TRIM(ERA%CO2) .EQ. "static1979") THEN
    !CO2air = 286.42   ! CO2 in ppm for 1860
    CO2air = 335.82   ! CO2 in ppm for 1979, see ftp://ftp.lrz.de/transfer/ERA5Forcing/global_co2_ann_1700_2018.txt
    RETURN

! If not S0, varying CO2 values will be used...
  ELSE

! On the first call, allocate the ERA%CO2VALS array to store the entire history of annual CO2 
! values, open the (ascii) CO2 file and read the values into the array. 
    IF (CALL1) THEN
      ALLOCATE( ERA%CO2VALS( 1700:2017 ) )
      CO2FILE = TRIM(ERA%BasePath)//"/co2/global_co2_ann_1700_2018.txt"
      CALL GET_UNIT(iunit)
      OPEN (iunit, FILE=TRIM(CO2FILE), STATUS="OLD", ACTION="READ")
      DO WHILE( IOS .EQ. 0 )
        READ(iunit, FMT=*, IOSTAT=IOS) iyear, ERA%CO2VALS(iyear)
      END DO
      CLOSE(iunit)
      
      CALL1 = .FALSE.

    END IF

! In all varying CO2 cases, return the element of the array for the current year
! as a single CO2 value.
! 
    CO2air = ERA%CO2VALS( ERA%CYEAR ) 

  END IF

  END SUBROUTINE GET_ERA_CO2

  
  
!**************************************************************************************************

  SUBROUTINE GET_ERA_Ndep( ERA, month )

! Get Ndep values for use with a ERA run. Assign a static 1979 value if specified otherwise
! on the first call read all the annual values from a file into the ERA%CO2VALS array.  

  IMPLICIT NONE
  
  TYPE(ERA_TYPE), INTENT(INOUT) :: ERA           ! All the info needed for ERA met runs
  INTEGER          :: month
  REAL,ALLOCATABLE :: tmparr(:,:) 
  INTEGER          :: i, iunit, iyear, IOS = 0, k, t
  INTEGER          :: xds, yds        ! Ndep file dimensions of long (x), lat (y)
 
  LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine 
  CHARACTER(400) :: NdepFILE
  ! Abbreviate dimensions for readability.
  xds = ERA%xdimsize
  yds = ERA%ydimsize
  allocate(tmparr(xds,yds))
  ! For S0_TRENDY, use only static 1860 CO2 value and return immediately



  ! On the first call, allocate the ERA%CO2VALS array to store the entire history of annual CO2 
  ! values, open the (ascii) CO2 file and read the values into the array. 
  IF (CALL1) THEN

     ! CRU
     !NdepFILE = TRIM(ERA%BasePath)// &
     !     "/ndep/NOy_plus_NHx_dry_plus_wet_deposition_hist_1850_2015_annual_1deg.nc"  
     ! ERA
     NdepFILE = TRIM(ERA%BasePath)// &
           "/ndep/CCMI_ndep_sum_nhx_noy_1860-2016_0.5deg.nc"

     
     ! Open the NDep and access the variables by their name and variable id.
     WRITE(*   ,*) 'Opening ndep data file: ', NdepFILE
     WRITE(logn,*) 'Opening ndep data file: ', NdepFILE


     ErrStatus = NF90_OPEN(TRIM(NdepFILE), NF90_NOWRITE, ERA%NdepF_ID)  
     CALL HANDLE_ERR(ErrStatus, "Opening ERA file "//NdepFILE )
     ErrStatus = NF90_INQ_VARID(ERA%NdepF_ID,'N_deposition', ERA%NdepV_ID)
     CALL HANDLE_ERR(ErrStatus, "Inquiring ERA var "//"N_deposition"// &
          " in "//NdepFILE )

     ! Set internal counter
     ERA%Ndep_CTSTEP = 1

     IF ( TRIM(ERA%Ndep) .EQ. "static1979" .OR. ERA%CYEAR<=1979) THEN
       ! read Ndep at year 1979
        ERA%Ndep_CTSTEP = 1440 + month  ! timestep of control (210 * 12 + 1)
        t =  ERA%Ndep_CTSTEP
        ErrStatus = NF90_GET_VAR(ERA%NdepF_ID, ERA%NdepV_ID, tmparr, &
             start=(/1,1,t/),count=(/xds,yds,1/) )
        CALL HANDLE_ERR(ErrStatus, "Reading from "//NdepFILE )
        DO k = 1, ERA%mland
           ERA%NdepVALS(k) = tmparr( land_x(k), land_y(k) )
        END DO
      
       
     END IF
     CALL1 = .FALSE.
  END IF

  IF ( TRIM(ERA%Ndep) .NE. "static1979" .and.  ERA%CYEAR>1979) THEN
  
     ! read Ndep at current year (noting that file starts at 1850 and ends in 2015)
     ERA%Ndep_CTSTEP = (min(ERA%CYEAR, 2016) - 1860 + 1)*12 + month
     t =  ERA%Ndep_CTSTEP
     ErrStatus = NF90_GET_VAR(ERA%NdepF_ID, ERA%NdepV_ID, tmparr, &
          start=(/1,1,t/),count=(/xds,yds,1/) )
     CALL HANDLE_ERR(ErrStatus, "Reading from "//NdepFILE )
     DO k = 1, ERA%mland
        ERA%NdepVALS(k) = tmparr( land_x(k), land_y(k) )
     END DO
    
  END IF

  
END SUBROUTINE GET_ERA_Ndep


!**************************************************************************************************

  SUBROUTINE OPEN_ERA_MET( ERA )

! Opens each of the met files required for one year. This is where the distinction is made between
! the nominal run year (CYEAR) and the year of met required (MetYear), which is different for 
! S0_TRENDY and S1_TRENDY than for a standard run (S2_TRENDY). 

  USE cable_IO_vars_module, ONLY: timeunits ! (Char33) Name of time units read from nc file 

  IMPLICIT NONE

  TYPE( ERA_TYPE ), INTENT(INOUT) :: ERA ! All ERA related quantities and flags

  INTEGER             :: iunit, iyear, IOS = 0
  INTEGER             :: iVar            ! Loop counter through met variables
  INTEGER             :: tID             ! Numerical variable identifier returned by NetCDF routines,
                                         ! in this case for time. Needed to retrieve the time units.
  !INTEGER             :: MetYear         ! Year of met to access. Equals CYEAR for normal runs, but 
  !                                       ! must be calculated for S0_TRENDY and initialisation runs.
  INTEGER, SAVE       :: RunStartYear    ! The value of ERA%CYEAR on the first call, also equals syear.
                                         ! Allows the calculation of MetYear during S0_TRENDY and init runs. 
  LOGICAL, SAVE       :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine
  CHARACTER(300)      :: SEQFILE         ! Sequence of years for the S0 run
  CHARACTER(4)        :: CY              ! Character representation of ERA%MetYear
  REAL :: tmp

! Keep the initial value of CYEAR for calculation of different MetYear if required. 
  !IF (CALL1) RunStartYear = 1710 ! edit vh !
  !IF (CALL1) RunStartYear = 1691 ! edit vh !
  
  !DO iVar = 1, ERA%NMET  ! For each met variable (not needed for era forcing, variable loop starts at the end of the subroutine)

! Determine MetYear:
! For S0_TRENDY and initialisation, calculate the required met year for repeatedly cycling through the 
! 30 years of 1901-1930 spinup meteorology. For normal runs 1901-2015, MetYear = CYEAR.
     
  IF (TRIM(ERA%Run) .EQ. 'S0') THEN  
    IF (CALL1) THEN  ! if first call, read sequence of years

      ALLOCATE( ERA%SeqYears( 1970:2018 ) )
      SEQFILE = TRIM(ERA%BasePath)//"/S0_SequeceYears.txt"
	  
      CALL GET_UNIT(iunit)
      OPEN (iunit, FILE=TRIM(SEQFILE), STATUS="OLD", ACTION="READ")
      DO WHILE( IOS .EQ. 0 )
         READ(iunit, FMT=*, IOSTAT=IOS) iyear, tmp
         ERA%SeqYears(iyear) = tmp
      END DO
      CLOSE(iunit)
	  
      CALL1 = .FALSE.
	  
    END IF
    ERA%MetYear = INT(ERA%SeqYears( ERA%CYEAR ))
    
  ELSE
	
    ERA%MetYear = ERA%CYEAR
	
  END IF
	

  !CALL CRU_GET_FILENAME( CRU, MetYear, iVar, CRU%MetFile(iVar) ) ! Call routine to build the filenames.
  WRITE(CY,FMT='(I4)')ERA%MetYear
  
  ERA%MetFile = TRIM(ERA%MetPath)//"era5_europe_S2.1_"//CY//"_0.5deg.nc"
  
  ! Open the new met files and access the variables by their name and variable id.
  WRITE(*   ,*) 'Opening met data file: ', ERA%MetFile
  WRITE(logn,*) 'Opening met data file: ', ERA%MetFile

  ErrStatus = NF90_OPEN(TRIM(ERA%MetFile), NF90_NOWRITE, ERA%F_ID)  
  CALL HANDLE_ERR(ErrStatus, "Opening ERA file "//ERA%MetFile)
  
  DO iVar = 1, ERA%NMET
    
	ErrStatus = NF90_INQ_VARID(ERA%F_ID,TRIM(ERA%VAR_NAME(iVar)), ERA%V_ID(iVar))
    CALL HANDLE_ERR(ErrStatus, "Inquiring ERA var "//TRIM(ERA%VAR_NAME(iVar))// &
         " in "//ERA%MetFile )
  
  END DO

  ! Set internal counter
  ERA%CTSTEP = 1

  CALL1 = .FALSE. ! No longer the first call (saved).

  END SUBROUTINE OPEN_ERA_MET

!**************************************************************************************************

  SUBROUTINE ERA_GET_SUBDIURNAL_MET(ERA, MET, CurYear, ktau, kend, LastYearOfMet )

! Obtain one day of CRU-NCEP meteorology, subdiurnalise it using a weather
! and return the result to the CABLE driver.

  USE cable_def_types_mod,   ONLY: MET_TYPE
  USE cable_IO_vars_module,  ONLY: LANDPT, latitude
  USE cable_common_module,   ONLY: DOYSOD2YMDHMS
  USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
                                   WGEN_DAILY_CONSTANTS
  USE cable_checks_module,   ONLY: rh_sh
 

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: CurYear, ktau, kend
  LOGICAL, INTENT(IN)  :: LastYearOfMet

  TYPE(ERA_TYPE) :: ERA

! Define MET the CABLE version, different from the MET defined and used within the CRU variable. 
! The structure of MET_TYPE is defined in cable_define_types.F90 
  TYPE(MET_TYPE) :: MET   
 
! Local variables
  LOGICAL   :: newday, LastDayOfYear  ! Flags for occurence of a new day (0 hrs) and the last day of the year.
  INTEGER   :: iVar, ii, t, k, x, y, realk
  INTEGER   :: iland                  ! Loop counter through 'land' cells (cells in the spatial domain)
  INTEGER   :: itimestep              ! Loop counter through subdiurnal timesteps in a day
  INTEGER   :: imetvar                ! Loop counter through met variables
  INTEGER   :: dM, dD                 ! Met date as year, month, and day returned from DOYSOD2YMDHMS
  INTEGER   :: xds, yds               ! Met dimensions of long (x) and lat (y)
  INTEGER   :: is, ie                 ! Starting and ending vegetation type per land cell
  REAL      :: dt                     ! Timestep in seconds
  REAL      :: CO2air                 ! CO2 concentration in ppm
  REAL      :: etime
  CHARACTER :: LandMaskFile*200       ! Name of the land mask file
  ! variables needed for the calculation of coszen
  REAL      :: itime
  REAL      :: rntime
  REAL      :: ritime
  REAL      :: TimeNoon
  REAL      :: TimeRad
  
  TYPE(WEATHER_GENERATOR_TYPE), SAVE :: WG
  LOGICAL,    SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine
  REAL,ALLOCATABLE :: tmparr(:,:)  ! packing into ERA%MET(iVar)%METVALS(k)
  
! Purely for readability...
  dt = ERA%DTsecs

! On first step read and check CRU settings and read land-mask (only needed for coszen now)
  IF ( CALL1 ) CALL WGEN_INIT( WG, ERA%mland, latitude, dt )

! Pass time-step information to CRU 
  ERA%CYEAR = CurYear
  ERA%ktau  = ktau     ! ktau is the current timestep in the year.

!!!!  this only works with CANBERRA cable_driver, as ktau    !!!!
!!!!  restarts on Jan 1                                      !!!!
! Based on the ktau timestep, calculate date and time information (the same for the whole spatial dimension.)
  met%hod (:) = REAL(MOD( (ktau-1) * NINT(dt), INT(SecDay)) ) / 3600.  ! Hour of the day
  met%doy (:) = INT(REAL(ktau-1) * dt / SecDay ) + 1                   ! Day of Year = days since 0 hr 1st Jan 
  met%year(:) = CurYear                                                ! Current year

! Using the day-of-year and seconds-of-day calculate the month and day-of-month, using the time information
! for the first land cell only (because they will be the same across the domain).
!                            In      In         In        Out       Out       Optional Out 
! SUBROUTINE DOYSOD2YMDHMS( Year, Yearday, SecondsOfDay, Month, DayOfMonth, [Hour, Min, Sec])
 
  CALL DOYSOD2YMDHMS(CurYear, INT(met%doy(1)), INT(met%hod(1)) * 3600, dM, dD)

  met%moy (:) = dM     ! Record the month

! It's a new day if the hour of the day is zero. 
  newday = ( met%hod(landpt(1)%cstart).EQ. 0 )

! Beginning-of-year accounting
  IF (ktau .EQ. 1) THEN  ! ktau is always reset to 1 at the start of the year.

! Read a new annual CO2 value and convert it from ppm to mol/mol
    CALL GET_ERA_CO2( ERA, CO2air )
    met%ca(:) = CO2air / 1.e+6  ! 

    CALL GET_ERA_Ndep( ERA, dM )
    DO iland = 1, ERA%mland
       met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
            ERA%NdepVALS(iland)  ! kg/m2/s > g/m2/d (1000.*3600.*24.), JK: right units here already
    END DO
    
! Open a new annual ERA met file.
    CALL OPEN_ERA_MET( ERA )

  ENDIF

! %%%%%% PRB to add his own comments from here down for this routine.
! Now get the Met data for this day
  IF ( newday ) THEN
     

!print *, CRU%CTSTEP, ktau, kend
!   CALL CPU_TIME(etime)
!   PRINT *, 'b4 daily ', etime, ' seconds needed '

    LastDayOfYear = (ktau .EQ. kend-((SecDay/dt)-1))
    CALL WGEN_DAILY_CONSTANTS( WG, ERA%mland, INT(met%doy(1))+1 )

  ENDIF


  ! calculate coszen (code taken from cable_weathergenerator.f90)
  ritime = REAL(itime)     * dt/3600.  ! Convert the current time to real
  rntime = REAL(NINT(REAL(SecDay))/dt) * dt/3600.  ! Convert ntime to real

  TimeNoon = ritime/rntime - 0.5
  TimeRad  = 2.0*Pi*TimeNoon

  WG%coszen = MAX(0.0, ( SIN(WG%DecRad)*SIN(WG%LatRad) + COS(WG%DecRad)*COS(WG%LatRad)*COS(TimeRad) ))
  

  
  ! JK: Main routine from GET_ERA_DAILY_MET is called here directly 
  !CALL GET_ERA_DAILY_MET( ERA, LastDayOfYear, LastYearOfMet )
  xds = ERA%xdimsize
  yds = ERA%ydimsize
  allocate(tmparr(xds,yds))
  
  DO iVar = 1, ERA%NMET

    ! iVar is not Tmin or Tmax so the variable index and timestep index are unchanged.
    ii = iVar
    t  = ERA%CTSTEP

   

    ! Standard read of the current variable, for the current timestep:
    ! Directly read the current points into the met vector (more efficient for small domains), 
    ! or read the whole grid into tmparr and extract them from there. 
    IF ( ERA%DirectRead ) THEN

      DO k = 1, ERA%mland
        ErrStatus = NF90_GET_VAR(ERA%F_ID, ERA%V_ID(iVar), ERA%MET(ii)%METVALS(k), &
                    start=(/land_x(k),land_y(k),t/) )
        CALL HANDLE_ERR(ErrStatus, "Reading directly from "//ERA%MetFile)
      END DO

    ELSE

      ErrStatus = NF90_GET_VAR(ERA%F_ID, ERA%V_ID(iVar), tmparr, &
                  start=(/1,1,t/),count=(/xds,yds,1/) )
      CALL HANDLE_ERR(ErrStatus, "Reading from "//ERA%MetFile)
      DO k = 1, ERA%mland
        ERA%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
      END DO
  
   ENDIF

  END DO  ! End of loop through all met variables
  
  

  
  ! JK: for ERA forcing, only convert variables, no weather generator needed
  ! assign ERA variables to cable met variables
  DO iland = 1, ERA%mland
    is = landpt(iland)%cstart
    ie = landpt(iland)%cend   
	
    met%ua(is:ie)        = ERA%MET(wind)%METVALS(iland)
    met%precip(is:ie)    = (ERA%MET(rainf)%METVALS(iland) + ERA%MET(snowf)%METVALS(iland)) * ERA%DTsecs
    met%precip_sn(is:ie) = ERA%MET(snowf)%METVALS(iland)* ERA%DTsecs
    met%tk(is:ie)        = ERA%MET(tair)%METVALS(iland)
    met%fld(is:ie)       = ERA%MET(lwdown)%METVALS(iland)
    met%qv(is:ie)        = ERA%MET(qair)%METVALS(iland)
    met%pmb(is:ie)       = ERA%MET(psurf)%METVALS(iland) / 100.  ! Convert pressure Pa -> hPa (formerly in GET_ERA_DAILY_MET)
	  
    ! Cable's swdown is split into two components, visible and nir, which 
    ! get half of the CRU-NCEP swdown each.
    met%fsd(is:ie,1) = ERA%MET(swdown)%METVALS(iland) * 0.5  ! Visible 
    met%fsd(is:ie,2) = ERA%MET(swdown)%METVALS(iland) * 0.5  ! NIR


 END DO

 where (met%fsd .lt. 1.e-2)
    met%fsd = 0.0
 endwhere

 where (met%precip .lt. 0.0)
    met%precip = 0.0
 endwhere

 where (met%precip_sn .lt. 0.0)
    met%precip_sn = 0.0
 endwhere

  ! initialise within canopy air temp
  met%tvair = met%tk
  met%tvrad = met%tk

  ! If this is the end of the year or the end of the met, close the current met files.
  !print *, "ktau, kend, LastYearOfMet as close test:", ktau, kend
  IF (ktau .EQ. kend) THEN
    !DO imetvar=1, ERA%NMET
    !print *, 'Close CRU%MetFile(imetvar)', CRU%MetFile(imetvar)
    ErrStatus = NF90_CLOSE(ERA%F_ID)
    CALL HANDLE_ERR(ErrStatus, "Closing ERA file"//ERA%MetFile)
    !END DO
  END IF

 ! Increment the internal timestep counter
  ERA%CTSTEP = ERA%CTSTEP + 1

  ! CALL1 is over...
  CALL1 = .FALSE.

  END SUBROUTINE ERA_GET_SUBDIURNAL_MET

END MODULE CABLE_ERA
