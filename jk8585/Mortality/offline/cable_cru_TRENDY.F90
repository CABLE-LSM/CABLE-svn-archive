MODULE CABLE_CRU

  USE netcdf                        ! Access to netcdf routines
  USE CABLE_COMMON_MODULE, ONLY: &  ! Selected cable_common.f90 routines:
      HANDLE_ERR,  &                ! Print error status info returned by netcdf file operations
      GET_UNIT                      ! Finds an unused unit number for file opens
  USE cable_IO_vars_module, ONLY: & ! Selected cable_iovars.F90 variables:
      logn,            &            ! Log file unit number
      land_x, land_y,  &            ! Col (x) & row (y) indices of each land point in land mask (dimension mland)
      exists                        ! Only for exists%Snowf, which we will set to .FALSE. because there is no snow
                                    ! in CRU-NCEP. Setting this ensures snow will be determined in CABLE from temperature.

  IMPLICIT NONE

  integer, parameter :: sp = kind(1.0)

  ! Define a type for CRU-NCEP information, and the subtype METVALS

  TYPE CRU_MET_TYPE
    REAL, DIMENSION(:), ALLOCATABLE :: METVALS  ! Define a spatial vector of meteorology for one timestep
  END TYPE CRU_MET_TYPE

  TYPE CRU_TYPE
    INTEGER :: mland              ! Number of land cells
    INTEGER :: NMET               ! Number of met variable types (rain, lwdn etc) NOT INCLUDING prevTmax and nextTmin
    INTEGER :: xdimsize, ydimsize ! Landmask grid size dimensions (x=cols, y=rows)
    INTEGER :: tdimsize           ! Time dimension of metfiles (met data timesteps per annual file)
    INTEGER :: CYEAR              ! Current Run Year, Same As Curyear, Not Necessarily The Same As Metyear
    Integer :: Metstart           ! First Year Of Met
    Integer :: Metend             ! Last Year Of Met
    Integer :: Ctstep             ! Current Met Data Timestep (1 To Tdimsize, I.E. 365 For Cru-Ncep Annual Daily Files)
    Integer :: Dtsecs             ! Model Timestep In Seconds, Converted From Namelist Value In Hours
    Integer :: Ktau               ! Current model timestep, reset at the start of a new year of met
    INTEGER :: metrecyc=32        ! number of years for the met recycling
    INTEGER, DIMENSION(9) :: F_ID, V_ID ! NetCDF object id's for files and variables (NetCDF bookkeeping stuff)
    ! Avg of one day's diurnal cycle of lwdn calculated by Swinbank. AVG_LWDN
    ! is used to rescale the diurnal cycle to match the day's CRUNCEP lwdn. (dim=mland)
    REAL, DIMENSION(:), ALLOCATABLE :: AVG_LWDN
    ! Global annual CO2 values (dim is the number of years of data, or 1 if time-invariant)
    REAL, DIMENSION(:), ALLOCATABLE :: CO2VALS
    LOGICAL :: DirectRead ! Flag to do with reading small numbers of points efficiently. Set true for small numbers of points
    LOGICAL :: LeapYears  ! Flag for whether leaps years occur, required by CABLE. Always false for CRUNCEP (no Feb 29th)
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: LandMask ! Logical landmask, true for land, false for non-land
    !
    CHARACTER(len=30)  :: Run            ! Where run type is      : "spinup", "P0" or "P1"
    CHARACTER(len=15)  :: CO2            ! CO2 takes value        : "static1860", "1860_1900", "1901_2015"
    CHARACTER(len=15)  :: Ndep           ! Ndep takes value        : "static1860", "1860_1900", "1901_2015"
    CHARACTER(len=15)  :: Forcing        ! Met Forcing takes value: "spinup",        "spinup", "1901_2015"
    CHARACTER(len=200) :: BasePath       ! Full path for the location of data used for CRU runs "/x/y"
    CHARACTER(len=200) :: MetPath        ! Full path for the location of the met files "/x/y"
    CHARACTER(len=50)  :: MetVersion     ! Met Forcing Version (currently CRUJRA_YEAR and VERIFY_2021)
    CHARACTER(len=200) :: LandMaskFile   ! Land mask filename, without path
    ! Netcdf variable 'Name' for each type of met (dim=# of met vars). Note: Name, not 'Title'
    CHARACTER(len=30),  DIMENSION(9)  :: VAR_NAME
    ! Met file names incl metpath, constructed in CRU_GET_FILENAME (dim=# of met vars)
    CHARACTER(len=200), DIMENSION(9)  :: MetFile
    ! Met data vectors (METVALS) for one timestep, dim=# of met vars + 2 for prev Tmax and next Tmin
    TYPE(CRU_MET_TYPE), DIMENSION(11) :: MET
    REAL, DIMENSION(:), ALLOCATABLE :: NdepVALS
    INTEGER :: NdepF_ID, NdepV_ID
    INTEGER :: Ndep_CTSTEP   ! counter for Ndep in input file
  END TYPE CRU_TYPE

  ! TYPE(CRU_TYPE) :: CRU  ! Define the variable CRU, of type CRU_TYPE

  ! Define local parameter names representing the position of each met var within variable MET.
  ! prevTmax and nextTmin are special cases of Tmax and Tmin that do not count as extra met variables per se.
  INTEGER, PRIVATE, PARAMETER :: &
       rain     =  1, &
       snow     =  2, &
       lwdn     =  3, &
       swdn     =  4, &
       pres     =  5, &
       qair     =  6, &
       tair     =  7, &
       wind     =  8

  ! Error status of various operations (mostly netcdf-related). Typically 0 means ok, > 0 means unexpected condition.
  INTEGER, PRIVATE :: ErrStatus

  REAL, PRIVATE, PARAMETER :: SecDay = 86400. ! Number of seconds in a day

  ! ! Filename prefix expected in the names of met files. Used by CRU_GET_FILENAME to construct met file names.
  ! CHARACTER(len=6), DIMENSION(9), PARAMETER, PRIVATE :: &
  !      !  PREF = (/ "rain  ", "lwdown", "swdown", "press ", "qair  ", "tmax  ", "tmin  ", "uwind ", "vwind " /)
  !      PREF = (/ "pre   ", "dlwrf ", "dswrf ", "pres  ", "spfh  ", "tmax  ", "tmin  ", "ugrd  ", "vgrd  " /)

CONTAINS

  !**************************************************************************************************

  SUBROUTINE CRU_INIT( CRU )

    ! Initialise the contents of the CRU defined type collection, from the CRU namelist file
    ! and by obtaining dimensions from the landmask

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
#ifdef __MPI__
    use mpi,                  only: MPI_Abort
#endif

    IMPLICIT NONE

    TYPE(CRU_TYPE), intent(inout) :: CRU

    INTEGER :: ErrStatus  ! Error status returned by nc routines (zero=ok, non-zero=error)
    INTEGER :: nmlunit    ! Unit number for reading namelist file
    INTEGER :: FID        ! NetCDF id for the landmask file
    INTEGER :: latID, lonID  ! NetCDF ids for dimensions in the landmask file
    INTEGER :: landID     ! NetCDF id for the landmask variable in the landmask file
    INTEGER :: landcnt    ! Manually incremented counter for the number of land cells
    INTEGER :: xcol, yrow ! Column and row position in the data file grids
    INTEGER :: imetvar    ! loop counter through met variables

    ! Temporary local names for CRU% variables as they are read from the namelist file.
    ! Note that CRU%CO2 and CRU%Forcing are assigned based on the value of Run, not read as options from the namelist file.
    LOGICAL            :: DirectRead = .FALSE.
    CHARACTER(len=30)  :: Run
    CHARACTER(len=200) :: BasePath
    CHARACTER(len=200) :: MetPath
    CHARACTER(len=50)  :: MetVersion
    CHARACTER(len=200) :: LandMaskFile
    REAL               :: DThrs   ! CABLE timestep (hrs), converted immediately to integer seconds for CRU%DTsecs
    REAL,    DIMENSION(:),   ALLOCATABLE :: CRU_lats, CRU_lons  ! Lat/long values for each grid rows/cols from landmask.
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: landmask

    ! Flag for errors
    LOGICAL :: ERR = .FALSE.
#ifdef __MPI__
    integer :: ierr
#endif

    NAMELIST /CRUNML/ BasePath, MetPath, MetVersion, LandMaskFile, Run, DThrs, DirectRead

    ! Read CRU namelist settings
    CALL GET_UNIT(nmlunit)  ! CABLE routine finds spare unit number
    OPEN(nmlunit, FILE="cru.nml", STATUS='OLD', ACTION='READ')
    READ(nmlunit, NML=CRUNML)
    CLOSE(nmlunit)

    ! Assign namelist settings to corresponding CRU defined-type elements
    CRU%BasePath     = BasePath
    CRU%MetPath      = MetPath
    CRU%MetVersion   = trim(MetVersion)
    CRU%LandMaskFile = trim(LandMaskFile)
    CRU%Run          = Run
    CRU%DTsecs       = int(DThrs * 3600.)  ! in seconds
    CRU%DirectRead   = DirectRead

    ! Assign Forcing and CO2 labels based only on the value of CRU%Run
    SELECT CASE(TRIM(CRU%Run))
    CASE( "spinup" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1901"
       CRU%Ndep    = "static1901"
       write(*,'(a)') "Run = 'spinup': Therefore Forcing = 'spinup', CO2 and Ndep = 'static1901'"
       WRITE(logn,*)  "Run = 'spinup': Therefore Forcing = 'spinup', CO2 and Ndep = 'static1901'"
    CASE( "P0" )
       CRU%Forcing = "spinup"
       CRU%CO2     = "static1901"
       CRU%Ndep    = "static1901"
       write(*,'(a)') "Run = 'P0': Therefore Forcing = 'spinup', CO2 and Ndep = 'static1901'"
       WRITE(logn,*)  "Run = 'P0': Therefore Forcing = 'spinup', CO2 and Ndep = 'static1901'"
    CASE( "P1" )
       CRU%Forcing = "1901_2014"
       CRU%CO2     = "1901_2014"
       CRU%Ndep    = "static1901"
       write(*,'(a)') "Run = 'P1': Therefore Forcing = '1901_2014', CO2 = '1901_2014'"
       WRITE(logn,*)  "Run = 'P1': Therefore Forcing = '1901_2014', CO2 = '1901_2014'"
    END SELECT

    ! Print settings
    write(*,'(a)') "========================================= CRU ============"
    write(*,'(a)') "CRU settings chosen:"
    WRITE(*,'(a)') "  BasePath: "//trim(CRU%BasePath)
    write(*,'(a)') "  LandMask: "//trim(CRU%LandMaskFile)
    write(*,'(a)') "  Run               : "//trim(CRU%Run)
    write(*,'(a)') "  Forcing (assigned): "//trim(CRU%Forcing)
    write(*,'(a)') "  CO2     (assigned): "//trim(CRU%CO2)
    write(*,'(a)') "  Ndep    (assigned): "//trim(CRU%Ndep)
    write(*,'(a,i5)') "  DT(secs): ", CRU%DTsecs
    WRITE(logn,*) "========================================= CRU ============"
    WRITE(logn,*) "CRU settings chosen:"
    WRITE(logn,*) " BasePath: ",TRIM(CRU%BasePath)
    WRITE(logn,*) " LandMask: ",TRIM(CRU%LandMaskFile)
    WRITE(logn,*) " Run               : ",TRIM(CRU%Run)
    WRITE(logn,*) " Forcing (assigned): ",TRIM(CRU%Forcing)
    WRITE(logn,*) " CO2     (assigned): ",TRIM(CRU%CO2)
    WRITE(logn,*) " Ndep    (assigned): ",TRIM(CRU%Ndep)
    WRITE(logn,*) " DT(secs): ",CRU%DTsecs

    ! Error trap for bad namelist.
    IF ( ERR ) THEN
       WRITE(logn,*) "Invalid settings in CRU_INIT"
       write(*,*) "Invalid settings in CRU_INIT"
#ifdef __MPI__
       call MPI_Abort(0, 5, ierr) ! Do not know comm nor rank here
#else
       stop 5
#endif
    ENDIF

    ! ! If this is a S0_TRENDY run look for met data in the spinup directory instead.
    ! IF (TRIM(CRU%Run) .EQ. "S0_TRENDY") THEN
    !   CRU%MetPath = TRIM(CRU%MetPath)//"/spinup_data"
    ! ENDIF

    ! ! Set variable names to their NetCDF 'Names' (i.e. not their 'Titles')
    ! CRU%NMET = 9
    ! CRU%VAR_NAME(rain)  = "Total_Precipitation"
    ! CRU%VAR_NAME(lwdn)  = "Incoming_Long_Wave_Radiation"
    ! CRU%VAR_NAME(swdn)  = "Incoming_Short_Wave_Radiation"
    ! CRU%VAR_NAME(pres)  = "Pression"
    ! CRU%VAR_NAME(qair)  = "Air_Specific_Humidity"
    ! ! CRU%VAR_NAME(tmax)  = "maximum_6h_air_temperature"
    ! ! CRU%VAR_NAME(tmin)  = "minimum_6h_air_temperature"
    ! CRU%VAR_NAME(tmax)  = "maximum_air_temperature"
    ! CRU%VAR_NAME(tmin)  = "minimum_air_temperature"
    ! CRU%VAR_NAME(uwind) = "U_wind_component"
    ! CRU%VAR_NAME(vwind) = "V_wind_component"

    ! GSWP3
    CRU%NMET = 8
    CRU%VAR_NAME(rain)  = "Rainf"
    CRU%VAR_NAME(snow)  = "Snowf"
    CRU%VAR_NAME(lwdn)  = "LWdown"
    CRU%VAR_NAME(swdn)  = "SWdown"
    CRU%VAR_NAME(pres)  = "PSurf"
    CRU%VAR_NAME(qair)  = "Qair"
    CRU%VAR_NAME(tair)  = "Tair"
    CRU%VAR_NAME(wind)  = "Wind"

    CRU%Metstart = 1901

    write(*,'(a)') "========================================= CRU ============"
    WRITE(logn,*)  "========================================= CRU ============"

    ! Now read landmask file
    ! Landmask file into init! Get LAt, LON etc. from there
    ! LMFILE = TRIM(CRU%LandMaskFile)
    write(*,'(a)') 'Opening CRU landmask file: '//TRIM(LandMaskFile)
    WRITE(logn,*)  'Opening CRU landmask file: ',TRIM(LandMaskFile)

    ! Open the land mask file
    ErrStatus = NF90_OPEN(TRIM(LandMaskFile), NF90_NOWRITE, FID)
    CALL HANDLE_ERR(ErrStatus, "Opening CRU Land-mask file"//TRIM(LandMaskFile))

    ! Latitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID,'latitude',latID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,latID,len=ydimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lat'"//TRIM(LandMaskFile))
    CRU%ydimsize = ydimsize

    ! Collect the latitudes into CRU_lats
    ALLOCATE( CRU_lats ( ydimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'latitude',latID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'latitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,latID,CRU_lats)
    CALL HANDLE_ERR(ErrStatus, "Reading 'latitudes'"//TRIM(LandMaskFile))

    ! Longitude: Get the dimension ID, find the size of the dimension, assign it to CRU.
    ErrStatus = NF90_INQ_DIMID(FID,'longitude',lonID)
    ErrStatus = NF90_INQUIRE_DIMENSION(FID,lonID,len=xdimsize)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'lon'"//TRIM(LandMaskFile))
    CRU%xdimsize = xdimsize

    ! Collect the longitudes into CRU_lons
    ALLOCATE( CRU_lons ( xdimsize ) )
    ErrStatus = NF90_INQ_VARID(FID,'longitude',lonID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'longitudes'"//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,lonID,CRU_lons)
    CALL HANDLE_ERR(ErrStatus, "Reading 'longitudes'"//TRIM(LandMaskFile))
    where (CRU_lons.gt.180.0)
       CRU_lons = CRU_lons - 360.0
    end where

    ! Allocate the landmask arrays for...
    ALLOCATE( CRU%landmask ( xdimsize, ydimsize) )  ! Passing out to other CRU routines (logical)
    ALLOCATE( landmask ( xdimsize, ydimsize) )      ! Local use in this routine (integer)
    ALLOCATE ( mask( xdimsize, ydimsize) )          ! Use by CABLE

    ! Check that the land mask variable is called "land" in the land mask file,
    ! and read it into local variable landmask
    ErrStatus = NF90_INQ_VARID(FID,'land',landID)
    CALL HANDLE_ERR(ErrStatus, "Inquiring 'land' "//TRIM(LandMaskFile))
    ErrStatus = NF90_GET_VAR(FID,landID,landmask)
    CALL HANDLE_ERR(ErrStatus, "Reading 'land' "//TRIM(LandMaskFile))

    ! Convert the integer landmask into the logical CRU%landmask
    WHERE ( landmask .GT. 0 )
      CRU%landmask = .TRUE.
      mask           = 1
    ELSEWHERE
      CRU%landmask = .FALSE.
      mask           = 0
    END WHERE

    ! Count the number of land cells -> mland
    CRU%mland = COUNT(CRU%landmask)

    ! Allocate CABLE land-only vectors for lat/long and row/col values/indices.
    ALLOCATE( latitude(CRU%mland), longitude(CRU%mland) )
    ALLOCATE( land_y  (CRU%mland), land_x   (CRU%mland) )

    ! Allocate vectors for each of the different met quantities, including extra
    ! prev/next temperatures for the Cesarracio temperature calculations in the
    ! weather generator.
    DO imetvar = 1, CRU%NMET
      ALLOCATE( CRU%MET(imetvar)%METVALS(CRU%mland) )
    END DO

    ! allocate array for Nitrogen deposition input data
    ALLOCATE( CRU%NdepVALS(CRU%mland) )

    ! Copy the col/row and lat/long positions of each land cell into the corresponding
    ! land only CABLE vectors. Q: We know mland at this point. Why not use landcnt to confirm
    ! the correct value of mland?
    landcnt = 1
    DO yrow = 1, ydimsize
       DO xcol = 1, xdimsize
          IF ( .NOT. CRU%landmask(xcol,yrow) ) CYCLE   ! Go to next iteration if not a land cell
          ! WRITE(6,FMT='(A15,I5,2(1X,F8.2),2(1x,I3))') "i, lo,la, xcol,yrow", landcnt, &
          !    CRU_lons(xcol),CRU_lats(yrow),xcol, yrow
          ! C
          land_x(landcnt)    = xcol
          land_y(landcnt)    = yrow
          longitude(landcnt) = CRU_lons(xcol)
          latitude(landcnt)  = CRU_lats(yrow)
          landcnt = landcnt + 1
       END DO
    END DO

    ! Set global CABLE variables
    metGrid     = "mask"
    ALLOCATE( mask(xdimsize, ydimsize) )
    mask        = landmask
    mland       = CRU%mland
    nmetpatches = 1
    ALLOCATE( lat_all(xdimsize, ydimsize), lon_all(xdimsize, ydimsize) )
    DO xcol = 1, xdimsize
       lat_all(xcol,:) = CRU_lats
    END DO
    DO yrow = 1, ydimsize
       lon_all(:,yrow) = CRU_lons
    END DO

    ! CABLE TIME-UNITS needed by load-parameters (only on CABLE_init)
    shod        = 0.
    sdoy        = 1
    smoy        = 1
    syear       = CRU%CYEAR

    ! Used to rescale the diurnal cycle from Swinbank calculation to match CRU-NCEP provided value.
    ALLOCATE( CRU%AVG_LWDN(mland) )

    DEALLOCATE ( landmask, CRU_lats, CRU_lons )

    ErrStatus = NF90_CLOSE(FID)
    FID = -1
    CALL HANDLE_ERR(ErrStatus, "Closing mask-file"//TRIM(LandMaskFile))

    ! set units to -1
    CRU%f_id = -1
    CRU%Ndepf_id = -1

  END SUBROUTINE CRU_INIT

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE CRU_GET_FILENAME(CRU, cyear, par, fn)

    ! Build the filename FN: One annual file of daily met for one met quantity.

    IMPLICIT NONE

    type(cru_type),     intent(in)  :: CRU    ! Information about CRU
    integer,            intent(in)  :: cyear  ! Current year as an integer
    integer,            intent(in)  :: par    ! Index (1-9) of which met quantity will be sought
    character(len=200), intent(out) :: fn     ! Met filename (outgoing)

    character(len=4)   :: cy     ! Character representation of cyear
    character(len=200) :: metp   ! Local repr of met path
    character(len=50)  :: cruver ! cru version as in filename

    ! Create a character version of the year for building that part of the filename.
    write(cy, fmt='(i4)') cyear

    ! Initialise the filename with the met path
    metp   = trim(CRU%MetPath)
    fn     = trim(metp)
    cruver = trim(CRU%MetVersion)

    ! Build the rest of the filename according to the value of par, which references 11 possible
    ! types of met through the parameter names rain, lwdn, etc.

    IF (trim(CRU%MetVersion) == "GSWP3") THEN
       SELECT CASE(par)
          CASE(rain) ; FN = TRIM(FN)//"/Rainf_"//trim(cruver)//"_"//cy//".nc"
          CASE(snow) ; FN = TRIM(FN)//"/Snowf_"//trim(cruver)//"_"//cy//".nc"
          CASE(lwdn) ; FN = TRIM(FN)//"/LWdown_"//trim(cruver)//"_"//cy//".nc"
          CASE(swdn) ; FN = TRIM(FN)//"/SWdown_"//trim(cruver)//"_"//cy//".nc"
          CASE(pres) ; FN = TRIM(FN)//"/PSurf_"//trim(cruver)//"_"//cy//".nc"
          CASE(qair) ; FN = TRIM(FN)//"/Qair_"//trim(cruver)//"_"//cy//".nc"
          CASE(tair) ; FN = TRIM(FN)//"/Tair_"//trim(cruver)//"_"//cy//".nc"
          CASE(wind) ; FN = TRIM(FN)//"/Wind_"//trim(cruver)//"_"//cy//".nc"
       END SELECT
   ENDIF

  END SUBROUTINE CRU_GET_FILENAME

  !**************************************************************************************************

  SUBROUTINE GET_CRU_CO2( CRU, CO2air )

    ! Get CO2 values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    IMPLICIT NONE

    TYPE(CRU_TYPE), intent(inout) :: CRU    ! All the info needed for CRU met runs
    REAL,           INTENT(OUT)   :: CO2air ! A single annual value of CO2air in ppm for the current year.

    INTEGER              :: iunit, iyear, IOS = 0
    CHARACTER(len=200)   :: CO2FILE
    LOGICAL,        SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

    ! For S0_TRENDY, use only static 1860 CO2 value and return immediately
    IF ( TRIM(CRU%CO2) .EQ. "static1860") THEN
       !CO2air = 286.42   ! CO2 in ppm for 1860
       CO2air = 276.59   ! CO2 in ppm for 1700
       RETURN
    ELSE IF ( TRIM(CRU%CO2) .EQ. "static1901") THEN
       CO2air = 296.80   ! CO2 in ppm for 1901
       RETURN

    ELSE ! If not S0_TRENDY, varying CO2 values will be used...

       ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
       ! values, open the (ascii) CO2 file and read the values into the array.
       IF (CALL1) THEN
          SELECT CASE (trim(CRU%MetVersion))
             CASE("GSWP3")
                ALLOCATE(CRU%CO2VALS(1700:2020))
                CO2FILE = TRIM(CRU%BasePath)//"/co2/global_co2_ann_1700_2020.txt"
             CASE("CRUJRA_2019")
                ALLOCATE(CRU%CO2VALS(1700:2018))
                CO2FILE = TRIM(CRU%BasePath)//"/co2/global_co2_ann_1700_2018.csv"
          END SELECT

          CALL GET_UNIT(iunit)
          OPEN(iunit, FILE=TRIM(CO2FILE), STATUS="OLD", ACTION="READ")
          DO WHILE( IOS .EQ. 0 )
             READ(iunit, FMT=*, IOSTAT=IOS) iyear, CRU%CO2VALS(iyear)
          END DO
          CLOSE(iunit)

          CALL1 = .FALSE.

       END IF

       ! In all varying CO2 cases, return the element of the array for the current year
       ! as a single CO2 value.
       !
       CO2air = CRU%CO2VALS( CRU%CYEAR )

    END IF

  END SUBROUTINE GET_CRU_CO2

  ! **************************************************************************************************

  SUBROUTINE GET_CRU_Ndep( CRU )

    ! Get Ndep values for use with a CRU-NCEP run. Assign a static 1860 value if specified otherwise
    ! on the first call read all the annual values from a file into the CRU%CO2VALS array. On the first
    ! and subsequent

    IMPLICIT NONE

    TYPE(CRU_TYPE), INTENT(INOUT) :: CRU ! All the info needed for CRU met runs

    INTEGER :: k
    LOGICAL, SAVE  :: CALL1 = .TRUE. ! A *local* variable recording the first call of this routine


    ! On the first call, allocate the CRU%CO2VALS array to store the entire history of annual CO2
    ! values, open the (ascii) CO2 file and read the values into the array.
    IF (CALL1) THEN

       IF ( TRIM(CRU%Ndep) .EQ. "static1901" .OR. CRU%CYEAR<=1850) THEN
          ! set Ndep to 2kg/ha/year 
          DO k = 1, CRU%mland
             CRU%NdepVALS(k) = 2.0
          END DO
       END IF
       CALL1 = .FALSE.
    END IF

  END SUBROUTINE GET_CRU_Ndep

  !**************************************************************************************************

  SUBROUTINE OPEN_CRU_MET( CRU )

    ! Opens each of the met files required for one year. This is where the distinction is made between
    ! the nominal run year (CYEAR) and the year of met required (MetYear), which is different for
    ! S0_TRENDY and S1_TRENDY than for a standard run (S2_TRENDY).

  IMPLICIT NONE

  TYPE(CRU_TYPE), INTENT(INOUT) :: CRU ! All CRU-NCEP related quantities and flags

  INTEGER       :: iVar            ! Loop counter through met variables
  INTEGER       :: MetYear         ! Year of met to access. Equals CYEAR for normal runs, but
                                   ! must be calculated for S0_TRENDY and initialisation runs.
  INTEGER, SAVE :: RunStartYear    ! The value of CRU%CYEAR on the first call, also equals syear.
                                   ! Allows the calculation of MetYear during S0_TRENDY and init runs.
  LOGICAL, SAVE :: CALL1 = .TRUE.  ! A *local* variable recording the first call of this routine

  ! Keep the initial value of CYEAR for calculation of different MetYear if required.
  !IF (CALL1) RunStartYear = 1710 ! edit vh !
  !IF (CALL1) RunStartYear = 1691 ! edit vh !
  IF (CALL1) THEN
       RunStartYear = 41 
  ENDIF

  DO iVar = 1, CRU%NMET  ! For each met variable

     ! For S0_TRENDY and initialisation, calculate the required met year for repeatedly cycling through the
     ! 30 years of 1901-1930 spinup meteorology. For normal runs 1901-2015, MetYear = CYEAR.
     ! JK: according to v9 protocol, met should be recycled from 1901-1920 (i.e. 20 years)
     !!$    IF ( TRIM(CRU%Run) .EQ. 'S0_TRENDY' .OR.  ( TRIM(CRU%Run) .EQ. 'S1_TRENDY' )) THEN
     !!$      MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear,30)
     !!$    ELSE IF ( TRIM(CRU%Run) .EQ. 'S2_TRENDY' ) THEN
     !!$      MetYear = CRU%CYEAR
     !!$    ENDIF
     IF (TRIM(CRU%Run) .EQ. 'spinup' .OR. TRIM(CRU%Run) .EQ. 'P0') THEN

        MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear, CRU%metrecyc)

     ELSE IF (TRIM(CRU%Run) .EQ. 'P1') THEN

        MetYear = CRU%CYEAR
     
     ENDIF

    CALL CRU_GET_FILENAME( CRU, MetYear, iVar, CRU%MetFile(iVar) ) ! Call routine to build the filenames.

    ! Open the new met files and access the variables by their name and variable id.
    write(*,'(a)') 'Opening met data file: '//trim(CRU%MetFile(iVar))
    WRITE(logn,*)  'Opening met data file: ', CRU%MetFile(iVar)

    ErrStatus = NF90_OPEN(TRIM(CRU%MetFile(iVar)), NF90_NOWRITE, CRU%F_ID(iVar))
    CALL HANDLE_ERR(ErrStatus, "Opening CRU file "//trim(CRU%MetFile(iVar)) )
    ErrStatus = NF90_INQ_VARID(CRU%F_ID(iVar),TRIM(CRU%VAR_NAME(iVar)), CRU%V_ID(iVar))
    CALL HANDLE_ERR(ErrStatus, "Inquiring CRU var "//TRIM(CRU%VAR_NAME(iVar))// &
         " in "//trim(CRU%MetFile(iVar)) )
  END DO

  ! Set internal counter
  CRU%CTSTEP = 1

  CALL1 = .FALSE. ! No longer the first call (saved).

  END SUBROUTINE OPEN_CRU_MET

  !**************************************************************************************************

  SUBROUTINE CRU_GET_MET( CRU, MET, Curyear, ktau, kend )

    USE cable_def_types_mod,   ONLY: MET_TYPE, r_2
    USE cable_IO_vars_module,  ONLY: LANDPT, latitude
    USE cable_common_module,   ONLY: DOYSOD2YMDHMS
    USE cable_weathergenerator,ONLY: WEATHER_GENERATOR_TYPE, WGEN_INIT, &
         WGEN_DAILY_CONSTANTS, WGEN_SUBDIURNAL_MET
    use mo_utils,              only: eq
    
    IMPLICIT NONE

    TYPE(CRU_TYPE), intent(inout) :: CRU
    TYPE(MET_TYPE), intent(inout) :: MET
    INTEGER,        INTENT(IN)    :: CurYear, ktau, kend
    !LOGICAL,        INTENT(IN)    :: LastYearOfMet

    REAL    :: tmp
    INTEGER :: iVar, ii, k
    LOGICAL :: newday                 ! flag for newday (0 hrs)
    INTEGER :: iland                  ! Loop counter through 'land' cells
    INTEGER :: itimestep              ! Loop counter through subdiurnal timesteps in a day
    INTEGER :: is, ie                 ! Starting and ending vegetation type per land cell
    REAL    :: dt                     ! Timestep in seconds
    REAL    :: CO2air                 ! CO2 concentration in ppm
    INTEGER :: t, tplus1              ! The current and next timestep
    INTEGER :: fid, vid, tid          ! Netcdf id's for file, variable, and time
    INTEGER :: xds, yds, tds          ! Metfile dimensions of long (x), lat (y), and time (t)
    INTEGER :: MetYear                ! Year of meteorology currently in use
    INTEGER :: NextMetYear            ! Next met year: Where to look for the nextTmin on Dec 31st
    INTEGER :: dM, dD                 ! Met date as year, month, and day returned from DOYSOD2YMDHMS
    CHARACTER(LEN=200) :: filename

    TYPE(WEATHER_GENERATOR_TYPE), SAVE :: WG
    INTEGER, SAVE     :: RunStartYear   ! The value of CRU%CYEAR on the first call, also equals syear.
    ! Allows the calculation of MetYear during S0_TRENDY and init runs.
    LOGICAL, SAVE     :: CALL1 = .TRUE. ! A *local* variable recording the first call of this routine
    REAL, ALLOCATABLE :: tmparr(:)    ! packing into CRU%MET(iVar)%METVALS(k)
    REAL, ALLOCATABLE :: ea(:)          ! water vapour pressure in mb

    ! If first call...
    ! Keep the initial value of CYEAR for calculation of different MetYear if required.
    IF (CALL1) THEN
       !RunStartYear = CRU%CYEAR
       RunStartYear = 41
    ENDIF

    allocate(ea(CRU%mland))
    allocate(tmparr(CRU%mland))

    !xds = CRU%xdimsize
    !yds = CRU%ydimsize
    !allocate(tmparr(xds,yds))

    ! Abbreviate dimensions for readability.
    !xds = CRU%xdimsize
    !yds = CRU%ydimsize

    ! Purely for readability...
    dt = CRU%DTsecs


    
    ! On first step read and check CRU settings and read land-mask
    IF ( CALL1 ) CALL WGEN_INIT( WG, CRU%mland, latitude, dt )
    
    CRU%CYEAR = CurYear
    CRU%ktau  = ktau ! ktau is the current timestep in the year.

    ! Based on the ktau timestep, calculate date and time information (the same for the whole spatial dimension.)
    met%hod (:) = REAL(MOD( (ktau-1) * NINT(dt), INT(SecDay)) ) / 3600.  ! Hour of the day
    met%doy (:) = INT(REAL(ktau-1) * dt / SecDay ) + 1                   ! Day of Year = days since 0 hr 1st Jan
    met%year(:) = CurYear                                                ! Current year

    CALL DOYSOD2YMDHMS(CurYear, INT(met%doy(1)), INT(met%hod(1)) * 3600, dM, dD)

    met%moy (:) = dM     ! Record the month

    ! It's a new day if the hour of the day is zero.
    newday = eq(met%hod(landpt(1)%cstart), 0.0)


    write(86,*) "Curyear:", Curyear
    write(87,*) "ktau:", ktau
    ! Beginning-of-year accounting
    IF (ktau .EQ. 1) THEN  ! ktau is always reset to 1 at the start of the year.

       ! Read a new annual CO2 value and convert it from ppm to mol/mol
       CALL GET_CRU_CO2( CRU, CO2air )
       met%ca(:) = CO2air / 1.e+6  !

       CALL GET_CRU_Ndep( CRU )
       DO iland = 1, CRU%mland
          !met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
          !     CRU%NdepVALS(iland)*86400000.  ! kg/m2/s > g/m2/d (1000.*3600.*24.)
          met%Ndep(landpt(iland)%cstart:landpt(iland)%cend) = &
               CRU%NdepVALS(iland) * (1000.0 / 10000.0 / 365.0)  ! kgN ha-1 year-1 -> gN m-2 d-1
       END DO
       ! Open a new annual CRU-NCEP met file.
       CALL OPEN_CRU_MET( CRU )

    ENDIF



    ! Beginning of day
     IF ( newday ) THEN

        !!! The following lines are pseudo-code that are only needed to calculate coszen within
        !!! the weather generator!
        ! Convert wind from u and v components to wind speed by Pythagorean Theorem
	WG%WindDay        = real(CRU%MET(wind)%METVALS, r_2)
	! Convert all temperatures from K to C
	WG%TempMinDay     = real(CRU%MET(tair)%METVALS - 273.15, r_2)
        WG%TempMaxDay     = real(CRU%MET(tair)%METVALS - 268.15, r_2)
        WG%TempMinDayNext = real(CRU%MET(tair)%METVALS - 273.15, r_2)
        WG%TempMaxDayPrev = real(CRU%MET(tair)%METVALS - 268.15, r_2)
        ! Convert solar radiation from J /m2/s to MJ/m2/d
	WG%SolarMJDay     = real(CRU%MET(swdn)%METVALS * 1.e-6 * SecDay, r_2) ! ->[MJ/m2/d]
        ! Convert precip from mm to m/day
	WG%PrecipDay      = real(max(CRU%MET(rain)%METVALS  * 1000., 0.0), r_2) ! ->[m/d]
        !WG%PrecipDay      = max(CRU%MET(  rain  )%METVALS  / 1000., 0.0)/2.0 ! ->[m/d] ! test vh !
        WG%SnowDay        = 0.0_r_2
        WG%VapPmbDay = real(esatf( real(WG%TempMinDay,sp) ), r_2)


        CALL WGEN_DAILY_CONSTANTS( WG, CRU%mland, INT(met%doy(1))+1 )

        DO itimestep = 1, NINT(SecDay/dt)
	     CALL WGEN_SUBDIURNAL_MET( WG, CRU%mland, itimestep-1 )
	    ! CRU%AVG_LWDN = CRU%AVG_LWDN + real(WG%PhiLD)
	END DO
	
	! To get the diurnal cycle for lwdn get whole day and scale with
	! LWDN from file later
     ENDIF

     CALL WGEN_SUBDIURNAL_MET( WG, CRU%mland, NINT(met%hod(1)*3600./dt) )
    
     ! Loop through all 9 met variables
     DO iVar = 1, CRU%NMET
        IF (TRIM(CRU%Run) .EQ. 'spinup' .OR. TRIM(CRU%Run) .EQ. 'P0') THEN

             MetYear = 1901 + MOD(CRU%CYEAR-RunStartYear, CRU%metrecyc)

        ELSE IF (TRIM(CRU%Run) .EQ. 'P1') THEN

             MetYear = CRU%CYEAR
     
        ENDIF

      ! iVar is not Tmin or Tmax so the variable index and timestep index are unchanged.
      ii = iVar
      t  = CRU%CTSTEP

      ! Standard read of the current variable, for the current timestep:
      ! Directly read the current points into the met vector (more efficient for small domains),
      ! or read the whole grid into tmparr and extract them from there.
      IF ( CRU%DirectRead ) THEN

         DO k = 1, CRU%mland
                ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), CRU%MET(ii)%METVALS(k), &
                            start=(/1,1,t/) )
                CALL HANDLE_ERR(ErrStatus, "Reading directly from "//trim(CRU%MetFile(iVar)))
         END DO

      ELSE
         
          ! ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
          !      start=(/1,1,t/),count=(/xds,yds,1/) )
          !ErrStatus = NF90_GET_VAR(CRU%F_ID(iVar), CRU%V_ID(iVar), tmparr, &
          !          start=t,count=1)
          CALL HANDLE_ERR(ErrStatus, "Reading from "//trim(CRU%MetFile(iVar)))
          DO k = 1, CRU%mland
             !CRU%MET(ii)%METVALS(k) = tmparr( land_x(k), land_y(k) )
             CRU%MET(ii)%METVALS(k) = tmparr(1)
          END DO
      ENDIF

    END DO  ! End of loop through all met variables

    ! Convert pressure Pa -> hPa

    CRU%MET(pres)%METVALS(:) = CRU%MET(pres)%METVALS(:) / 100.

    !! write met variables into the met structure and convert units if necessary
    DO iland = 1, CRU%mland

       is = landpt(iland)%cstart
       ie = landpt(iland)%cend

       met%pmb(is:ie)       = CRU%MET(pres)%METVALS(iland)
       met%precip(is:ie)    = real(CRU%MET(rain)%METVALS(iland) * dt)  ! precip in mm/time step
       met%precip_sn(is:ie) = real(CRU%MET(snow)%METVALS(iland) * dt)  ! precip in mm/time step
       met%tk(is:ie)        = CRU%MET(tair)%METVALS(iland)
       met%ua(is:ie)        = CRU%MET(wind)%METVALS(iland)
       met%fld(is:ie)       = CRU%MET(lwdn)%METVALS(iland)
       met%qv(is:ie)        = CRU%MET(qair)%METVALS(iland)

       ! Cable's swdown is split into two components, visible and nir, which
       ! get half of the CRU-NCEP swdown each.
       met%fsd(is:ie,1) = real(CRU%MET(swdn)%METVALS(iland) * 0.5_r_2)  ! Visible
       met%fsd(is:ie,2) = real(CRU%MET(swdn)%METVALS(iland) * 0.5_r_2)  ! NIR

       met%u10(is:ie)   = met%ua(is:ie)

       ! initialise within canopy air temp
       met%tvair(is:ie) = met%tk(is:ie)
       met%tvrad(is:ie) = met%tk(is:ie)
       met%pdep = 0.0

       ! At 3 hourly resolution, only coszen comes from the weather generator 
       met%coszen(is:ie) = real(WG%coszen(iland))

       ! Quick and dirty calculation of relative humidity (from bigleaf R package)
       ea(is:ie)        = met%qv(is:ie) * (met%pmb(is:ie) / 10.0) / ((1 - 0.622) * met%qv(is:ie) + 0.622)  ! kPa
       ea(is:ie)        = ea(is:ie) * 10.0 ! -> mbar (1kPa = 10 mb) 
       met%rhum(is:ie)  = real(ea(is:ie)/esatf(real(met%tk(is:ie),sp))) * 100.0 ! rel humidity (%)

    END DO


    ! Increment the internal timestep counter
    CRU%CTSTEP = CRU%CTSTEP + 1

    ! If this is the end of the year or the end of the met, close the current met files
    !print *, "ktau, kend, LastYearOfMet as close test:", ktau, kend
    IF (ktau .EQ. kend) THEN
       DO iVar=1, CRU%NMET
          ErrStatus = NF90_CLOSE(CRU%F_ID(iVar))
          CRU%F_ID(iVar) = -1
          CALL HANDLE_ERR(ErrStatus, "Closing CRU file"//trim(CRU%MetFile(iVar)))
       END DO
    END IF

    
    ! CALL1 can only happen once!
    IF (CALL1) CALL1 = .FALSE.

    !  print *, 'CRU%MET(tmin)%METVALS(:)'
    !  print *, CRU%MET(tmin)%METVALS(:)
    !  print *, 'CRU%MET(nexttmin)%METVALS(:)'
    !  print *, CRU%MET(nexttmin)%METVALS(:)
    !  print *, 'CRU%MET(prevtmax)%METVALS(:)'
    !  print *, CRU%MET(prevtmax)%METVALS(:)
    !  print *, 'CRU%MET(tmax)%METVALS(:)'
    !  print *, CRU%MET(tmax)%METVALS(:)
    !  print *, 'CRU%MET(pres)%METVALS(:)'
    !  print *, CRU%MET(pres)%METVALS(:)
    !  print *, 'CRU%MET(rain)%METVALS(:)'
    !  print *, CRU%MET(rain)%METVALS(:)

      CONTAINS

    ELEMENTAL FUNCTION Esatf(TC)
      ! ------------------------------------------------------------------------------
      ! At temperature TC [deg C], return saturation water vapour pressure Esatf [mb]
      ! from Teten formula.
      ! MRR, xx/1987
      ! PRB, 09/1999:   Convert to F95 elemental function; works on scalars and arrays
      !                 just like intrinsic functions.
      ! MRR, 12-mar-02: Convert Qsatf (specific humidity routine) to Esatf
      ! ------------------------------------------------------------------------------
      implicit none

      real(sp), intent(in) :: TC          ! temp [deg C]
      real(sp)             :: Esatf       ! saturation vapour pressure [mb]

      real(sp) :: TCtmp                   ! local
      real(sp),parameter:: A = 6.106      ! Teten coefficients
      real(sp),parameter:: B = 17.27      ! Teten coefficients
      real(sp),parameter:: C = 237.3      ! Teten coefficients

      TCtmp = TC                          ! preserve TC
      if (TCtmp .gt. 100.0) TCtmp = 100.0   ! constrain TC to (-40.0,100.0)
      if (TCtmp .lt. -40.0) TCtmp = -40.0

      Esatf = A*EXP(B*TCtmp/(C+TCtmp))    ! sat vapour pressure (mb)

    END FUNCTION Esatf

  END SUBROUTINE CRU_GET_MET

  !**************************************************************************************************

  subroutine cru_close(CRU)

    implicit none

    type(cru_type), intent(inout) :: CRU

    integer :: i

    write(*,*) 'Closing CRU files.'
    do i=1, CRU%nmet
       if (CRU%f_id(i) > -1) then
          errstatus = nf90_close(CRU%f_id(i))
          call handle_err(errstatus, "Closing CRU met file "//trim(CRU%MetFile(i)))
       end if
    end do

    if (CRU%Ndepf_id > -1) then
       errstatus = nf90_close(CRU%Ndepf_id)
       call handle_err(errstatus, "Closing CRU Ndep file.")
    end if

  end subroutine cru_close

END MODULE CABLE_CRU
