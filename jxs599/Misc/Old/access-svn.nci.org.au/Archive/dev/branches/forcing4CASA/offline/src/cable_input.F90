! cable_input.f90
!
! Input module for CABLE land surface scheme; 
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! The subroutines in this module, which are used both online and offline:
!   allocate_cable_vars - allocating all the cable variables with mp which
!                         is the total number of patches to be calculated
!
! The subroutines in this module, which are used online only:
!   get_default_lai - reads all monthly LAI from default gridded netcdf file
!
! The subroutines in this module, which are used offline only:
!   open_met_file - opens netcdf met forcing file and checks for variables
!   get_met_data - reads met (and LAI if present) from met file each time step
!   close_met_file - closes the netcdf met forcing file
!   load_parameters - decides where CABLE's parameters and init should load from
!   get_parameters_met - looks for and loads any of CABLE's parameters from
!                        met file
!
MODULE input_module   
! Note that any precision changes from r_1 to REAL(4) enable running with -r8
!
  USE abort_module, ONLY: abort, nc_abort
  USE define_dimensions, ONLY: mland,mp,r_1,i_d
  USE define_types
  USE casadimension, ONLY: icycle
  USE casavariable
  USE phenvariable
  USE physical_constants
  USE parameter_module
  USE checks_module, ONLY: ranges, rh_sh
  USE radiation_module, ONLY: sinbet
  USE io_variables
  USE read_module, ONLY: readpar
  USE initialisation_module
  USE netcdf ! link must be made in cd to netcdf-x.x.x/src/f90/netcdf.mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC get_default_lai, open_met_file, close_met_file,load_parameters, &
       allocate_cable_vars, get_met_data
  INTEGER(i_d) :: ncid_met ! met data netcdf file ID
  INTEGER(i_d) :: ncid_rain ! following are netcdf file IDs for gswp run
  INTEGER(i_d) :: ncid_snow
  INTEGER(i_d) :: ncid_lw
  INTEGER(i_d) :: ncid_sw
  INTEGER(i_d) :: ncid_ps
  INTEGER(i_d) :: ncid_qa
  INTEGER(i_d) :: ncid_ta
  INTEGER(i_d) :: ncid_wd
  ! - see ALMA compress by gathering
  INTEGER(i_d),POINTER,DIMENSION(:) :: landGrid ! for ALMA compressed variables
  REAL(r_1),POINTER,DIMENSION(:)  :: elevation ! site/grid cell elevation
  REAL(r_1),POINTER,DIMENSION(:)  :: avPrecip ! site/grid cell average precip
  TYPE met_varID_type 
     INTEGER(i_d) :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair,Qair,Rainf, &
          Snowf,CO2air,Elev,LAI,avPrecip,iveg,isoil,patchfrac
  END TYPE met_varID_type
  TYPE(met_varID_type) :: id ! netcdf variable IDs for input met variables
  TYPE met_units_type
     CHARACTER(LEN=20) :: SWdown,LWdown,Wind,Wind_E,PSurf,Tair, &
          Qair,Rainf,Snowf,CO2air,Elev,avPrecip
  END TYPE met_units_type
  TYPE(met_units_type) :: metunits ! units for meteorological variables
  TYPE convert_units_type
     REAL(r_1) :: PSurf,Tair,Qair,Rainf,CO2air,Elev
  END TYPE convert_units_type
  TYPE(convert_units_type) :: convert ! units change factors for met variables
  INTEGER(i_d) :: ok   ! netcdf error status

  !$OMP THREADPRIVATE(ok,exists)

CONTAINS

  ! =================================== LAI ====================================
  SUBROUTINE get_default_lai
  ! Fetches all default LAI data from global netcdf file
  !
  ! Input variables:
  !   filename%LAI      - via io_variables
  !   landpt            - via io_variables (%nap,cstart,cend,ilon,ilat)
  !   exists%laiPatch   - via io_variables
  ! Output variables:
  !   defaultLAI(mp,12) - via io_variables

  ! New input structure using netcdf with 1x1 resolution so that
  ! no more interpolation or translation from CCAM grid (BP apr2010)

    INTEGER(i_d) :: ncid
    INTEGER(i_d) :: xID, yID, pID, tID, laiID
    INTEGER(i_d) :: nlon, nlat, nLaiPatches, ntime
    INTEGER(i_d) :: e, tt ! do loop counter
    REAL(r_1), DIMENSION(:,:,:),  ALLOCATABLE :: inLai3D
    REAL(r_1), DIMENSION(:,:,:,:),ALLOCATABLE :: inLai4D

    ! Allocate default LAI variable: changed mland to mp (BP apr2010)
    ALLOCATE(defaultLAI(mp,12))   ! mp = mp_global

    WRITE(logn,*) ' Loading LAI from default file ', TRIM(filename%LAI)
    ! Open netcdf file
    ok = NF90_OPEN(filename%LAI,0,ncid) 
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error opening default LAI file.')

    ok = NF90_INQ_DIMID(ncid,'x',xID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
      ok = NF90_INQ_DIMID(ncid,'longitude',xID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring x dimension.')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid,xID,LEN=nlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid,'y',yID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
      ok = NF90_INQ_DIMID(ncid,'latitude',yID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring y dimension.')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid,yID,LEN=nlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting y dimension.')
    ok = NF90_INQ_DIMID(ncid,'patch',pID)
    IF(ok/=NF90_NOERR) THEN ! if failed
      exists%laiPatch = .FALSE.
      WRITE(logn,*) ' **ALL patches will be given the same LAI**'
    ELSE
      exists%laiPatch = .TRUE.
      ok = NF90_INQUIRE_DIMENSION(ncid,pID,len=nLaiPatches)
      IF(ANY(landpt(:)%nap>nLaiPatches)) THEN
        WRITE(logn,*) ' **Some patches will be given the same LAI**'
      END IF
      IF(nLaiPatches>max_vegpatches) THEN ! input file can have more info
        WRITE(*,*) ' Note that LAI input file has ', nLaiPatches, ' patches.'
        WRITE(*,*) ' while the model has max at ', max_vegpatches, ' patches.'
      END IF
    END IF
    ok = NF90_INQ_DIMID(ncid,'time',tID)
    IF (ok /= NF90_NOERR) THEN   ! added to read some input files (BP mar2011)
      ok = NF90_INQ_DIMID(ncid,'month',tID)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error inquiring t dimension.')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid,tID,LEN=ntime)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error getting time dimension.')
    IF (ntime /= 12) CALL abort('Time dimension not 12 months.')

    ok = NF90_INQ_VARID(ncid,'LAI',laiID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error finding LAI variable.')

    ! Read LAI values:
    IF (exists%laiPatch) THEN
      ALLOCATE(inLai4D(nlon,nlat,nLaiPatches,ntime))
      ok = NF90_GET_VAR(ncid,laiID,inLai4D)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading 4D LAI variable.')
      DO e = 1, mland  ! over all land grid points
        DO tt = 1, ntime
          defaultLAI(landpt(e)%cstart:landpt(e)%cend,tt) = &
                inLai4D(landpt(e)%ilon,landpt(e)%ilat,1:landpt(e)%nap,tt)
        END DO
      END DO
    ELSE
      ALLOCATE(inLai3D(nlon,nlat,ntime))
      ok = NF90_GET_VAR(ncid,laiID,inLai3D)
      IF (ok /= NF90_NOERR) CALL nc_abort(ok,'Error reading 3D LAI variable.')
      DO e = 1, mland
        DO tt = 1, ntime
          defaultLAI(landpt(e)%cstart:landpt(e)%cend,tt) = &
                inLai3D(landpt(e)%ilon,landpt(e)%ilat,tt)
        END DO
      END DO
    END IF

    ! Close netcdf file
    ok = NF90_CLOSE(ncid)

  END SUBROUTINE get_default_lai
  !=============================================================================
  SUBROUTINE open_met_file(dels,kend,spinup)
    ! Opens netcdf file containing meteorological (LSM input) data
    ! and determines:
    ! 1. Spatial details - number of sites/grid cells, latitudes, longitudes
    ! 2. Timing details - time step size, number of timesteps, starting date,
    !    and whether time coordinate is local or GMT
    ! 3. Checks availability, including units issues, of all required
    !    meteorological input variables. Also checks whether or not LAI is 
    !    present, and fetches prescribed veg ans soil type if present.
    REAL(r_1), INTENT(OUT) :: dels ! time step size
    INTEGER(i_d), INTENT(OUT) :: kend ! number of time steps in simulation
    LOGICAL, INTENT(IN) :: spinup ! will a model spinup be performed?
    INTEGER(i_d) :: timevarID ! time variable ID number
    INTEGER(i_d),DIMENSION(1) :: timedimID ! time dimension ID number
    INTEGER(i_d) :: xdimID,ydimID    ! x and y dimension ID numbers
    INTEGER(i_d) :: patchdimID ! patch dimension ID
    INTEGER(i_d) :: monthlydimID ! month dimension ID for LAI info
    INTEGER(i_d) :: maskID    ! mask variable ID
    INTEGER(i_d) :: landID    ! land variable ID
    INTEGER(i_d) :: landdimID ! land dimension ID
    INTEGER(i_d) :: latitudeID, longitudeID ! lat and lon variable IDs
    REAL(r_1),POINTER, DIMENSION(:) :: lat_temp, lon_temp ! lat and lon
    INTEGER,POINTER,DIMENSION(:) ::land_xtmp,land_ytmp ! temp indicies
    REAL(r_1)    :: tshod        ! temporary variable
    INTEGER(i_d) :: tsdoy,tsyear ! temporary variables
    REAL(r_1)    :: ehod ! end time hour-of-day
    INTEGER(i_d) :: edoy,eyear ! end time day-of-year and year
    INTEGER(i_d) :: jump_days ! days made by first "time" entry
    INTEGER(i_d) :: sdoytmp ! used to determine start time hour-of-day
    INTEGER(i_d) :: mland_ctr ! counter for number of land points read from file
    INTEGER(i_d) :: mland_fromfile ! number of land points in file 
    INTEGER(i_d) :: lai_dims ! number of dims of LAI var if in met file
    INTEGER(i_d) :: iveg_dims ! number of dims of iveg var if in met file
    INTEGER(i_d) :: isoil_dims ! number of dims of isoil var if in met file
    LOGICAL :: all_met ! ALL required met in met file (no synthesis)?
    CHARACTER(LEN=10) :: todaydate, nowtime ! used to timestamp log file
    REAL(4),DIMENSION(1) :: data1 ! temp variable for netcdf reading
    INTEGER(i_d),DIMENSION(1) :: data1i ! temp variable for netcdf reading
    REAL(4),DIMENSION(1,1) :: data2 ! temp variable for netcdf reading
    INTEGER(i_d),DIMENSION(1,1) :: data2i ! temp variable for netcdf reading
    REAL(4),POINTER,DIMENSION(:,:,:) :: tempPrecip3 ! used for spinup adj
    REAL(4),POINTER,DIMENSION(:,:) :: tempPrecip2 ! used for spinup adj
    REAL(4),POINTER,DIMENSION(:) :: temparray1 ! temp read in variable
    REAL(4),POINTER,DIMENSION(:,:) :: temparray2 ! temp read in variable
    INTEGER(i_d),DIMENSION(4) :: laidimids ! for checking lai variable
    REAL(r_1) :: precipTot ! used for spinup adj
    REAL(r_1) :: avPrecipInMet ! used for spinup adj
    INTEGER(i_d) :: x,y,i,j ! do loop counters
    INTEGER(i_d) :: tempmonth
    logical :: smod100, smod400a, smod400b, smod400

    ! Initialise parameter loading switch - will be set to TRUE when 
    ! parameters are loaded:
    exists%parameters = .FALSE. ! initialise
    ! Initialise initialisation loading switch - will be set to TRUE when 
    ! initialisation data are loaded:
    exists%initial = .FALSE. ! initialise

    ! Write filename to log file:
    WRITE(logn,*) '============================================================'
    WRITE(logn,*) 'Log file for offline CABLE run:'
    CALL DATE_AND_TIME(todaydate, nowtime)
    todaydate=todaydate(1:4)//'/'//todaydate(5:6)//'/'//todaydate(7:8)
    nowtime=nowtime(1:2)//':'//nowtime(3:4)//':'//nowtime(5:6)
    WRITE(logn,*) TRIM(nowtime),' ',TRIM(todaydate)
    WRITE(logn,*) '============================================================'

    ! Open netcdf file:
  IF (ncciy > 0) THEN
    WRITE(logn,*) 'Opening met data file: ', TRIM(gswpfile%rainf), ' and 7 more'
    ok = NF90_OPEN(gswpfile%rainf,0,ncid_rain)
    ok = NF90_OPEN(gswpfile%snowf,0,ncid_snow)
    ok = NF90_OPEN(gswpfile%LWdown,0,ncid_lw)
    ok = NF90_OPEN(gswpfile%SWdown,0,ncid_sw)
    ok = NF90_OPEN(gswpfile%PSurf,0,ncid_ps)
    ok = NF90_OPEN(gswpfile%Qair,0,ncid_qa)
    ok = NF90_OPEN(gswpfile%Tair,0,ncid_ta)
    ok = NF90_OPEN(gswpfile%wind,0,ncid_wd)
    ncid_met = ncid_rain
  ELSE
    WRITE(logn,*) 'Opening met data file: ', TRIM(filename%met)
    ok = NF90_OPEN(filename%met,0,ncid_met) ! open met data file
    IF (ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error opening netcdf met forcing file '//TRIM(filename%met)// &
         ' (SUBROUTINE open_met_file)') 
  ENDIF

    !!=====================VV Determine spatial details VV=================
    ! Determine number of sites/gridcells.
    ! Find size of 'x' or 'lat' dimension:
    ok = NF90_INQ_DIMID(ncid_met,'x', xdimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       ! Try 'lat' instead of x
       ok = NF90_INQ_DIMID(ncid_met,'lat', xdimID)
       IF(ok/=NF90_NOERR) CALL nc_abort &
            (ok,'Error finding x dimension in '&
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid_met,xdimID,len=xdimsize)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining size of x dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Find size of 'y' dimension:
    ok = NF90_INQ_DIMID(ncid_met,'y', ydimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       ! Try 'lon' instead of y
       ok = NF90_INQ_DIMID(ncid_met,'lon', ydimID)
       IF(ok/=NF90_NOERR) CALL nc_abort &
            (ok,'Error finding y dimension in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ok = NF90_INQUIRE_DIMENSION(ncid_met,ydimID,len=ydimsize)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining size of y dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Determine number of gridcells in netcdf file:
    ngridcells = xdimsize*ydimsize
    WRITE(logn,'(A28,I7)') 'Total number of gridcells: ', ngridcells

    ! Get all latitude and longitude values.
    ! Find latitude variable (try 'latitude' and 'nav_lat'(ALMA)):
    ok = NF90_INQ_VARID(ncid_met, 'latitude', latitudeID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_VARID(ncid_met, 'nav_lat', latitudeID)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding latitude variable in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ! Allocate space for lat_all variable and its temp counterpart:
    ALLOCATE(lat_all(xdimsize,ydimsize))
    ALLOCATE(temparray2(xdimsize,ydimsize))
    ! Get latitude values for entire region:
    ok= NF90_GET_VAR(ncid_met,latitudeID,temparray2)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading latitude variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Needed since r_1 will be double precision with -r8:
    lat_all = REAL(temparray2,r_1)
    ! Find longitude variable (try 'longitude' and 'nav_lon'(ALMA)):
    ok = NF90_INQ_VARID(ncid_met, 'longitude', longitudeID)
    IF(ok /= NF90_NOERR) THEN
       ok = NF90_INQ_VARID(ncid_met, 'nav_lon', longitudeID)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding longitude variable in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    END IF
    ! Allocate space for lon_all variable:
    ALLOCATE(lon_all(xdimsize,ydimsize))
    ! Get longitude values for entire region:
    ok= NF90_GET_VAR(ncid_met,longitudeID,temparray2)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading longitude variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Needed since r_1 will be double precision with -r8:
    lon_all = REAL(temparray2,r_1)

    ! Check for "mask" variable or "land" variable to tell grid type
    ! (and allow neither if only one gridpoint). "mask" is a 2D variable
    ! with dims x,y and "land" is a 1D variable.
    ok = NF90_INQ_VARID(ncid_met, 'mask', maskID) ! check for "mask"
    IF(ok /= NF90_NOERR) THEN ! if error, i.e. no "mask" variable:
       ! Check for "land" variable:
       ok = NF90_INQ_VARID(ncid_met, 'land', landID)
       IF(ok /= NF90_NOERR) THEN ! ie no "land" or "mask"
          IF(ngridcells==1) THEN 
             ! Allow no explicit grid system if only one gridpoint
             ALLOCATE(mask(xdimsize,ydimsize)) ! Allocate "mask" variable
             metGrid='mask' ! Use mask system, one gridpoint.
             mask = 1
             ALLOCATE(latitude(1),longitude(1))
             latitude = lat_all(1,1)
             longitude = lon_all(1,1)
             mland_fromfile=1
             ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
             land_x = 1
             land_y = 1
          ELSE
             ! Call abort if more than one gridcell and no
             ! recognised grid system:
             CALL nc_abort &
                  (ok,'Error finding grid system ("mask" or "land") variable in ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          END IF
       ELSE ! i.e. "land" variable exists
          metGrid='land'
          ! Check size of "land" dimension:
          ok = NF90_INQ_DIMID(ncid_met,'land', landdimID)
          IF(ok/=NF90_NOERR) CALL nc_abort &
               (ok,'Error finding land dimension in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ok = NF90_INQUIRE_DIMENSION(ncid_met,landdimID,len=mland_fromfile)
          IF(ok/=NF90_NOERR) CALL nc_abort &
               (ok,'Error determining size of land dimension in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Allocate landGrid variable and its temporary counterpart:
          ALLOCATE(landGrid(mland_fromfile))
          ALLOCATE(temparray1(mland_fromfile))
          ! Get values of "land" variable from file:
          ok= NF90_GET_VAR(ncid_met,landID,temparray1)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading "land" variable in ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Needed since r_1 will be double precision with -r8:
          landGrid = REAL(temparray1,r_1)
          DEALLOCATE(temparray1)
          ! Allocate latitude and longitude variables:
          ALLOCATE(latitude(mland_fromfile),longitude(mland_fromfile))
          ! Write to indicies of points in all-grid which are land
          ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
          ! Allocate "mask" variable:
          ALLOCATE(mask(xdimsize,ydimsize))
          ! Initialise all gridpoints as sea:
          mask = 0
          DO j=1, mland_fromfile ! over all land points
             ! Find x and y coords of current land point
             y = INT((landGrid(j)-1)/xdimsize)
             x = landGrid(j) - y * xdimsize
             y=y+1
             ! Write lat and lon to land-only lat/lon vars:
             latitude(j) = lat_all(x,y)
             longitude(j) = lon_all(x,y)
             ! Write to mask variable:
             mask(x,y)=1
             ! Save indicies:
             land_x(j) = x
             land_y(j) = y
          END DO
       END IF ! does "land" variable exist 
    ELSE ! i.e. "mask" variable exists
       ! Allocate "mask" variable:
       ALLOCATE(mask(xdimsize,ydimsize))
       metGrid='mask' ! Use mask system
       ! Get mask values from file:
       ok= NF90_GET_VAR(ncid_met,maskID,mask)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error reading "mask" variable in ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       ! Allocate space for extracting land lat/lon values:
       ALLOCATE(lat_temp(ngridcells),lon_temp(ngridcells))
       ! Allocate space for extracting index of mask which is land
       ALLOCATE(land_xtmp(ngridcells),land_ytmp(ngridcells))
       ! Cycle through all gridsquares:
       mland_ctr = 0 ! initialise
       DO y=1,ydimsize
          DO x=1,xdimsize
             IF(mask(x,y)==1) THEN ! If land
                mland_ctr = mland_ctr + 1
                ! Store lat and lon for land points
                lat_temp(mland_ctr) = lat_all(x,y)
                lon_temp(mland_ctr) = lon_all(x,y)
                ! Store indicies of points in mask which are land
                land_xtmp(mland_ctr) = x
                land_ytmp(mland_ctr) = y
             END IF
          END DO
       END DO
       ! Record number of land points
       mland_fromfile = mland_ctr
       ! Allocate latitude and longitude variables:
       ALLOCATE(latitude(mland_fromfile),longitude(mland_fromfile))
       ! Write to latitude and longitude variables:
       latitude = lat_temp(1:mland_fromfile)
       longitude = lon_temp(1:mland_fromfile)
       ! Write to indicies of points in mask which are land
       ALLOCATE(land_x(mland_fromfile),land_y(mland_fromfile))
       land_x = land_xtmp(1:mland_fromfile)
       land_y = land_ytmp(1:mland_fromfile)
       ! Clear lon_temp, lat_temp,land_xtmp,land_ytmp
       DEALLOCATE(lat_temp,lon_temp,land_xtmp,land_ytmp)
    END IF ! "mask" variable or no "mask" variable

    ! Set global mland value (number of land points), used to allocate
    ! all of CABLE's arrays:
    mland = mland_fromfile

    ! Write number of land points to log file:
    WRITE(logn,'(24X,I7,A29)') mland_fromfile, ' of which are land grid cells'

    ! Check if veg/soil patch dimension exists (could have
    ! parameters with patch dimension)
    ok = NF90_INQ_DIMID(ncid_met,'patch', patchdimID)
    IF(ok/=NF90_NOERR) THEN ! if failed
       exists%patch = .FALSE.
       nmetpatches = 1       ! initialised so that old met files without patch
       ! data can still be run correctly (BP apr08)
    ELSE ! met file does have patch dimension
       exists%patch = .TRUE.
       ok = NF90_INQUIRE_DIMENSION(ncid_met,patchdimID,len=nmetpatches)
    END IF

    ! Check if monthly dimension exists for LAI info
    ok = NF90_INQ_DIMID(ncid_met,'monthly', monthlydimID)
    IF(ok==NF90_NOERR) THEN ! if found
       ok = NF90_INQUIRE_DIMENSION(ncid_met,monthlydimID,len=tempmonth)
       IF(tempmonth/=12) CALL abort ('Number of months in met file /= 12.')
    END IF

    ! Set longitudes to be [-180,180]:
    WHERE(longitude>180.0) 
       longitude = longitude - 360.0
    END WHERE
    ! Check ranges for latitude and longitude:
    IF(ANY(longitude>180.0).OR.ANY(longitude<-180.0)) &
         CALL abort('Longitudes read from '//TRIM(filename%met)// &
         ' are not [-180,180] or [0,360]! Please set.')
    IF(ANY(latitude>90.0).OR.ANY(latitude<-90.0)) &
         CALL abort('Latitudes read from '//TRIM(filename%met)// &
         ' are not [-90,90]! Please set.')

    !!=================^^ End spatial details ^^========================

    !!=========VV Determine simulation timing details VV================
    ! Inquire 'time' variable's ID:
    ok = NF90_INQ_VARID(ncid_met, 'time', timevarID)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding time variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get ID for dimension upon which time depends:
    ok = NF90_INQUIRE_VARIABLE(ncid_met,timevarID,dimids=timedimID)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining "time" dimension dimension in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Determine number of time steps:
    ok = NF90_INQUIRE_DIMENSION(ncid_met,timedimID(1),len=kend)
    IF(ok/=NF90_NOERR) CALL nc_abort &
         (ok,'Error determining number of timesteps in ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Allocate time variable:
    ALLOCATE(timevar(kend))
    ! Fetch 'time' variable:
    ok= NF90_GET_VAR(ncid_met,timevarID,timevar)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error reading time variable in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Set time step size:
    dels = REAL(timevar(2) - timevar(1),r_1)
    WRITE(logn,'(1X,A29,I8,A3,F10.3,A5)') 'Number of time steps in run: ',&
         kend,' = ', REAL(kend)/(3600/dels*24),' days'
    !********* gswp input file has bug in timevar **************
    IF (ncciy > 0) THEN
      PRINT *, 'original timevar(kend) = ', timevar(kend)
      DO i = 1, kend - 1
        timevar(i+1) = timevar(i) + dels
      ENDDO
      PRINT *, 'New      timevar(kend) = ', timevar(kend)
!! hacking (BP feb2011)
!      kend = 16   ! 2 days for 1986
!!      kend = 480  ! 2 months for 1986
!      PRINT *, 'Hacked   timevar(kend) = ', timevar(kend)
!      PRINT *, 'Hacked kend = ', kend
!! end hacking
    END IF
    !********* done bug fixing for timevar in gswp input file **
    ! Write time step size to log file:
    WRITE(logn,'(1X,A17,F8.1,1X,A7)') 'Time step size:  ', dels, 'seconds'
    ! Get units for 'time' variable:
    ok = NF90_GET_ATT(ncid_met,timevarID,'units',timeunits)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding time variable units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    !********* gswp input file has bug in timeunits ************
    IF (ncciy > 0) WRITE(timeunits(26:27),'(i2.2)') 0
    !********* done bug fixing for timeunits in gwsp file ******
    WRITE(logn,*) 'Time variable units: ', timeunits
    ! Get coordinate field:
    ok = NF90_GET_ATT(ncid_met,timevarID,'coordinate',time_coord)
    ! If error getting coordinate field (i.e. it doesn't exist):
    IF(ok /= NF90_NOERR) THEN
       ! Assume default time coordinate:
       IF(mland_fromfile==1) THEN ! If single site, this is local time
          time_coord = 'LOC' ! 12am is 12am local time, at site/gridcell
       ELSE ! If multiple/global/regional, use GMT
          time_coord = 'GMT' ! 12am is GMT time, local time set by longitude
       END IF
    ELSE IF((ok==NF90_NOERR.AND.time_coord=='LOC'.AND.mland_fromfile>1)) THEN
       ! Else if local time is selected for regional simulation, abort:
       CALL abort('"time" variable must be GMT for multiple site simulation!' &
            //' Check "coordinate" field in time variable.' &
            //' (SUBROUTINE open_met_file)')
    ELSE IF(time_coord/='LOC'.AND.time_coord/='GMT') THEN
       CALL abort('Meaningless time coordinate in met data file!' &
            // ' (SUBROUTINE open_met_file)')
    END IF

    ! Use internal files to convert "time" variable units (giving the run's 
    ! start time) from character to integer; calculate starting hour-of-day,
    ! day-of-year, year:
    READ(timeunits(15:18),*) syear
    READ(timeunits(20:21),*) smoy ! integer month
    READ(timeunits(23:24),*) sdoytmp ! integer day of that month
    READ(timeunits(26:27),*) shod  ! starting hour of day 
    ! Decide day-of-year for non-leap year:
    SELECT CASE(smoy)
    CASE(1) ! Jan
       sdoy=sdoytmp
    CASE(2) ! Feb
       sdoy=sdoytmp+lastday(1)
    CASE(3) ! Mar
       sdoy=sdoytmp+lastday(2)
    CASE(4)
       sdoy=sdoytmp+lastday(3)
    CASE(5)
       sdoy=sdoytmp+lastday(4)
    CASE(6)
       sdoy=sdoytmp+lastday(5)
    CASE(7)
       sdoy=sdoytmp+lastday(6)
    CASE(8)
       sdoy=sdoytmp+lastday(7)
    CASE(9)
       sdoy=sdoytmp+lastday(8)
    CASE(10)
       sdoy=sdoytmp+lastday(9)
    CASE(11)
       sdoy=sdoytmp+lastday(10)
    CASE(12) 
       sdoy=sdoytmp+lastday(11)
    CASE DEFAULT
       CALL abort('Could not interpret month in "time" units from ' &
            //TRIM(filename%met)// '(SUBROUTINE open_met_file)')
    END SELECT
    IF(leaps) THEN ! If we're using leap year timing:
       ! If start year is a leap year and start month > Feb, add a day:
       IF( MOD(syear,4)==0.AND.MOD(syear,100)/=0 ) smod100 = .true.
       IF( MOD(syear,4)==0.AND.MOD(syear,400)==0 ) smod400a = .true.
       if (smoy>2) smod400b = .true.
       if (smod400a .and. smod400b) smod400 = .true.       
       if(smod100 .or. smod400) sdoy = sdoy + 1
                 
       ! Number of days between start position and 1st timestep:
       jump_days = INT((timevar(1)/3600.0 + shod)/24.0)
       ! Cycle through days to find leap year inclusive starting date:
       DO i=1,jump_days
          sdoy = sdoy + 1
          IF((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)) THEN
             ! Set month of year for leap year:
             SELECT CASE(sdoy)
             CASE(1) ! Jan
                smoy = 1
             CASE(32) ! Feb
                smoy = 2
             CASE(61) ! Mar
                smoy = 3
             CASE(92)
                smoy = 4
             CASE(122)
                smoy = 5
             CASE(153)
                smoy = 6
             CASE(183)
                smoy = 7
             CASE(214)
                smoy = 8
             CASE(245)
                smoy = 9
             CASE(275)
                smoy = 10
             CASE(306)
                smoy = 11
             CASE(336) 
                smoy = 12
             CASE(367)! end of year; increment
                syear = syear + 1 
                smoy = 1
                sdoy = 1
             END SELECT
          ELSE 
             ! Set month of year for non-leap year:
             SELECT CASE(sdoy)
             CASE(1) ! Jan
                smoy = 1
             CASE(32) ! Feb
                smoy = 2
             CASE(60) ! Mar
                smoy = 3
             CASE(91)
                smoy = 4
             CASE(121)
                smoy = 5
             CASE(152)
                smoy = 6
             CASE(182)
                smoy = 7
             CASE(213)
                smoy = 8
             CASE(244)
                smoy = 9
             CASE(274)
                smoy = 10
             CASE(305)
                smoy = 11
             CASE(335) 
                smoy = 12
             CASE(366) ! end of year; increment
                syear = syear + 1 
                smoy = 1
                sdoy = 1
             END SELECT
          END IF
       END DO
       ! Update starting hour-of-day for first time step's value
       shod = MOD(REAL(timevar(1)/3600.0 + shod),24.0)
    ELSE ! If not using leap year timing,
       ! simply update starting times for first value of "time":
       tshod = MOD(REAL(timevar(1)/3600.0 + shod),24.0)
       tsdoy = MOD(INT((timevar(1)/3600.0 + shod)/24.0) + sdoy, 365)
       tsyear = INT(REAL(INT((timevar(1)/3600.0+shod)/24.0)+sdoy)/365.0)+syear
       shod=tshod  ! real valued
       sdoy=tsdoy  ! integer valued
       syear=tsyear ! integer valued
       ! Set moy:
       SELECT CASE(sdoy)
       CASE(1:31) ! Jan
          smoy = 1
       CASE(32:59) ! Feb
          smoy = 2
       CASE(60:90) ! Mar
          smoy = 3
       CASE(91:120)
          smoy = 4
       CASE(121:151)
          smoy = 5
       CASE(152:181)
          smoy = 6
       CASE(182:212)
          smoy = 7
       CASE(213:243)
          smoy = 8
       CASE(244:273)
          smoy = 9
       CASE(274:304)
          smoy = 10
       CASE(305:334)
          smoy = 11
       CASE(335:365) 
          smoy = 12
       END SELECT
    END IF
    ! Now all start time variables established, report to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') &
         'Run begins: ',shod,' hour-of-day, ',sdoy, ' day-of-year, ',&
         syear, time_coord, 'time'
    ! Determine ending time of run...
    IF(leaps) THEN ! If we're using leap year timing...
       ! Number of days between beginning and end of run:
       jump_days = INT(((timevar(kend)-timevar(1)+dels)/3600.0 + shod)/24.0)
       ! initialise:
       ehod = shod
       edoy = sdoy
       eyear = syear
       ! Cycle through days to find leap year inclusive ending date:
       DO i=1,jump_days
          edoy = edoy + 1
          IF((MOD(eyear,4)==0.AND.MOD(eyear,100)/=0).OR. & 
               (MOD(eyear,4)==0.AND.MOD(eyear,400)==0)) THEN
             ! Set moy for leap year:
             SELECT CASE(edoy)
             CASE(367)! end of year; increment
                eyear = eyear + 1 
                edoy = 1
             END SELECT
          ELSE 
             ! Set moy for non-leap year:
             SELECT CASE(edoy)
             CASE(366) ! end of year; increment
                eyear = eyear + 1 
                edoy = 1
             END SELECT
          END IF
       END DO
       ! Update starting hour-of-day fot first time step's value
       ehod = MOD(REAL((timevar(kend)-timevar(1)+dels)/3600.0 + shod),24.0)
    ELSE ! if not using leap year timing
       ! Update shod, sdoy, syear for first "time" value:
       ehod = MOD(REAL((timevar(kend)-timevar(1)+dels)/3600.0 + shod),24.0)
       edoy = MOD(INT(((timevar(kend)-timevar(1)+dels)/3600.0 + shod)/24.0) &
            + sdoy, 365)
       eyear = INT(REAL(INT(((timevar(kend)-timevar(1)+dels) &
            /3600.0+shod)/24.0)+sdoy)/365.0)+syear
    END IF
    ! Report finishing time to log file:
    WRITE(logn,'(1X,A12,F5.2,A14,I3,A14,I4,2X,A3,1X,A4)') 'Run ends:   ',&
         ehod,' hour-of-day, ',edoy, &
         ' day-of-year, ', eyear, time_coord, 'time'
    !!===================^^ End timing details ^^==========================

    !!===================VV Look for met variables VV======================
    all_met = .TRUE. ! initialise
    ! Look for SWdown (essential):- - - - - - - - - - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_sw
    ok = NF90_INQ_VARID(ncid_met,'SWdown',id%SWdown)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding SWdown in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get SWdown units and check okay:
    ok = NF90_GET_ATT(ncid_met,id%SWdown,'units',metunits%SWdown)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding SWdown units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%SWdown(1:4)/='W/m2'.AND.metunits%SWdown(1:5) &
         /='W/m^2'.AND.metunits%SWdown(1:5)/='Wm^-2' &
         .AND.metunits%SWdown(1:4)/='Wm-2') THEN
       WRITE(*,*) metunits%SWdown
       CALL abort('Unknown units for SWdown'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Tair (essential):- - - - - - - - - - - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_ta
    ok = NF90_INQ_VARID(ncid_met,'Tair',id%Tair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Tair in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Tair units and check okay:
    ok = NF90_GET_ATT(ncid_met,id%Tair,'units',metunits%Tair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Tair units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Tair(1:1)=='C'.OR.metunits%Tair(1:1)=='c') THEN
       ! Change from celsius to kelvin:
       convert%Tair = tfrz
       WRITE(logn,*) 'Temperature will be converted from C to K'
    ELSE IF(metunits%Tair(1:1)=='K'.OR.metunits%Tair(1:1)=='k') THEN
       ! Units are correct
       convert%Tair = 0.0
    ELSE
       WRITE(*,*) metunits%Tair
       CALL abort('Unknown units for Tair'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Qair (essential):- - - - - - - - - - - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_qa
    ok = NF90_INQ_VARID(ncid_met,'Qair',id%Qair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Qair in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Qair units:
    ok = NF90_GET_ATT(ncid_met,id%Qair,'units',metunits%Qair)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Qair units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Qair(1:1)=='%'.OR.metunits%Qair(1:1)=='-') THEN
       ! Change from relative humidity to specific humidity:
       convert%Qair = -999.0
       WRITE(logn,*) 'Humidity will be converted from relative to specific'
    ELSE IF(metunits%Qair(1:3)=='g/g'.OR.metunits%Qair(1:5)=='kg/kg' &
         .OR.metunits%Qair(1:3)=='G/G'.OR.metunits%Qair(1:5)=='KG/KG') THEN
       ! Units are correct
       convert%Qair=1.0
    ELSE
       WRITE(*,*) metunits%Qair
       CALL abort('Unknown units for Qair'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Look for Rainf (essential):- - - - - - - - - - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_rain
    ok = NF90_INQ_VARID(ncid_met,'Rainf',id%Rainf)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Rainf in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    ! Get Rainf units:
    ok = NF90_GET_ATT(ncid_met,id%Rainf,'units',metunits%Rainf)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Rainf units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Rainf(1:8)=='kg/m^2/s'.OR.metunits%Rainf(1:10)== &
         'kgm^-2s^-1'.OR.metunits%Rainf(1:4)=='mm/s'.OR. &
         metunits%Rainf(1:6)=='mms^-1'.OR. &
         metunits%Rainf(1:7)=='kg/m^2s') THEN
       ! Change from mm/s to mm/time step:
       convert%Rainf = dels
    ELSE IF(metunits%Rainf(1:4)=='mm/h'.OR.metunits%Rainf(1:6)== &
         'mmh^-1') THEN
       ! Change from mm/h to mm/time step:
       convert%Rainf = dels/3600.0
    ELSE
       WRITE(*,*) metunits%Rainf
       CALL abort('Unknown units for Rainf'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Multiply acceptable Rainf ranges by time step size:
    ranges%Rainf = ranges%Rainf*dels ! range therefore depends on dels
    ! Look for Wind (essential):- - - - - - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_wd
    ok = NF90_INQ_VARID(ncid_met,'Wind',id%Wind)
    IF(ok /= NF90_NOERR) THEN
       ! Look for vector wind:
       ok = NF90_INQ_VARID(ncid_met,'Wind_N',id%Wind)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Wind in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       ok = NF90_INQ_VARID(ncid_met,'Wind_E',id%Wind_E)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Wind_E in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       exists%Wind = .FALSE. ! Use vector wind when reading met
    ELSE
       exists%Wind = .TRUE. ! 'Wind' variable exists
    END IF
    ! Get Wind units:
    ok = NF90_GET_ATT(ncid_met,id%Wind,'units',metunits%Wind)
    IF(ok /= NF90_NOERR) CALL nc_abort &
         (ok,'Error finding Wind units in met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
    IF(metunits%Wind(1:3)/='m/s'.AND.metunits%Wind(1:2)/='ms') THEN
       WRITE(*,*) metunits%Wind
       CALL abort('Unknown units for Wind'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    END IF
    ! Now "optional" variables:
    ! Look for LWdown (can be synthesised):- - - - - - - - - - - - - - -
    IF (ncciy > 0) ncid_met = ncid_lw
    ok = NF90_INQ_VARID(ncid_met,'LWdown',id%LWdown)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%LWdown = .TRUE. ! LWdown is present in met file
       ! Get LWdown units and check okay:
       ok = NF90_GET_ATT(ncid_met,id%LWdown,'units',metunits%LWdown)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding LWdown units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       IF(metunits%LWdown(1:4)/='W/m2'.AND.metunits%LWdown(1:5) &
            /='W/m^2'.AND.metunits%LWdown(1:5)/='Wm^-2' &
            .AND.metunits%LWdown(1:4)/='Wm-2') THEN
          WRITE(*,*) metunits%LWdown
          CALL abort('Unknown units for LWdown'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE
       exists%LWdown = .FALSE. ! LWdown is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,*) 'LWdown not present in met file; ', &
            'values will be synthesised based on air temperature.'
    END IF
    ! Look for PSurf (can be synthesised):- - - - - - - - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_ps
    ok = NF90_INQ_VARID(ncid_met,'PSurf',id%PSurf)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%PSurf = .TRUE. ! PSurf is present in met file
       ! Get PSurf units and check:
       ok = NF90_GET_ATT(ncid_met,id%PSurf,'units',metunits%PSurf)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding PSurf units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       IF(metunits%PSurf(1:2)=='Pa'.OR.metunits%PSurf(1:2)=='pa'.OR. &
            metunits%PSurf(1:2)=='PA' ) THEN
          ! Change from pa to mbar (cable uses mbar):
          convert%PSurf = 0.01
          WRITE(logn,*) 'Pressure will be converted from Pa to mb'
       ELSE IF(metunits%PSurf(1:2)=='KP'.OR.metunits%PSurf(1:2)=='kP' &
            .OR.metunits%PSurf(1:2)=='Kp'.OR.metunits%PSurf(1:2)=='kp') THEN
          ! convert from kPa to mb
          convert%PSurf = 10.0
          WRITE(logn,*) 'Pressure will be converted from kPa to mb'
       ELSE IF(metunits%PSurf(1:2)=='MB'.OR.metunits%PSurf(1:2)=='mB' &
            .OR.metunits%PSurf(1:2)=='Mb'.OR.metunits%PSurf(1:2)=='mb') THEN
          ! Units are correct
          convert%PSurf = 1.0
       ELSE
          WRITE(*,*) metunits%PSurf
          CALL abort('Unknown units for PSurf'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! If PSurf not present
       exists%PSurf = .FALSE. ! PSurf is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Look for "elevation" variable to approximate pressure based
       ! on elevation and temperature:
       ok = NF90_INQ_VARID(ncid_met,'Elevation',id%Elev)
       IF(ok == NF90_NOERR) THEN ! elevation present
          ! Get elevation units:
          ok = NF90_GET_ATT(ncid_met,id%Elev,'units',metunits%Elev)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error finding elevation units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Units should be metres or feet:
          IF(metunits%Elev(1:1)=='m'.OR.metunits%Elev(1:1)=='M') THEN
             ! This is the expected unit - metres
             convert%Elev = 1.0
          ELSE IF(metunits%Elev(1:1)=='f'.OR.metunits%Elev(1:1)=='F') THEN
             ! Convert from feet to metres:
             convert%Elev = 0.3048
          ELSE
             CALL abort('Unknown units for Elevation'// &
                  ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
          END IF
          ! Allocate space for elevation variable:
          ALLOCATE(elevation(mland))
          ! Get site elevations:
          IF(metGrid=='mask') THEN
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%Elev,data2, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading elevation in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                elevation(i)=REAL(data2(1,1),r_1)*convert%Elev
             END DO
          ELSE IF(metGrid=='land') THEN
             ! Collect data from land only grid in netcdf file:
             ok= NF90_GET_VAR(ncid_met,id%Elev,data1)
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading elevation in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             elevation = REAL(data1,r_1) * convert%Elev
          END IF
       ELSE ! If both PSurf and elevation aren't present, abort:
          CALL abort &
               ('Error finding PSurf or Elevation in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       END IF
       ! Note static pressure based on elevation in log file:
       WRITE(logn,*) 'PSurf not present in met file; values will be ', &
            'synthesised based on elevation and temperature.'
    END IF
    ! Look for CO2air (can be assumed to be static):- - - - - - - - - - -
    ok = NF90_INQ_VARID(ncid_met,'CO2air',id%CO2air)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%CO2air = .TRUE. ! CO2air is present in met file
       ! Get CO2air units:
       ok = NF90_GET_ATT(ncid_met,id%CO2air,'units',metunits%CO2air)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding CO2air units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       IF(metunits%CO2air(1:3)/='ppm') THEN
          WRITE(*,*) metunits%CO2air
          CALL abort('Unknown units for CO2air'// &
               ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
       END IF
    ELSE ! CO2 not present
       exists%CO2air = .FALSE. ! CO2air is not present in met file
       all_met=.FALSE. ! not all met variables are present in file
       ! Note this in log file:
       WRITE(logn,'(A33,A24,I4,A5)') ' CO2air not present in met file; ', &
            'values will be fixed at ',INT(fixedCO2),' ppmv'
    END IF
    ! Look for Snowf (could be part of Rainf variable):- - - - - - - - - - 
    IF (ncciy > 0) ncid_met = ncid_snow
    ok = NF90_INQ_VARID(ncid_met,'Snowf',id%Snowf)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%Snowf = .TRUE. ! Snowf is present in met file
       ! Get Snowf units:
       ok = NF90_GET_ATT(ncid_met,id%Snowf,'units',metunits%Snowf)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding Snowf units in met data file ' &
            //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
       ! Make sure Snowf units are the same as Rainf units:
       IF(metunits%Rainf/=metunits%Snowf) CALL abort &
            ('Please ensure Rainf and Snowf units are the same'// &
            ' in '//TRIM(filename%met)//' (SUBROUTINE open_met_data)')
    ELSE
       exists%Snowf = .FALSE. ! Snowf is not present in met file
       !  all_met=.FALSE. not required; Snowf assumed to be in Rainf
       ! Note this in log file:
       WRITE(logn,*) 'Snowf not present in met file; ', &
            'Assumed to be contained in Rainf variable'
    END IF
    ! Look for LAI - - - - - - - - - - - - - - - - - - - - - - - - -
    ok = NF90_INQ_VARID(ncid_met,'LAI',id%LAI)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       exists%LAI = .TRUE. ! LAI is present in met file
       ! LAI will be read in which ever land grid is used
       ! Check dimension of LAI variable:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%LAI, &
            ndims=lai_dims,dimids=laidimids)
       ! If any of LAI's dimensions are the time dimension
       IF(ANY(laidimids==timedimID(1))) THEN
          exists%LAI_T = .TRUE. ! i.e. time varying LAI
          WRITE(logn,*) 'LAI found in met file - time dependent;'
       ELSE
          exists%LAI_T = .FALSE. ! i.e. not time varying LAI
       END IF
       IF(ANY(laidimids==monthlydimID)) THEN
          exists%LAI_M = .TRUE. ! i.e. time varying LAI, but monthly only
          WRITE(logn,*) 'LAI found in met file - monthly values;'
       ELSE
          exists%LAI_M = .FALSE.
       END IF
       IF(ANY(laidimids==patchdimID)) THEN
          exists%LAI_P = .TRUE. ! i.e. patch varying LAI
          WRITE(logn,*) 'LAI found in met file - patch-specific values'
       ELSE
          exists%LAI_P = .FALSE. ! i.e. not patch varying LAI
       END IF
    ELSE
       exists%LAI = .FALSE. ! LAI is not present in met file
       ! Report to log file
       WRITE(logn,*) 'LAI not present in met file; ', &
            'Will use MODIS coarse grid monthly LAI'
    END IF
    ! If a spinup is to be performed:
    IF(spinup) THEN
       ! Look for avPrecip variable (time invariant - used for spinup):
       ok = NF90_INQ_VARID(ncid_met,'avPrecip',id%avPrecip)
       IF(ok == NF90_NOERR) THEN ! If inquiry is okay and avPrecip exists
          ! Report to log file than modified spinup will be used:
          WRITE(logn,*) 'Spinup will use modified precip - avPrecip variable found'
          WRITE(logn,*) '  precip will be rescaled to match these values during spinup:'
          WRITE(*,*) 'Spinup will use modified precip - avPrecip variable found'
          WRITE(*,*) '  precip will be rescaled to match these values during spinup'
          ! Spinup will modify precip values:
          exists%avPrecip = .TRUE.
          ! Get avPrecip units:
          ok = NF90_GET_ATT(ncid_met,id%avPrecip,'units',metunits%avPrecip)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error finding avPrecip units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          IF(metunits%avPrecip(1:2)/='mm') CALL abort( &
               'Unknown avPrecip units in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
          ! Allocate space for avPrecip variable:
          ALLOCATE(avPrecip(mland))
          ! Get avPrecip from met file:
          IF(metGrid=='mask') THEN
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%avPrecip,data2, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading avPrecip in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                avPrecip(i)=REAL(data2(1,1),r_1)
             END DO
          ELSE IF(metGrid=='land') THEN
             ! Allocate single preciaion temporary variable:
             ALLOCATE(temparray1(mland))
             ! Collect data from land only grid in netcdf file:
             ok= NF90_GET_VAR(ncid_met,id%avPrecip,temparray1)
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading avPrecip in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             ! Needed since r_1 will be double precision with -r8:
             avPrecip = REAL(temparray1,r_1)
          END IF
          ! Now find average precip from met data, and create rescaling
          ! factor for spinup:
          ALLOCATE(PrecipScale(mland))
          DO i = 1, mland
             IF(metGrid=='mask') THEN
                ! Allocate space for temporary precip variable:
                ALLOCATE(tempPrecip3(1,1,kend))
                ! Get rainfall data for this grid cell:
                ok= NF90_GET_VAR(ncid_met,id%Rainf,tempPrecip3, &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,kend/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading Rainf in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Store total Rainf for this grid cell:
                PrecipTot = REAL(SUM(SUM(SUM(tempPrecip3,3),2)),r_1) &
                     * convert%Rainf 
                ! Get snowfall data for this grid cell:
                IF(exists%Snowf) THEN
                   ok= NF90_GET_VAR(ncid_met,id%Snowf,tempPrecip3, &
                        start=(/land_x(i),land_y(i),1/),count=(/1,1,kend/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading Snowf in met data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                   ! Add total Snowf to this grid cell total:
                   PrecipTot = PrecipTot + &
                        (REAL(SUM(SUM(SUM(tempPrecip3,3),2)),r_1) &
                        * convert%Rainf)
                END IF
                DEALLOCATE(tempPrecip3)
             ELSE IF(metGrid=='land') THEN
                ! Allocate space for temporary precip variable:
                ALLOCATE(tempPrecip2(1,kend))
                ! Get rainfall data for this land grid cell:
                ok= NF90_GET_VAR(ncid_met,id%Rainf,tempPrecip2, &
                     start=(/i,1/),count=(/1,kend/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading Rainf in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Store total Rainf for this land grid cell:
                PrecipTot = REAL(SUM(SUM(tempPrecip2,2)),r_1)*convert%Rainf 
                IF(exists%Snowf) THEN
                   ok= NF90_GET_VAR(ncid_met,id%Snowf,tempPrecip2, &
                        start=(/i,1/),count=(/1,kend/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading Snowf in met data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                   ! Add total Snowf to this land grid cell total:
                   PrecipTot = PrecipTot + (REAL(SUM(SUM(tempPrecip2,2)),r_1) &
                        * convert%Rainf)
                END IF
                DEALLOCATE(tempPrecip2)
             END IF
             ! Create rescaling factor for this grid cell to ensure spinup
             ! rainfall/snowfall is closer to average rainfall:
             ! First calculate annual average precip in met data:
             avPrecipInMet = PrecipTot/REAL(kend) * 3600.0/dels * 365 * 24
             PrecipScale(i) = avPrecipInMet/avPrecip(i)
             WRITE(logn,*) '  Site number:',i
             WRITE(logn,*) '  average precip quoted in avPrecip variable:', &
                  avPrecip(i)
             WRITE(logn,*) '  average precip in met data:',avPrecipInMet
          END DO ! over each land grid cell 
          DEALLOCATE(avPrecip)
       ELSE ! avPrecip doesn't exist in met file
          ! Spinup will not modify precip values:
          exists%avPrecip = .FALSE.
          WRITE(logn,*) 'Spinup will repeat entire data set until states converge'
          WRITE(logn,*) '  (see below for convergence criteria);'
          WRITE(*,*) 'Spinup will repeat entire data set until states converge:'
       END IF
    END IF  ! if a spinup is to be performed

    ! Look for veg type - - - - - - - - - - - - - - - - -:
    ok = NF90_INQ_VARID(ncid_met,'iveg',id%iveg)
    IF(ok == NF90_NOERR) THEN ! If 'iveg' exists in the met file
       ! Note existence of at least one model parameter in the met file:
       exists%parameters = .TRUE.
       ! Allocate space for user-defined veg type variable:
       ALLOCATE(vegtype_metfile(mland,nmetpatches))
       ! Check dimension of veg type:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%iveg,ndims=iveg_dims)
       IF(metGrid=='mask') THEN ! i.e. at least two spatial dimensions
          IF(iveg_dims==2) THEN ! no patch specific iveg information, just x,y
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%iveg,data2i, & ! get iveg data
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                vegtype_metfile(i,:)=data2i(1,1)
             END DO
          ELSE IF(iveg_dims==3) THEN ! i.e. patch specific iveg information
             ! Patch-specific iveg variable MUST be accompanied by 
             ! patchfrac variable with the same dimensions. So,
             ! Make sure that the patchfrac variable exists:
             ok = NF90_INQ_VARID(ncid_met,'patchfrac',id%patchfrac)
             IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                  (ok,'Patch-specific vegetation type (iveg) must be accompanied '// &
                  'by a patchfrac variable - this was not found in met data file '&
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             DO i = 1, mland
                ! Then, get the patch specific iveg data:
                ok= NF90_GET_VAR(ncid_met,id%iveg,vegtype_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       ELSE IF(metGrid=='land') THEN
          ! Collect data from land only grid in netcdf file:
          IF(iveg_dims==1) THEN ! i.e. no patch specific iveg information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%iveg,data1i, &
                     start=(/i/),count=(/1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                vegtype_metfile(i,:) = data1i(1)
             END DO
          ELSE IF(iveg_dims==2) THEN ! i.e. patch specific iveg information
             ! Patch-specific iveg variable MUST be accompanied by 
             ! patchfrac variable with same dimensions. So,
             ! Make sure that the patchfrac variable exists:
             ok = NF90_INQ_VARID(ncid_met,'patchfrac',id%patchfrac)
             IF(ok /= NF90_NOERR) CALL nc_abort & ! check read ok
                  (ok,'Patch-specific vegetation type (iveg) must be accompanied'// &
                  'by a patchfrac variable - this was not found in met data file '&
                  //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             DO i = 1, mland
                ! Then, get the patch specific iveg data:
                ok= NF90_GET_VAR(ncid_met, id%iveg, &
                     vegtype_metfile(i,:),&
                     start=(/i,1/), count=(/1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading iveg in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       END IF
    END IF

    ! Look for soil type:
    ok = NF90_INQ_VARID(ncid_met,'isoil',id%isoil)
    IF(ok == NF90_NOERR) THEN ! If inquiry is okay
       ! Note existence of at least one model parameter in the met file:
       exists%parameters = .TRUE.
       ! Check dimension of soil type:
       ok=NF90_INQUIRE_VARIABLE(ncid_met,id%isoil,ndims=isoil_dims)
       ! Allocate space for user-defined soil type variable:
       ALLOCATE(soiltype_metfile(mland,nmetpatches))
       ! Get soil type from met file:
       IF(metGrid=='mask') THEN
          IF(isoil_dims==2) THEN ! i.e. no patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil,data2i, &
                     start=(/land_x(i),land_y(i)/),count=(/1,1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all soil patches in grid cell to be this single type
                soiltype_metfile(i,:)=data2i(1,1)
             END DO
          ELSE IF(isoil_dims==3) THEN ! i.e. patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil, &
                     soiltype_metfile(i,:), &
                     start=(/land_x(i),land_y(i),1/),count=(/1,1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       ELSE IF(metGrid=='land') THEN
          IF(isoil_dims==1) THEN ! i.e. no patch specific isoil information
             ! Collect data from land only grid in netcdf file:
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met,id%isoil,data1i, &
                     start=(/i/),count=(/1/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
                ! Set all veg patches in grid cell to be this single type
                soiltype_metfile(i,:) = data1i(1)
             END DO
          ELSE IF(isoil_dims==2) THEN ! i.e. patch specific isoil information
             DO i = 1, mland
                ok= NF90_GET_VAR(ncid_met, id%isoil, &
                     soiltype_metfile(i,:), &
                     start=(/i,1/), count=(/1,nmetpatches/))
                IF(ok /= NF90_NOERR) CALL nc_abort &
                     (ok,'Error reading isoil in met data file ' &
                     //TRIM(filename%met)//' (SUBROUTINE open_met_file)')
             END DO
          END IF
       END IF
    END IF
    ! Deallocate read in arrays:
    IF(ASSOCIATED(temparray1)) DEALLOCATE(temparray1)
    IF(ASSOCIATED(temparray2)) DEALLOCATE(temparray2)
    
    ! Report finding met variables to log file:
    IF(all_met) THEN
       WRITE(logn,*) 'Found all met variables in met file.'
    ELSE
       WRITE(logn,*) 'Found all ESSENTIAL met variables in met file,', &
            ' some synthesised (as above).'
    END IF
    !!=================^^ End met variables search^^=======================
  END SUBROUTINE open_met_file
  !========================================================================
  SUBROUTINE get_met_data(spinup,spinConv,met,soil,rad,veg,kend,dels) 
    ! Fetches meteorological forcing data from the netcdf met forcing file
    ! for a single time step, including LAI if it exists.
    ! Note that currently met forcing is duplicated for every 
    ! vegetated patch in a single gridcell.
    ! Precision changes from REAL(4) to r_1 enable running with -r8
    USE cable_common_module, ONLY: ktau_gl
    LOGICAL, INTENT(IN) :: spinup ! are we performing a spinup?
    LOGICAL, INTENT(IN) :: spinConv ! has model spinup converged?
    TYPE(met_type),INTENT(OUT):: met ! meteorological data
    TYPE (soil_parameter_type),INTENT(IN) :: soil 
    TYPE (radiation_type),INTENT(IN) :: rad
    TYPE(veg_parameter_type),INTENT(INOUT) :: veg ! LAI retrieved from file
    INTEGER(i_d), INTENT(IN) :: kend ! total number of timesteps in run
    REAL(r_1),INTENT(IN) :: dels ! time step size
    REAL(KIND=4),DIMENSION(1,1,1) :: data3 ! temp variable for netcdf reading
    REAL(KIND=4),DIMENSION(1,1,1,1) :: data4 !  " " "
    REAL(KIND=4),DIMENSION(1,1)    :: data2 ! " "
    REAL(KIND=4),DIMENSION(1)    :: data1 ! " "
    INTEGER(i_d) :: i,j ! do loop counter

    DO i=1,mland ! over all land points/grid cells
       ! First set timing variables:
       ! All timing details below are initially written to the first patch
       ! of each gridcell, then dumped to all patches for the gridcell.
       IF(ktau_gl==1) THEN ! initialise...
          SELECT CASE(time_coord)
          CASE('LOC')! i.e. use local time by default
             ! hour-of-day = starting hod 
             met%hod(landpt(i)%cstart) = shod 
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE('GMT')! use GMT
             ! hour-of-day = starting hod + offset from GMT time:
             met%hod(landpt(i)%cstart) = shod + (longitude(i)/180.0)*12.0
             ! Note above that all met%* vars have dim mp,
             ! while longitude and latitude have dimension mland.
             met%doy(landpt(i)%cstart) = sdoy
             met%moy(landpt(i)%cstart) = smoy
             met%year(landpt(i)%cstart) = syear
          CASE DEFAULT
             CALL abort('Unknown time coordinate! ' &
                  //' (SUBROUTINE get_met_data)')
          END SELECT
       ELSE
          ! increment hour-of-day by time step size:
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + dels/3600.0
       END IF
       ! 
       IF(met%hod(landpt(i)%cstart)<0.0) THEN ! may be -ve since longitude
          ! has range [-180,180]
          ! Reduce day-of-year by one and ammend hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) - 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) + 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                met%doy(landpt(i)%cstart) = 365 ! prev year not leap year as this is
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1
             CASE(60) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(91) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(121)
                met%moy(landpt(i)%cstart) = 4
             CASE(152)
                met%moy(landpt(i)%cstart) = 5
             CASE(182)
                met%moy(landpt(i)%cstart) = 6
             CASE(213)
                met%moy(landpt(i)%cstart) = 7
             CASE(244)
                met%moy(landpt(i)%cstart) = 8
             CASE(274)
                met%moy(landpt(i)%cstart) = 9
             CASE(305)
                met%moy(landpt(i)%cstart) = 10
             CASE(335)
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          ELSE ! not a leap year or not using leap year timing
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(0) ! ie Dec previous year
                met%moy(landpt(i)%cstart) = 12
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) - 1
                ! If previous year is a leap year
                IF((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
                     (MOD(syear,4)==0.AND.MOD(syear,400)==0)) THEN
                   met%doy(landpt(i)%cstart) = 366
                ELSE
                   met%doy(landpt(i)%cstart) = 365
                END IF
             CASE(31) ! Jan
                met%moy(landpt(i)%cstart) = 1 
             CASE(59) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(90)
                met%moy(landpt(i)%cstart) = 3
             CASE(120)
                met%moy(landpt(i)%cstart) = 4
             CASE(151)
                met%moy(landpt(i)%cstart) = 5
             CASE(181)
                met%moy(landpt(i)%cstart) = 6
             CASE(212)
                met%moy(landpt(i)%cstart) = 7
             CASE(243)
                met%moy(landpt(i)%cstart) = 8
             CASE(273)
                met%moy(landpt(i)%cstart) = 9
             CASE(304)
                met%moy(landpt(i)%cstart) = 10
             CASE(334) 
                met%moy(landpt(i)%cstart) = 11
             END SELECT
          END IF ! if leap year or not
       ELSE IF(met%hod(landpt(i)%cstart)>=24.0) THEN
          ! increment or GMT adj has shifted day
          ! Adjust day-of-year and hour-of-day:
          met%doy(landpt(i)%cstart) = met%doy(landpt(i)%cstart) + 1
          met%hod(landpt(i)%cstart) = met%hod(landpt(i)%cstart) - 24.0
          ! If a leap year AND we're using leap year timing:
          IF(((MOD(syear,4)==0.AND.MOD(syear,100)/=0).OR. & 
               (MOD(syear,4)==0.AND.MOD(syear,400)==0)).AND.leaps) THEN
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(61) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(92)
                met%moy(landpt(i)%cstart) = 4
             CASE(122)
                met%moy(landpt(i)%cstart) = 5
             CASE(153)
                met%moy(landpt(i)%cstart) = 6
             CASE(183)
                met%moy(landpt(i)%cstart) = 7
             CASE(214)
                met%moy(landpt(i)%cstart) = 8
             CASE(245)
                met%moy(landpt(i)%cstart) = 9
             CASE(275)
                met%moy(landpt(i)%cstart) = 10
             CASE(306)
                met%moy(landpt(i)%cstart) = 11
             CASE(336) 
                met%moy(landpt(i)%cstart) = 12
             CASE(367)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1 
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
             ! ELSE IF not leap year and Dec 31st, increment year
          ELSE 
             SELECT CASE(INT(met%doy(landpt(i)%cstart)))
             CASE(32) ! Feb
                met%moy(landpt(i)%cstart) = 2
             CASE(60) ! Mar
                met%moy(landpt(i)%cstart) = 3
             CASE(91)
                met%moy(landpt(i)%cstart) = 4
             CASE(121)
                met%moy(landpt(i)%cstart) = 5
             CASE(152)
                met%moy(landpt(i)%cstart) = 6
             CASE(182)
                met%moy(landpt(i)%cstart) = 7
             CASE(213)
                met%moy(landpt(i)%cstart) = 8
             CASE(244)
                met%moy(landpt(i)%cstart) = 9
             CASE(274)
                met%moy(landpt(i)%cstart) = 10
             CASE(305)
                met%moy(landpt(i)%cstart) = 11
             CASE(335) 
                met%moy(landpt(i)%cstart) = 12
             CASE(366)! end of year; increment
                met%year(landpt(i)%cstart) = met%year(landpt(i)%cstart) + 1 
                met%moy(landpt(i)%cstart) = 1
                met%doy(landpt(i)%cstart) = 1
             END SELECT
          END IF ! if leap year or not
       END IF ! if increment has pushed hod to a different day
       ! Now copy these values to all veg/soil patches in the current grid cell:
       met%hod(landpt(i)%cstart:landpt(i)%cend) = met%hod(landpt(i)%cstart)
       met%doy(landpt(i)%cstart:landpt(i)%cend) = met%doy(landpt(i)%cstart)
       met%moy(landpt(i)%cstart:landpt(i)%cend) = met%moy(landpt(i)%cstart)
       met%year(landpt(i)%cstart:landpt(i)%cend) = met%year(landpt(i)%cstart)

       IF(metGrid=='mask') THEN
          ! Get SWdown data for mask grid:
          ok= NF90_GET_VAR(ncid_met,id%SWdown,data3, &
               start=(/land_x(i),land_y(i),ktau_gl/),count=(/1,1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading SWdown in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable (no units change required):
          !jhan:quick fix, use (1/2) dimension 1 here arbitrarily
          met%fsd(landpt(i)%cstart:landpt(i)%cend,1) = 0.5 * REAL(data3(1,1,1),r_1)
          met%fsd(landpt(i)%cstart:landpt(i)%cend,2) = 0.5 * REAL(data3(1,1,1),r_1)
          ! Get Tair data for mask grid:- - - - - - - - - - - - - - - - - -
          ok= NF90_GET_VAR(ncid_met,id%Tair,data4, &
               start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Tair in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable with units change:
          met%tk(landpt(i)%cstart:landpt(i)%cend) = REAL(data4(1,1,1,1),r_1) + convert%Tair
!          met%tc(landpt(i)%cstart:landpt(i)%cend) = met%tk(landpt(i)%cstart)-tfrz
          ! Get PSurf data for mask grid:- - - - - - - - - - - - - - - - - -
          IF(exists%PSurf) THEN ! IF PSurf is in met file:
             ok= NF90_GET_VAR(ncid_met,id%PSurf,data4, &
                  start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading PSurf in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%pmb(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(data4(1,1,1,1),r_1) * convert%PSurf
          ELSE ! PSurf must be fixed as a function of site elevation and T:
             met%pmb(landpt(i)%cstart:landpt(i)%cend)=1013.25* &
                  (met%tk(landpt(i)%cstart)/(met%tk(landpt(i)%cstart) + 0.0065* &
                  elevation(i)))**(9.80665/287.04/0.0065)
          END IF
          ! Get Qair data for mask grid: - - - - - - - - - - - - - - - - - -
          ok= NF90_GET_VAR(ncid_met,id%Qair,data4, &
               start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Qair in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          IF(convert%Qair==-999.0) THEN
             ! Convert relative value using only first veg/soil patch values
             ! (identical)
             CALL rh_sh(REAL(data4(1,1,1,1),r_1), met%tk(landpt(i)%cstart), &
                  met%pmb(landpt(i)%cstart),met%qv(landpt(i)%cstart))
          ELSE
             met%qv(landpt(i)%cstart:landpt(i)%cend) = REAL(data4(1,1,1,1),r_1)
          END IF
          ! Get Wind data for mask grid: - - - - - - - - - - - - - - - - - -
          IF(exists%Wind) THEN ! Scalar Wind
             ok= NF90_GET_VAR(ncid_met,id%Wind,data4, &
                  start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             ! Assign value to met data variable (no units change required):
             met%ua(landpt(i)%cstart:landpt(i)%cend) = REAL(data4(1,1,1,1),r_1)
          ELSE ! Vector wind
             ! Get Wind_N:
             ok= NF90_GET_VAR(ncid_met,id%Wind,data4, &
                  start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind_N in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%ua(landpt(i)%cstart) = REAL(data4(1,1,1,1),r_1) ! only part of wind variable
             ok= NF90_GET_VAR(ncid_met,id%Wind_E,data4, &
                  start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind_E in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             ! Write final scalar Wind value:
             met%ua(landpt(i)%cstart:landpt(i)%cend) &
                  = SQRT(met%ua(landpt(i)%cstart)**2 + REAL(data4(1,1,1,1),r_1)**2)
          END IF
          ! Get Rainf and Snowf data for mask grid:- - - - - - - - - - - - -
          ok= NF90_GET_VAR(ncid_met,id%Rainf,data3, &
               start=(/land_x(i),land_y(i),ktau_gl/),count=(/1,1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Rainf in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          met%precip(landpt(i)%cstart:landpt(i)%cend) = REAL(data3(1,1,1),r_1) ! store Rainf
          IF(exists%Snowf) THEN
             ok= NF90_GET_VAR(ncid_met,id%Snowf,data3, &
                  start=(/land_x(i),land_y(i),ktau_gl/),count=(/1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Snowf in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             ! store Snowf value (EK nov2007)
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = REAL(data3(1,1,1),r_1)
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart:landpt(i)%cend) &
                  + met%precip_sn(landpt(i)%cstart:landpt(i)%cend)
          END IF
          ! Convert units:
          met%precip(landpt(i)%cstart:landpt(i)%cend) = &
               met%precip(landpt(i)%cstart) * convert%Rainf
          met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
               met%precip_sn(landpt(i)%cstart) * convert%Rainf  ! (EK nov2007)
          ! If we're performing a spinup, the spinup hasn't converged, 
          ! and an avPrecip variable has been found, modify precip to 
          ! ensure reasonable equilibration:
          IF(spinup.AND.(.NOT.spinConv).AND.exists%avPrecip) THEN
             ! Rescale precip to average rainfall for this site:
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart:landpt(i)%cend) / PrecipScale(i)
             ! Added for snow (EK nov2007)
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip_sn(landpt(i)%cstart:landpt(i)%cend) / PrecipScale(i)
          END IF
          ! Get LWdown data for mask grid: - - - - - - - - - - - - - - - - - 
          IF(exists%LWdown) THEN ! If LWdown exists in met file
             ok= NF90_GET_VAR(ncid_met,id%LWdown,data3, &
                  start=(/land_x(i),land_y(i),ktau_gl/),count=(/1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading LWdown in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%fld(landpt(i)%cstart:landpt(i)%cend)=REAL(data3(1,1,1),r_1)
          ELSE ! Synthesise LWdown based on temperature
             ! Use Swinbank formula:
             met%fld(landpt(i)%cstart:landpt(i)%cend) = &
                  0.0000094*0.0000000567*(met%tk(landpt(i)%cstart)**6.0)
          END IF
          ! Get CO2air data for mask grid:- - - - - - - - - - - - - - - - - -
          IF(exists%CO2air) THEN ! If CO2air exists in met file
             ok= NF90_GET_VAR(ncid_met,id%CO2air,data4, &
                  start=(/land_x(i),land_y(i),1,ktau_gl/),count=(/1,1,1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading CO2air in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%ca(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(data4(1,1,1,1),r_1)/1000000.0
          ELSE 
             ! Fix CO2 air concentration:
             met%ca(landpt(i)%cstart:landpt(i)%cend) = fixedCO2 /1000000.0
          END IF
          ! Get LAI, if it's present, for mask grid:- - - - - - - - - - - - -
          IF(exists%LAI) THEN ! If LAI exists in met file
             IF(exists%LAI_T) THEN ! i.e. time dependent LAI
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1,nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data4, &
                           start=(/land_x(i),land_y(i),j,ktau_gl/),count=(/1,1,1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met1 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      !                veg%vlai(landpt(i)%cstart:landpt(i)%cend+j-1) = data4(1,1,1,1)
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data4(1,1,1,1),r_1)     ! BP apr08
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                        start=(/land_x(i),land_y(i),ktau_gl/),count=(/1,1,1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met2 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = REAL(data3(1,1,1),r_1)
                END IF
             ELSEIF(exists%LAI_M) THEN ! i.e. monthly LAI (BP apr08)
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1,nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data4, &
                           start=(/land_x(i),land_y(i),j,met%moy/),count=(/1,1,1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met3 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data4(1,1,1,1),r_1)
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                        start=(/land_x(i),land_y(i),met%moy/),count=(/1,1,1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met4 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = REAL(data3(1,1,1),r_1)
                END IF
             ELSE ! i.e. time independent LAI
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1,nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                           start=(/land_x(i),land_y(i),j/),count=(/1,1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met5 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data3(1,1,1),r_1)
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data2, &
                        start=(/land_x(i),land_y(i)/),count=(/1,1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met6 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = REAL(data2(1,1),r_1)
                END IF
             END IF
          ELSE 
             ! If not in met file, use default LAI value:
             veg%vlai(landpt(i)%cstart:landpt(i)%cend) =  &
                  defaultLAI(i,met%moy(landpt(i)%cstart))
          END IF

       ELSE IF(metGrid=='land') THEN
          ! Collect data from land only grid in netcdf file:
          ! Get SWdown data for land-only grid: - - - - - - - - - - - - -
          IF (ncciy > 0) ncid_met = ncid_sw
          ok= NF90_GET_VAR(ncid_met,id%SWdown,data2, &
               start=(/i,ktau_gl/),count=(/1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading SWdown in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable (no units change required):
          !jhan:quick fix, use dimension 1 here arbitrarily
          met%fsd(landpt(i)%cstart:landpt(i)%cend,1) = 0.5 *  REAL(data2(1,1),r_1)
          met%fsd(landpt(i)%cstart:landpt(i)%cend,2) = 0.5 * REAL(data2(1,1),r_1)
          ! Get Tair data for land-only grid:- - - - - - - - - - - - - - -
          IF (ncciy > 0) ncid_met = ncid_ta
          ok= NF90_GET_VAR(ncid_met,id%Tair,data2, &
               start=(/i,ktau_gl/),count=(/1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Tair in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          ! Assign value to met data variable with units change:
          met%tk(landpt(i)%cstart:landpt(i)%cend) = &
               REAL(data2(1,1),r_1) + convert%Tair
!          met%tc(landpt(i)%cstart:landpt(i)%cend) =  &
!               met%tk(landpt(i)%cstart)-tfrz
          ! Get PSurf data for land-only grid:- -- - - - - - - - - - - - - -
          IF(exists%PSurf) THEN ! IF PSurf is in met file:
             IF (ncciy > 0) ncid_met = ncid_ps
             ok= NF90_GET_VAR(ncid_met,id%PSurf,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
!hzz to fix the problem of ps data on time step 2184, (no pressure for 120-180E)
!    using previous day ps instead
             IF ((ncciy == 1986).and.(ktau_gl == 2184)) then
             ok= NF90_GET_VAR(ncid_met,id%PSurf,data2, &
                  start=(/i,2176/),count=(/1,1/))
             endif
!hzz end
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading PSurf in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%pmb(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(data2(1,1),r_1) * convert%PSurf
          ELSE ! PSurf must be fixed as a function of site elevation and T:
             met%pmb(landpt(i)%cstart:landpt(i)%cend) = 1013.25 &
                  *(met%tk(landpt(i)%cstart)/(met%tk(landpt(i)%cstart) &
                  + 0.0065*elevation(i)))**(9.80665/287.04/0.0065)
          END IF
          ! Get Qair data for land-only grid:- - - - - - - - - - - - - - - -
          IF (ncciy > 0) ncid_met = ncid_qa
          ok= NF90_GET_VAR(ncid_met,id%Qair,data2, &
               start=(/i,ktau_gl/),count=(/1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Qair in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          IF(convert%Qair==-999.0) THEN
             CALL rh_sh(REAL(data2(1,1),r_1), met%tk(landpt(i)%cstart), &
                  met%pmb(landpt(i)%cstart),met%qv(landpt(i)%cstart))
             met%qv(landpt(i)%cstart:landpt(i)%cend)=met%qv(landpt(i)%cstart)
          ELSE
             met%qv(landpt(i)%cstart:landpt(i)%cend) = REAL(data2(1,1),r_1)
          END IF
          ! Get Wind data for land-only grid: - - - - - - - - - - - - - - - -
          IF (ncciy > 0) ncid_met = ncid_wd
          IF(exists%Wind) THEN ! Scalar Wind
             ok= NF90_GET_VAR(ncid_met,id%Wind,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             ! Assign value to met data variable (no units change required):
             met%ua(landpt(i)%cstart:landpt(i)%cend) = REAL(data2(1,1),r_1)
          ELSE ! Vector wind
             ! Get Wind_N:
             ok= NF90_GET_VAR(ncid_met,id%Wind,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind_N in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%ua(landpt(i)%cstart) = REAL(data2(1,1),r_1) ! only part of the wind variable
             ok= NF90_GET_VAR(ncid_met,id%Wind_E,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Wind_E in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             ! Write final scalar Wind value:
             met%ua(landpt(i)%cstart:landpt(i)%cend) = &
                  SQRT(met%ua(landpt(i)%cstart)**2 + REAL(data2(1,1),r_1)**2)
          END IF
          ! Get Rainf and Snowf data for land-only grid: - - - - - - - - - - -
          IF (ncciy > 0) ncid_met = ncid_rain
          ok= NF90_GET_VAR(ncid_met,id%Rainf,data2, &
               start=(/i,ktau_gl/),count=(/1,1/))
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error reading Rainf in met data file ' &
               //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
          met%precip(landpt(i)%cstart:landpt(i)%cend) = REAL(data2(1,1),r_1) ! store Rainf
          IF (ncciy > 0) ncid_met = ncid_snow
          IF(exists%Snowf) THEN
             ok= NF90_GET_VAR(ncid_met,id%Snowf,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading Snowf in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = REAL(data2(1,1),r_1)
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart:landpt(i)%cend) &
                  + met%precip_sn(landpt(i)%cstart:landpt(i)%cend)
          END IF
          ! Convert units:
          met%precip(landpt(i)%cstart:landpt(i)%cend) = &
               met%precip(landpt(i)%cstart:landpt(i)%cend) * convert%Rainf
          met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = &
               met%precip_sn(landpt(i)%cstart:landpt(i)%cend) * convert%Rainf
          ! If we're performing a spinup, the spinup hasn't converged, 
          ! and an avPrecip variable has been found, modify precip to 
          ! ensure reasonable equilibration:
          IF(spinup.AND.(.NOT.spinConv).AND.exists%avPrecip) THEN
             ! Rescale precip to average rainfall for this site:
             met%precip(landpt(i)%cstart:landpt(i)%cend) = &
                  met%precip(landpt(i)%cstart) / PrecipScale(i)
          END IF
          ! Get LWdown data for land-only grid: - - - - - - - - - - - - - - 
          IF (ncciy > 0) ncid_met = ncid_lw
          IF(exists%LWdown) THEN ! If LWdown exists in met file
             ok= NF90_GET_VAR(ncid_met,id%LWdown,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading LWdown in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%fld(landpt(i)%cstart:landpt(i)%cend)=REAL(data2(1,1),r_1)
          ELSE ! Synthesise LWdown based on temperature
             ! Use Swinbank formula:
             met%fld(landpt(i)%cstart:landpt(i)%cend) = &
                  0.0000094*0.0000000567*(met%tk(landpt(i)%cstart)**6.0)
          END IF
          ! Get CO2air data for land-only grid:- - - - - - - - - - - - - -
          IF(exists%CO2air) THEN ! If CO2air exists in met file
             ok= NF90_GET_VAR(ncid_met,id%CO2air,data2, &
                  start=(/i,ktau_gl/),count=(/1,1/))
             IF(ok /= NF90_NOERR) CALL nc_abort &
                  (ok,'Error reading CO2air in met data file ' &
                  //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
             met%ca(landpt(i)%cstart:landpt(i)%cend) = &
                  REAL(data2(1,1),r_1)/1000000.0
          ELSE 
             ! Fix CO2 air concentration:
             met%ca(landpt(i)%cstart:landpt(i)%cend) = fixedCO2 /1000000.0
          END IF
          ! Get LAI data, if it exists, for land-only grid:- - - - - - - - -
          IF(exists%LAI) THEN ! If LAI exists in met file
             IF(exists%LAI_T) THEN ! i.e. time dependent LAI
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1, nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                           start=(/i,j,ktau_gl/),count=(/1,1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met7 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data3(1,1,1),r_1)
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data2, &
                        start=(/i,ktau_gl/),count=(/1,1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met8 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(data2(1,1),r_1)
                END IF
             ELSEIF(exists%LAI_M) THEN ! i.e. monthly LAI (BP apr08)
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1, nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data3, &
                           start=(/i,j,met%moy/),count=(/1,1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met9 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data3(1,1,1),r_1)
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data2, &
                        start=(/i,met%moy/),count=(/1,1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met10 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(data2(1,1),r_1)
                END IF
             ELSE ! LAI time independent
                IF(exists%LAI_P) THEN ! i.e. patch dependent LAI
                   DO j=1, nmetpatches
                      ok= NF90_GET_VAR(ncid_met,id%LAI,data2, &
                           start=(/i,j/),count=(/1,1/))
                      IF(ok /= NF90_NOERR) CALL nc_abort &
                           (ok,'Error reading LAI in met11 data file ' &
                           //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                      veg%vlai(landpt(i)%cstart+j-1) = REAL(data2(1,1),r_1)
                   END DO
                ELSE ! i.e. patch independent LAI
                   ok= NF90_GET_VAR(ncid_met,id%LAI,data1, &
                        start=(/i/),count=(/1/))
                   IF(ok /= NF90_NOERR) CALL nc_abort &
                        (ok,'Error reading LAI in met12 data file ' &
                        //TRIM(filename%met)//' (SUBROUTINE get_met_data)')
                   veg%vlai(landpt(i)%cstart:landpt(i)%cend) = &
                        REAL(data1(1),r_1)
                END IF
             END IF
          ELSE 
             ! If not in met file, use default LAI value:
             veg%vlai(landpt(i)%cstart:landpt(i)%cend) =  &
                  defaultLAI(i,met%moy(landpt(i)%cstart))
          END IF
       ELSE
          CALL abort('Unrecognised grid type')
       END IF ! grid type

       ! Set solid precip based on temp
       met%precip_sn(landpt(i)%cstart:landpt(i)%cend) = 0.0 ! (EK nov2007)
       IF( met%tk(landpt(i)%cstart) <= tfrz ) &
            met%precip_sn(landpt(i)%cstart:landpt(i)%cend) &
            = met%precip(landpt(i)%cstart) ! (EK nov2007)

    END DO ! 1, mland over all land grid points

    ! Set cosine of zenith angle (provided by GCM when online):
    met%coszen = sinbet(met%doy, rad%latitude, met%hod)
    ! initialise within canopy air temp
    met%tvair = met%tk 
    met%tvrad = met%tk 
    IF(check%ranges) THEN
       ! Check ranges are okay:
          !jhan:quick fix, use dimension 1 here arbitrarily
       IF(ANY(met%fsd(:,1)<ranges%SWdown(1)).OR.ANY(met%fsd(:,1)>ranges%SWdown(2))) &
            CALL abort('SWdown out of specified ranges!')
       IF(ANY(met%fsd(:,2)<ranges%SWdown(1)).OR.ANY(met%fsd(:,2)>ranges%SWdown(2))) &
            CALL abort('SWdown out of specified ranges!')
       IF(ANY(met%fld<ranges%LWdown(1)).OR.ANY(met%fld>ranges%LWdown(2))) &
            CALL abort('LWdown out of specified ranges!')
       IF(ANY(met%qv<ranges%Qair(1)).OR.ANY(met%qv>ranges%Qair(2))) &
            CALL abort('Qair out of specified ranges!')
       IF(ANY(met%precip<ranges%Rainf(1)).OR.ANY(met%precip>ranges%Rainf(2))) &
            CALL abort('Rainf out of specified ranges!')
       IF(ANY(met%ua<ranges%Wind(1)).OR.ANY(met%ua>ranges%Wind(2))) &
            CALL abort('Wind out of specified ranges!')
       IF(ANY(met%tk<ranges%Tair(1)).OR.ANY(met%tk>ranges%Tair(2))) &
            CALL abort('Tair out of specified ranges!')
       IF(ANY(met%pmb<ranges%PSurf(1)).OR.ANY(met%pmb>ranges%PSurf(2))) &
            CALL abort('PSurf out of specified ranges!')
    END IF

  END SUBROUTINE get_met_data
  !============================================================================
  SUBROUTINE close_met_file
    ok=NF90_CLOSE(ncid_met)
    IF(ok /= NF90_NOERR) CALL nc_abort (ok,'Error closing met data file ' &
         //TRIM(filename%met)//' (SUBROUTINE close_met_file)')
    ! Clear lat_all and lon_all variables
    DEALLOCATE(lat_all,lon_all)
  END SUBROUTINE close_met_file
  !============================================================================
  SUBROUTINE load_parameters(met,air,ssoil,veg,bgc,soil,canopy, &
         rough,rad,sum_flux,bal,logn,vegparmnew, &
         casabiome,casapool,casaflux,casamet,casabal,phen)
  ! Checks where parameters and initialisations should be loaded from.
  ! If they can be found in either the met file or restart file, they will 
  ! load from there, with the met file taking precedence. Otherwise, they'll
  ! be chosen from a coarse global grid of veg and soil types, based on 
  ! the lat/lon coordinates.
  !
  ! Input variables not listed:
  !   mland          - via define_dimensions
  !   filename%type  - via io_variables
  !   exists%type    - via io_variables
  !   smoy           - via io_variables
  ! Output variables not listed:
  !   (determined here or from sub get_default_params <- countPatch)
  !   mp             - via define_dimensions
  !   landpt%type    - via io_variables (nap,cstart,cend,ilon,ilat)
  !   max_vegpatches - via io_variables

    IMPLICIT NONE
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    TYPE (casa_biome)  , INTENT(OUT) :: casabiome
    TYPE (casa_pool)   , INTENT(OUT) :: casapool
    TYPE (casa_flux)   , INTENT(OUT) :: casaflux
    TYPE (casa_met)    , INTENT(OUT) :: casamet
    TYPE (casa_balance), INTENT(OUT) :: casabal
    TYPE(phen_variable), INTENT(OUT) :: phen
    INTEGER,INTENT(IN) :: logn     ! log file unit number
    LOGICAL,INTENT(IN) :: vegparmnew  ! are we using the new format?
    REAL(r_1),POINTER,DIMENSION(:) :: pfractmp ! temp store of patch fraction
    LOGICAL :: completeSet ! was a complete parameter set found?
    INTEGER(i_d) :: mp_restart ! total number of patches in restart file
    INTEGER(i_d) :: mpID,napID
    INTEGER(i_d) :: i ! do loop variables

    ! Allocate spatial heterogeneity variables:
    ALLOCATE(landpt(mland))

    WRITE(logn,*) '-------------------------------------------------------'
    WRITE(logn,*) 'Looking for parameters and initial states....'  
    WRITE(logn,*) ' Loading initialisations from default grid.'

    ! Parameter values and some grid info are read in.
    ! They will be overwritten by values from the restart file, if present.
    ! Those variables found in the met file will again overwrite existing ones.

    CALL get_default_params(logn,vegparmnew)

    CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
            sum_flux,veg,mp)
    WRITE(logn,*) ' CABLE variables allocated with ', mp, ' patch(es).'
    IF (icycle > 0) THEN
      CALL alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
      CALL alloc_phenvariable(phen,mp)
    ENDIF

    ! Write parameter values to CABLE's parameter variables:
    CALL write_default_params(met,air,ssoil,veg,bgc,soil,canopy,rough, &
            rad,logn,vegparmnew,smoy)
    IF (icycle > 0) THEN
      CALL write_cnp_params(veg,casaflux,casamet)
      CALL casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
      CALL casa_readphen(veg,casamet,phen)
      CALL casa_init(casabiome,casamet,casapool,casabal,veg,phen)
    ENDIF

! removed get_default_inits and get_default_lai as they are already done
! in write_default_params
!    ! Load default initialisations from Mk3L climatology:
!    CALL get_default_inits(met,soil,ssoil,canopy,logn)
!
!    ! load default LAI values from global data:
!    CALL get_default_lai

    ! Look for restart file (which will have parameters):
    ok = NF90_OPEN(filename%restart_in,0,ncid_rin) ! open restart file
    IF (ok /= NF90_NOERR) THEN
      ! With no restart file, use default parameters already loaded
      WRITE(logn,*) ' Could not find restart file ', TRIM(filename%restart_in)
      WRITE(logn,*) ' Pre-loaded default initialisations are used.'
      WRITE(*,*)    ' Could not find restart file ', TRIM(filename%restart_in)
      WRITE(*,*)    ' Pre-loaded default initialisations are used.'

    ELSE
      ! Restart file exists, parameters and init will be loaded from it.
      WRITE(logn,*) ' Overwriting initialisations with values in ', &
                    'restart file: ', TRIM(filename%restart_in)
      WRITE(*,*)    ' Overwriting initialisations with values in ', &
                    'restart file: ', TRIM(filename%restart_in)
     
      ! Check total number of patches in restart file:
      ok = NF90_INQ_DIMID(ncid_rin,'mp',mpID)
      IF(ok /= NF90_NOERR) THEN
        ok = NF90_INQ_DIMID(ncid_rin,'mp_patch',mpID)
        IF(ok /= NF90_NOERR)  CALL nc_abort &
           (ok,'Error finding mp or mp_patch dimension in restart file ' &
           //TRIM(filename%restart_in)//' (SUBROUTINE load_parameters) ' &
           //'Recommend running without restart file.')
      END IF
      ok = NF90_INQUIRE_DIMENSION(ncid_rin,mpID,len=mp_restart)
      IF(ok /= NF90_NOERR) CALL nc_abort &
           (ok,'Error finding total number of patches in restart file ' &
           //TRIM(filename%restart_in)//' (SUBROUTINE load_parameters) ' &
           //'Recommend running without restart file.')
      ! Check that mp_restart = mp from default/met values
      IF(mp_restart /= mp) CALL abort('Number of patches in '// &
           'restart file '//TRIM(filename%restart_in)//' does not equal '// &
           'to number in default/met file settings. (SUB load_parameters) ' &
           //'Recommend running without restart file.')

      ! Load initialisations and parameters from restart file:
      CALL get_restart_data(logn,ssoil,canopy,rough,bgc,bal,veg, &
                            soil,rad,vegparmnew)

    END IF ! if restart file exists

    ! Overwrite default values by those available in met file:
    CALL get_parameters_met(soil,veg,bgc,rough,completeSet)

    ! Results of looking for parameters in the met file:
    WRITE(logn,*)
    IF(exists%parameters.AND.completeSet) THEN
      ! All pars were found in met file:
      WRITE(logn,*) ' Loaded all parameters from met input file: ', &
                    TRIM(filename%met)
    ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
      ! Only some pars were found in met file:
      WRITE(logn,*) ' Loaded some parameters from met input file: ', &
                    TRIM(filename%met), & ! write to log file
                    ' the rest are default values'
      WRITE(*,*)    ' Loaded some parameters from met input file: ', &
                    TRIM(filename%met), & ! write to screen
                    ' the rest are default values - check log file'
    ELSE
      ! No parameters were found in met file:
      WRITE(logn,*) ' Met file has no params; all parameters remain as default.'
    END IF
    WRITE(logn,*)

    ! Construct derived parameters and zero initialisations, regardless 
    ! of where parameters and other initialisations have loaded from:
    CALL derived_parameters(soil,sum_flux,bal,ssoil,veg,rough)

    ! Check for basic inconsistencies in parameter values:
    CALL check_parameter_values(soil,veg,ssoil)

    ! Write per-site parameter values to log file if requested:
    CALL report_parameters(logn,soil,veg,bgc,rough,ssoil,canopy, &
         vegparmnew,verbose)

  END SUBROUTINE load_parameters

  !============================================================================

  SUBROUTINE old_load_parameters(met,air,ssoil,veg,bgc,soil,canopy, &
       & rough,rad,sum_flux,bal,logn,vegparmnew)
  ! Checks where parameters and initialisations should be loaded from.
  ! If they can be found in either the met file or restart file, they will 
  ! load from there, with the met file taking precedence. Otherwise, they'll
  ! be chosen from a coarse global grid of veg and soil types, based on 
  ! the lat/lon coordinates.
  !
  ! Input variables not listed:
  !   mland          - via define_dimensions
  !   filename%type  - via io_variables
  !   exists%type    - via io_variables
  !   smoy           - via io_variables
  ! Output variables not listed:
  !   (determined here or from sub get_default_params <- countPatch)
  !   mp             - via define_dimensions
  !   landpt%type    - via io_variables (nap,cstart,cend,ilon,ilat)
  !   max_vegpatches - via io_variables
   use cable_common_module, only : get_type_parameters
    IMPLICIT NONE
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(OUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(OUT)  :: veg
    TYPE (bgc_pool_type), INTENT(OUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(OUT) :: soil
    TYPE (canopy_type), INTENT(OUT)    :: canopy
    TYPE (roughness_type), INTENT(OUT) :: rough
    TYPE (radiation_type),INTENT(OUT)  :: rad
    TYPE (sum_flux_type), INTENT(OUT)  :: sum_flux
    TYPE (balances_type), INTENT(OUT)  :: bal
    INTEGER,INTENT(IN) :: logn     ! log file unit number
    LOGICAL,INTENT(IN) :: vegparmnew  ! are we using the new format?
    REAL(r_1),POINTER,DIMENSION(:) :: pfractmp ! temp store of patch fraction
    LOGICAL :: completeSet ! was a complete parameter set found?
    INTEGER(i_d) :: mp_restart ! total number of patches in restart file
    INTEGER(i_d) :: mpID,napID
    INTEGER(i_d) :: i ! do loop variables

    ! Allocate spatial heterogeneity variables:
    ALLOCATE(landpt(mland))

    WRITE(logn,*) '-------------------------------------------------------'
    WRITE(logn,*) 'Looking for parameters and initial states....'

    ! Unless a restart file is found AND the met file has LAI,
    ! default_params will be called (although the initialisations 
    ! parameter values it writes may be overwritten by restart 
    ! and/or met file values, if found).

    ! Look for restart file (which will have parameters):
    ok = NF90_OPEN(filename%restart_in,0,ncid_rin) ! open restart file
    IF (ok /= NF90_NOERR) THEN
       ! If no restart file, use default_parameters to load default
       ! parameters and initialisations (and get default grid coords 
       ! if default LAI is required). Parameter values will be over written 
       ! by any found in the met file.
       WRITE(logn,*) ' Could not find restart file ',&
            TRIM(filename%restart_in)
       WRITE(logn,*) ' Loading initialisations ', &
            'from default grid. '

       ! Load default parameters and initialisations (this will also write 
       ! to gdpt, landpt%cstart,landpt%cend)
       CALL get_default_params(logn,vegparmnew)

       ! Allocate CABLE's main variables:
       CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
            sum_flux,veg,mp)
       WRITE(logn,*) ' CABLE variables allocated with ', mp, ' patch(es).'

       ! Write parameter values to CABLE's parameter variables:
       CALL write_default_params(met,air,ssoil,veg,bgc,soil,canopy,rough, &
            rad,logn,vegparmnew,smoy)

       ! Load default initialisations from Mk3L climatology:
       CALL get_default_inits(met,soil,ssoil,canopy,logn)

       ! Check if there are any parameters in the met file; if so, 
       ! overwrite default values for those that are available:
       CALL get_parameters_met(soil,veg,bgc,rough,completeSet)

       ! Results of looking for pars in the met file:
       WRITE(logn,*)
       IF(exists%parameters.AND.completeSet) THEN
          ! All pars were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from met input file: ', TRIM(filename%met)
       ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
          ! Only some pars were found in met file:
          WRITE(logn,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename%met), & ! write to log file
               ' the rest are default values'
          WRITE(*,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename%met), & ! write to screen
               ' the rest are default values - check log file'
       ELSE
          ! No parameters were found in met file:
          WRITE(logn,*) ' Loaded all parameters from default grid.'
       END IF

       ! If LAI does not exist in the met file, load default values for run:
       IF(.NOT.exists%LAI) THEN
          CALL get_default_lai
       END IF

    ELSE ! RESTART FILE EXISTS
       ! If met file does not have LAI, then default_params will still
       ! need to be called to load default grid to get LAI in the main
       ! time step loop. The parameters and initialisations it writes 
       ! will be overwritten by the restart (and possibly met file) values.

       ! Restart file exists, parameters and init will be loaded from it.
       WRITE(logn,*) ' Loading initialisations ', &
            'from restart file: ', TRIM(filename%restart_in)  

       ! Check total number of patches in restart file:
       ok = NF90_INQ_DIMID(ncid_rin,'mp',mpID)
       IF(ok /= NF90_NOERR) THEN
         ok = NF90_INQ_DIMID(ncid_rin,'mp_patch',mpID)
         IF(ok /= NF90_NOERR)  CALL nc_abort &
            (ok,'Error finding mp or mp_patch dimension in restart file ' &
            //TRIM(filename%restart_in)//' (SUBROUTINE load_parameters)')
       END IF
       ok = NF90_INQUIRE_DIMENSION(ncid_rin,mpID,len=mp_restart)
       IF(ok /= NF90_NOERR) CALL nc_abort &
            (ok,'Error finding total number of patches in restart file ' &
            //TRIM(filename%restart_in)//' (SUBROUTINE load_parameters)')

       IF(.NOT.exists%LAI) THEN ! i.e. no LAI in met file
          ! Call get_default_params to get default grid details;
          ! i.e. set 'gdpt' - this will be used to load LAI later.
          ! This also sets landpt(:)%nap, cstart and cend:
          WRITE(logn,*) ' Loading default grid to use default LAI:'
          CALL get_default_params(logn,vegparmnew)
          ! Check that mp_restart = mp from default/met values
          IF(mp_restart /= mp) CALL abort('Number of patches in '// &
               'restart file '//TRIM(filename%restart_in)//' is not equal '// &
               'to number in default/met file settings. (SUB load_parameters)')

          ! If LAI does not exist in the met file, load default values for run:
          IF(.NOT.exists%LAI) THEN
             CALL get_default_lai
          END IF

       ELSE ! i.e. LAI in met file - number of active patches per grid cell and 
          ! therefore CABLE array bounds will be decided from restart file alone

          ! get mvtype, mstype and allocate/get veg_desc
          CALL get_type_parameters(logn,vegparmnew)

          ! Set mp to be the value found in the restart file:
          mp = mp_restart

          ! Load number of active patches from restart file:
          ok = NF90_INQ_VARID(ncid_rin,'nap',napID)
          IF(ok /= NF90_NOERR) CALL nc_abort &
               (ok,'Error finding number of active patches in restart file ' &
               //TRIM(filename%restart_in)//' (SUBROUTINE load_parameters)')
          ok=NF90_GET_VAR(ncid_rin,napID,landpt(:)%nap)            
          IF(ok/=NF90_NOERR) CALL nc_abort(ok,'Error reading number of '// &
               'active patches in file ' &
               //TRIM(filename%restart_in)// '(SUBROUTINE load_parameters)')

          ! Establish landpt(:)%cstart and landpt(:)%cend
          ! (which would otherwise be established by get_default_params):
          DO i=1,mland
             IF(i==1) THEN
                landpt(i)%cstart = 1 ! cstart of the first landpt should start with 1
                landpt(i)%cend = landpt(i)%nap
             ELSE
                landpt(i)%cstart = landpt(i-1)%cend + 1
                landpt(i)%cend = landpt(i)%cstart + landpt(i)%nap - 1
             END IF
          END DO

          ! Set the maximum number of active patches in any grid cell:
          max_vegpatches = MAXVAL(landpt(:)%nap)
          WRITE(logn,*) ' max_vegpatches is now set at ', max_vegpatches

       END IF ! if LAI not in met file

       ! Allocate CABLE's main variables:
       CALL allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
            sum_flux,veg,mp)

       ! Load initialisations and parameters from restart file:
       CALL get_restart_data(logn,ssoil,canopy,rough,bgc,bal,veg, &
            soil,rad,vegparmnew)

       ! Save value of patch fractions from restart file 
       ! (potentially overwriten if patchfrac is in met file):
       ALLOCATE(pfractmp(mp))
       pfractmp = patch(:)%frac
       ! Initialise:
       exists%parameters=.FALSE.
       ! Overwrite any parameters found in met file:
       CALL get_parameters_met(soil,veg,bgc,rough,completeSet)
       ! If met file patchfrac suggests more active patches than the restart
       ! file, print a warning (restart will only contain initialisations and
       ! parameters for ACTIVE patches):
       DO i=1,mp
          IF(pfractmp(i)==0.0 .AND. patch(i)%frac/=0.0) THEN
             WRITE(logn,'(A47,A54)') &
                  'WARNING - met file contains more active patches',&
                  ' than restart file - SOME SITES MAY NOT BE INITIALISED'
             WRITE(*,'(A47,A54)') &
                  'WARNING - met file contains more active patches',&
                  ' than restart file - SOME SITES MAY NOT BE INITIALISED'
          END IF
          IF(pfractmp(i)/=0.0 .AND. patch(i)%frac==0.0) THEN
             WRITE(logn,'(A45,A51)') &
                  'ERROR - met file contains less active patches',&
                  ' than restart file - CHECK THE VALIDITY OF THE FILE'
             WRITE(*,'(A45,A51)') &
                  'ERROR - met file contains less active patches',&
                  ' than restart file - CHECK THE VALIDITY OF THE FILE'
             CALL abort('Check the patchfrac in both files.')
          END IF
       END DO
       DEALLOCATE(pfractmp)
       ! Results of looking for pars in the met file:
       WRITE(logn,*)
       IF(exists%parameters.AND.completeSet) THEN
          ! All pars were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from met input file: ', TRIM(filename%met)
       ELSE IF(exists%parameters.AND..NOT.completeSet) THEN
          ! Only some pars were found in met file:
          WRITE(logn,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename%met), &
               ' the rest are from restart file: ', &
               TRIM(filename%restart_in)
          WRITE(*,*) ' Loaded some parameters from met input file: ', &
               TRIM(filename%met), &
               ' the rest are from restart file: ', &
               TRIM(filename%restart_in),' - check log file'
       ELSE
          ! No parameters were found in met file:
          WRITE(logn,*) ' Loaded all parameters ', &
               'from restart file: ',TRIM(filename%restart_in)
       END IF
    END IF ! if restart file exists
    WRITE(logn,*)

    ! Construct derived parameters and zero initialisations, regardless 
    ! of where parameters and other initialisations have loaded from:
    CALL derived_parameters(soil,sum_flux,bal,ssoil,veg,rough)

    ! Check for basic inconsistencies in parameter values:
    CALL check_parameter_values(soil,veg,ssoil)

    ! Write per-site parameter values to log file if requested:
    CALL report_parameters(logn,soil,veg,bgc,rough,ssoil,canopy, &
         vegparmnew,verbose)

  END SUBROUTINE old_load_parameters

  !============================================================================
  SUBROUTINE get_parameters_met(soil,veg,bgc,rough,completeSet)
    ! This subroutine looks for parameters in the met file, and 
    ! loads those that are found.
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (roughness_type), INTENT(INOUT) :: rough
    LOGICAL, INTENT(OUT) :: completeSet ! were all pars found?
    INTEGER(i_d) :: parID ! parameter's netcdf ID

    ! removed the following section because already in IGBP types (BP apr08)
    !    ! First, if user defined surface type ratios are present in the 
    !    ! met file then use them:
    !    IF(ASSOCIATED(vegfrac_user)) THEN
    !       DO i=1,mland
    !          ! Overwrite landpt(i)%*%frac, which will be set either by restart
    !          ! or default values:
    !          landpt(i)%veg%frac = vegfrac_user(i)
    !          landpt(i)%urban%frac = urbanfrac_user(i)
    !          landpt(i)%lake%frac = lakefrac_user(i)
    !          landpt(i)%ice%frac = icefrac_user(i)
    !       END DO
    !    END IF

    completeSet=.TRUE. ! initialise (assume all param will load from met file)

    ! Get parameter values:
    ! Arguments: netcdf file ID; parameter name; complete set check;
    !   parameter value; filename for error messages; number of veg/soil patches
    !   in met file; switch to indicate size of dimensions of the parameter.
    ! ! Use 'defd' for single dim double precision.
! veg and soil types already obtained in sub open_met_file
!    CALL readpar(ncid_met,'iveg',completeSet,veg%iveg,filename%met, &
!         nmetpatches,'def')
    CALL readpar(ncid_met,'patchfrac',completeSet,patch(:)%frac,filename%met, &
         nmetpatches,'def')
!    CALL readpar(ncid_met,'isoil',completeSet,soil%isoilm,filename%met, &
!         nmetpatches,'def')
    CALL readpar(ncid_met,'clay',completeSet,soil%clay,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sand',completeSet,soil%sand,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'silt',completeSet,soil%silt,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'ssat',completeSet,soil%ssat,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sfc',completeSet,soil%sfc,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'swilt',completeSet,soil%swilt,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'bch',completeSet,soil%bch,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'hyds',completeSet,soil%hyds,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'sucs',completeSet,soil%sucs,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'css',completeSet,soil%css,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rhosoil',completeSet,soil%rhosoil,filename%met, &
         nmetpatches,'def')
!    CALL readpar(ncid_met,'rs20',completeSet,soil%rs20,filename%met, &
!         nmetpatches,'def')
    CALL readpar(ncid_met,'rs20',completeSet,veg%rs20,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'albsoil',completeSet,soil%albsoil,filename%met, &
         nmetpatches,'nrb')
    CALL readpar(ncid_met,'froot',completeSet,veg%froot,filename%met, &
         nmetpatches,'ms')
    CALL readpar(ncid_met,'hc',completeSet,veg%hc,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'canst1',completeSet,veg%canst1,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'dleaf',completeSet,veg%dleaf,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'frac4',completeSet,veg%frac4,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'ejmax',completeSet,veg%ejmax,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vcmax',completeSet,veg%vcmax,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rp20',completeSet,veg%rp20,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'rpcoef',completeSet,veg%rpcoef,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'shelrb',completeSet,veg%shelrb,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'xfang',completeSet,veg%xfang,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'wai',completeSet,veg%wai,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vegcf',completeSet,veg%vegcf,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'extkn',completeSet,veg%extkn,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'tminvj',completeSet,veg%tminvj,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'tmaxvj',completeSet,veg%tmaxvj,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'vbeta',completeSet,veg%vbeta,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'xalbnir',completeSet,veg%xalbnir,filename%met, &
         nmetpatches,'def')
    CALL readpar(ncid_met,'meth',completeSet,veg%meth,filename%met, &
         nmetpatches,'def')
    ok = NF90_INQ_VARID(ncid_met,'za',parID)
    IF(ok == NF90_NOERR) THEN ! if it does exist
      CALL readpar(ncid_met,'za',completeSet,rough%za_uv,filename%met, &
         nmetpatches,'def')
      CALL readpar(ncid_met,'za',completeSet,rough%za_tq,filename%met, &
         nmetpatches,'def')
    ELSE
      CALL readpar(ncid_met,'za_uv',completeSet,rough%za_uv,filename%met, &
         nmetpatches,'def')
      CALL readpar(ncid_met,'za_tq',completeSet,rough%za_tq,filename%met, &
         nmetpatches,'def')
    ENDIF
    CALL readpar(ncid_met,'zse',completeSet,soil%zse,filename%met, &
         nmetpatches,'ms')
    CALL readpar(ncid_met,'ratecp',completeSet,bgc%ratecp,filename%met, &
         nmetpatches,'ncp')
    CALL readpar(ncid_met,'ratecs',completeSet,bgc%ratecs,filename%met, &
         nmetpatches,'ncs')

  END SUBROUTINE get_parameters_met
  !===========================================================================
  SUBROUTINE allocate_cable_vars(air,bgc,canopy,met,bal,rad,rough,soil,ssoil, &
       sum_flux,veg,arraysize)
    TYPE (met_type), INTENT(INOUT) :: met
    TYPE (air_type), INTENT(INOUT) :: air
    TYPE (soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE (bgc_pool_type), INTENT(INOUT)  :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (canopy_type), INTENT(INOUT)    :: canopy
    TYPE (roughness_type), INTENT(INOUT) :: rough
    TYPE (radiation_type),INTENT(INOUT)  :: rad
    TYPE (sum_flux_type), INTENT(INOUT)  :: sum_flux
    TYPE (balances_type), INTENT(INOUT)  :: bal
    INTEGER, INTENT(IN) :: arraysize

    ! Allocate CABLE's main variables:
    CALL alloc_cbm_var(air, arraysize)
    CALL alloc_cbm_var(bgc, arraysize)
    CALL alloc_cbm_var(canopy, arraysize)
    CALL alloc_cbm_var(met, arraysize)
    CALL alloc_cbm_var(bal, arraysize)
    CALL alloc_cbm_var(rad, arraysize)
    CALL alloc_cbm_var(rough, arraysize)
    CALL alloc_cbm_var(soil, arraysize)
    CALL alloc_cbm_var(ssoil, arraysize)
    CALL alloc_cbm_var(sum_flux, arraysize)
    CALL alloc_cbm_var(veg, arraysize)

    ! Allocate patch fraction variable:
    ALLOCATE(patch(arraysize))

  END SUBROUTINE allocate_cable_vars

  !==========================================================================
END MODULE input_module
