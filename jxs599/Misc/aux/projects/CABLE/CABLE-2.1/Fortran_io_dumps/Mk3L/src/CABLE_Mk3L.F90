! CABLE_Mk3L.f90
!
! Written by Jude Ambrose, Mao Jiafu, Gab Abramowitz, Steve Phipps, UNSW 2009
! Modified and added casaCNP by Bernard Pak, 2010
!
!#########################################################################
! CABLE_Mk3L_module declares all the necessary variables
!
! CABLE_init subroutine;
!  1. set CABLE's spatial and temporal variables (i.e., mland, dels and kend)
!     Assign latitude and longitude values from Mk3L
!     Allocate spatial heterogeneity variables(lanpt%cstart and cend)
!  2. Use usual load parameter and initialization routines 
!       (i)   get the default_global_vales
!       (ii)  Allocate the complete set of CABLE derived type variables to store
!     CABLE data for the WHOLE GLOBE (i.e., _global variables). 
!       - For all grid cells with only active patches in the landpoints
!       (iii) load the deault lai values
!       (iv) check for restart file; if exists load from restart file
!       (v)  construct derived parameters, check parameters, 
!             report parameters 
!  3. Determine how many land points are in each execution set of Mk3L, and
!     find out the indices of those land points within the Mk3L grid, as 
!     well as within the CABLE grid 
! 
! Routines within the timestep loop
!   Execution step loop   
! CABLE_run subroutine 
!  4. Set mland within CABLE to match execution step
!  5. Point all CABLE's variables to the subset of the global CABLE 
!     variable set (_global) that correspond to the land points within the
!     execution loop.
!  6. Write Mk3L met data to CABLE met variables for land points within the
!     execution loop  
!       - call get_default_lai for each landpoint in the step      
!  7. Call CABLE 
!  8. Write CABLE outputs to Mk3L variables for the execution loop 
!
! CABLE_write_output subroutine
! 
! CABLE_close subroutine
! 9. Write global restart file
! 10. Deallocate CABLE _global variables
!########################################################################
! SUBROUTINES WITH INPUTS AND OUTPUTS:
! 1. CABLE_init:
!       inputs (arguments):
!            lat    =  number of points in latitude             | = 28
!            lon    =  number of points in half longitude       | = 64
!            ln2    =  number of points in longitude            | = 128 Mk3L variables defined in PARAM.F
!            imsl   =  an  array of all the points in the globe | Mk3L variable declared in LSMI.f
!            sia    =  for latitude computation                 | declared in GAUSL.f
!            mstep  =  size of time step in minutes             | Mk3L variable declared in TIMEX.f
!            rrvco2 =  CO2 concentration                        |  Mk3L variable in RADISW.f
!            iyear  =  local time year AD
!            month  =  local time month of year
!            mins   =   minutes since 0000:01:01  00:00
!       outputs:  NONE
! 
! 2. CABLE_run
!       inputs (arguments):
!            lg      =  the number of the execution step currently executing
!            sg      =  net downward shortwave surface
!            rgsav   =  longwave (W/m2)
!            rcondx  =  rainfall(liquid+solid)(mm/dels)
!            pg      =  surface air pressure (mbar)
!            ttg     =  surface air temperature (K)
!            vmod    =  surface wind speed (m/s): computed in hsflux.f
!            zmin    =  height of the lowest level
!            iyear   =  local time year AD                               | 
!            month   =  local time month of year                         |
!            kdays   =  local time day of year = days since 0 hr 1st Jan | included in TIMEX.f
!            mins    =  minutes since 0000:01:01  00:00                  |
!            ln2     =  number of points for each step
!            nl      =  number of layers?
!            lat     =  number of points in latitude
!            imsl    =  an  array of all the points in the globe
!            nstepsa =  the execution step counter
!            qcloud  =  used for precip_s computation
!            qtg     =  water vapour mixing ratio (kg/kg)
!            rpreci  =  for precip_s computation
!            mstep   =  size of time step in minutes | Mk3L variable declared in TIMEX.f
!        outputs (arguments):
!            rg      =  net longwave heating at ground, lwnet from CABLE
!            scalev  =  scaling evaporation in mm/timestep.
!            tg      =  combined soil-canopy temperature of the surface layer
!            fg      =  sensible heat flux
!            eg      =  latent heat flux
!            tb2     =  temperature of the second layer
!            tb3     =  temperature of the lowest layer
!            totpev  =  potential evaporation in mm/timestep.
!            mcmax   =  maximum skin reservation depth
!            tddd    =  total wet evaporation
!            tgf     =  vegetation covered ground temperature
!            tggsl   =  vegetation covered ground temperature
!            als     =  albedo
!            snowd   =  snow depth
!            runoff  =  ground runoff, passed onto radin.f through surfa.f
!            wb      =  soil moisture 
!            wg      =  soil moisture available for evapotranspiration at the top layer
!            wg2     =  soil moisture available for evapotranspiration at the lower layer
!
! 3. CABLE_write_output
!        inputs (arguments):
!           nstepsa  =  the execution step counter
!        outputs: NONE
!
! 4. CABLE_close
!        inputs : NONE
!        outputs : NONE 
!#########################################################################

MODULE CABLE_Mk3L_module
  USE output_module  
  USE cbm_module  
  USE define_types
  USE define_dimensions ! mvtype, mstype now moved here sep2010
  USE io_variables  ! including namelist variables
  USE netcdf
  USE radiation_module,ONLY:sinbet
  USE physical_constants
  USE input_module
  USE parameter_module
  USE initialisation_module
  USE casa_cnp   ! mplant, mlitter, msoil, mphase and icycle
                 ! are declared in casadimension which is used in casa_cnp;
                 ! casafile is declared in casavariable, also used in casa_cnp;
                 ! also used modules casaparm and phenvariable

  IMPLICIT NONE
  ! CABLE variables for the entire globe:
  TYPE (met_type)       :: met_global
  TYPE (air_type)       :: air_global
  TYPE (soil_snow_type) :: ssoil_global
  TYPE (veg_parameter_type) :: veg_global
  TYPE (bgc_pool_type)  :: bgc_global
  TYPE (soil_parameter_type) :: soil_global
  TYPE (canopy_type)    :: canopy_global
  TYPE (roughness_type) :: rough_global
  TYPE (radiation_type) :: rad_global
  TYPE (sum_flux_type)  :: sum_flux_global
  TYPE (balances_type)  :: bal_global
  ! casaCNP variables for the entire globe: (BP apr2010)
  TYPE (casa_pool)      :: casapool_g
  TYPE (casa_flux)      :: casaflux_g
  TYPE (casa_met)       :: casamet_g
  TYPE (casa_balance)   :: casabal_g 
  TYPE (casa_biome)     :: casabiome ! N.B. casabiome for veg type only
  TYPE (phen_variable)  :: phen_global

  INTEGER(i_d) :: mland_global ! number of landpoints globally
  TYPE(land_type),DIMENSION(:),POINTER :: landpt_global ! # active patches, cstart, cend globally
  INTEGER :: mp_global ! Number of patches globally; length of global cable arrays
  TYPE(patch_type),DIMENSION(:),POINTER :: patch_global  ! patch fractions for global arrays
  REAL(r_1),POINTER,DIMENSION(:,:)  :: defaultLAI_global ! (BP apr2010)
  LOGICAL :: vegparmnew ! using new format input file (BP dec 2007)
  ! For each Mk3L execution (latitude) set, idx_start/end gives index of land points
  INTEGER(i_d), ALLOCATABLE, DIMENSION(:) :: idx_start ! in CABLE's global arrays
  INTEGER(i_d), ALLOCATABLE, DIMENSION(:) :: idx_end   

  ! Variables needed from Mk3L are declared
  REAL(r_1) :: dels          ! size of the time step in seconds
  INTEGER(i_d) :: kend      ! number of timesteps in the run
  REAL(r_2) :: rrvco2_global ! co2 concentration passed from Mk3L
  INTEGER(i_d) :: kstart = 1 ! start of simulation #
  INTEGER (i_d) :: ktau        ! total run timestep
  INTEGER :: ok
  INTEGER :: threadID
  !$OMP THREADPRIVATE(ktau)
  !$OMP THREADPRIVATE(threadID)

END MODULE CABLE_Mk3L_module


!###########################################################
!####             CABLE_init subroutine            #########
!###########################################################

SUBROUTINE CABLE_init(imsl, lon, ln2, lat, sia, mstep, rrvco2, & ! INPUTS
     & iyear,month,mins)   ! INPUTS
  USE CABLE_Mk3L_module
  USE math_constants                               ! to access pi_c 
  IMPLICIT NONE
  ! #############    INPUTS FROM Mk3L   ########################
  INTEGER(i_d), INTENT(IN) :: lat ! number of points in latitude | = 28
  INTEGER(i_d), INTENT(IN) :: lon ! number of points in half longitude | = 64
  INTEGER(i_d), INTENT(IN) :: ln2 ! number of points in longitude | = 128 Mk3L variables defined in PARAM.F
  ! an  array of all the points in the globe | Mk3L variable declared in LSMI.f:
  INTEGER(i_d), INTENT(IN), DIMENSION(ln2, lat) :: imsl 
  REAL(r_2), INTENT(IN), DIMENSION(lat) :: sia ! for latitude computation | declared in GAUSL.f
  INTEGER(i_d), INTENT(IN) :: mstep ! size of time step in minutes | Mk3L variable declared in TIMEX.f
  REAL(r_2), INTENT(IN) :: rrvco2  ! CO2 concentration |  Mk3L variable in RADISW.f      
  INTEGER(i_d),INTENT(IN) :: iyear ! local time year AD
  INTEGER(i_d), INTENT(IN) :: month ! local time month of year
  INTEGER(i_d), INTENT(IN) :: mins  ! minutes since 0000:01:01  00:00  
  ! ########################################################### 

  ! ################### LOCAL VARIABLES #######################
  INTEGER(i_d) :: i,j                              ! temporary variables for loops
  INTEGER(i_d) :: tmp                              ! temporary variable
  CHARACTER(LEN=13) :: timestring                  ! will have a fixed value "seconds since"   
  ! ###########################################################
  ! CABLE namelist values are initialized
  NAMELIST/CABLE/filename,vegparmnew,check,verbose,leaps,         & 
       logn,fixedCO2,output,icycle,casafile
  !===================================================================!
  ! Open, read and close the CABLE namelist file.
  OPEN(10,FILE='cable_mk3l.nml')
  READ(10,NML=CABLE)
  CLOSE(10)
  !=====================================================================!

  ! CABLE grid dimensions based on the Mk3L 56x64:
  xdimsize = 64
  ydimsize = 56
  ! Set CABLE grid to mask - mask array will be used
  metGrid = 'mask' 
  ! rrvco2 value mapped to the global variable
  rrvco2_global = rrvco2
  ! define time variables
!!$      ALLOCATE(timevar(999999)) ! just allocated to get it working to the value of kend.. REMOVE AT THE END
  time_coord = 'GMT' ! set to GMT and the other option is LOC  
  ! timeunits should have this string "seconds since 2002-01-01 00:01:00"
  ! replacing the current date and time values 
  timestring = '                               '
  WRITE(timeunits,*) timestring
  timestring = 'seconds since'          
  WRITE(timeunits(1:13),'(A13)')  timestring
  WRITE(timeunits(15:18),'(I4)') iyear
  IF(iyear < 10) THEN  ! if the year is one digit concatenate with 0
     timestring = '000' // timeunits(18:18)
     WRITE(timeunits(15:18),'(A4)') timestring
  ELSE IF(iyear > 10 .AND. iyear < 100) THEN
     timestring = '00' // timeunits(17:18)
     WRITE(timeunits(15:18),'(A4)') timestring
  ELSE IF(iyear > 100 .AND. iyear < 1000) THEN
     timestring = '0' // timeunits(16:18)
     WRITE(timeunits(15:18),'(A4)') timestring
  END IF
  timestring = '-'
  WRITE(timeunits(19:19),'(A1)')  timestring
  WRITE(timeunits(20:21),'(I2)') month
  IF(month < 10) THEN  ! if the month is one digit concatenate with 0
     timestring = '0' // timeunits(21:21)
     WRITE(timeunits(20:21),'(A2)') timestring
  END IF
  timestring = '-'
  WRITE(timeunits(22:22),'(A1)')  timestring
  timestring = '01'
  WRITE(timeunits(23:24),'(A2)') timestring
  tmp = MOD(mins, 1440)/60  
  WRITE(timeunits(26:27),'(I2)') tmp
  IF(month < 10) THEN  ! if the month is one digit concatenate with 0
     timestring = '0' // timeunits(27:27)
     WRITE(timeunits(26:27),'(A2)') timestring
  END IF
  timestring = ':00:00'
  WRITE(timeunits(28:33),'(A6)')  timestring

  ! Open log file:
  OPEN(logn,FILE=filename%log)

  ! Allocate CABLE longitude and latitude arrays for entire grid:
  ALLOCATE(lon_all(xdimsize,ydimsize))
  ALLOCATE(lat_all(xdimsize,ydimsize))
  
  ! Allocate land/sea mask variable:
  ALLOCATE(mask(xdimsize,ydimsize))

  ! 1. Values needed from Mk3L are set
  kend = INT(365 * 24 * 60/mstep) ! set number of timesteps in a year run
  dels = REAL(mstep,r_1)*60.0  ! convert minutes to seconds
  ! Find the number of land points in the Mk3L global grid (mland_global):
  mland_global = 0
  DO j = 1, lat
     DO i = 1, ln2
        IF(imsl(i,j) == 4) THEN
           mland_global = mland_global + 1    ! increment the landpoint counter
        END IF
     END DO
  END DO
  ! Temporarily assign # land points (mland) to global # land points
  ! (mland will later represent # land points in a latitude band in Mk3L):
  mland = mland_global 

  ! Allocate latitude and longitude variables for land points:
  ALLOCATE(latitude(mland_global),longitude(mland_global))
  ! Allocate variables for x and y indicies of land points in global x-y grid
  ALLOCATE(land_x(mland_global),land_y(mland_global))
  mask = 0    ! Initialise all gridpoints as sea:
  tmp = 1     ! land point counter
  ! Find latitude/longitude, mask and x, y indicies of all land points:
  DO j = 1, lat 
     DO i = 1, ln2
        IF(imsl(i,j)== 4) THEN ! If land point
           IF(i <= lon) THEN        ! compute values when less or equal to lon
              latitude(tmp)  = REAL(ACOS(sia(j))*180.0/pi_c,r_1)
              longitude(tmp) = REAL(i-1)*360.0/lon
              mask(i,2*lat-j+1) = 1  ! setting up the mask variable to 1 only for landpts
              land_x(tmp)    = i     ! setting up the x and y indices of the landpt
              land_y(tmp)    = 2*lat-j+1          
           ELSE     ! compute values when greater than lon
              latitude(tmp)  = REAL(0.0-ACOS(sia(j))*180.0/pi_c,r_1)
              longitude(tmp) = REAL(i-lon-1)*360.0/lon
              mask(i-lon,j)  = 1     ! setting up the mask variable to 1 only for landpts
              land_x(tmp)    = i-lon ! setting up the x and y indices of the landpt
              land_y(tmp)    = j     
           END IF
           tmp = tmp + 1  ! increment counter for the next land point  
        END IF
        ! Setting values for all grid point latitudes and longitudes, where
        ! the grid of 28x128 from Mk3L is mapped to 56x64 for CABLE
        IF(i <= lon) THEN    
           lat_all(i,2*lat-j+1)  = REAL(ACOS(sia(j))*180.0/pi_c,r_1)
           lon_all(i,2*lat-j+1)  = REAL(i-1)*360.0/lon
        ELSE
           lat_all(i-lon,j)  = REAL(0.0-ACOS(sia(j))*180.0/pi_c,r_1)
           lon_all(i-lon,j)  = REAL(i-lon-1)*360.0/lon
        END IF
     END DO
  END DO

  ! Longitude values varied from 0 -360 from Mk3L are converted between
  ! -180 and 180 for CABLE   
  WHERE (longitude > 180.0)
     longitude = longitude - 360.0
  END WHERE

  ! Allocate spatial heterogeneity variables:
  ALLOCATE(landpt(mland_global))
  ALLOCATE(landpt_global(mland_global))

  ! Computations of landpt%cstart, cend and nap are performed in
  ! cable_parameters/get_default_params. That is, landpt
  ! is first set as a global array below.

  ! 2. Loading parameters and initialization routines
  ! default_params has to be called to load default grid 
  ! The parameters and initialisations it writes
  ! will be overwritten by the restart 

  ! Initialisations:
  nmetpatches = 0 ! i.e. there is no met forcing file for CABLE

  ! 2.(i)
  WRITE(logn,*) ' Loading initialisations ', &
       'from default grid. '
  CALL get_default_params(logn,vegparmnew,lon,ln2,lat)
  ! assign the size of mp to mp_global
  mp_global = mp   

  ! 2.(ii) spaces allocated for all globe arrays, including 'patch' from iovar
  CALL allocate_cable_vars(air_global,bgc_global,canopy_global,met_global, & 
        bal_global,rad_global,rough_global,soil_global,ssoil_global, &
        sum_flux_global,veg_global,mp)

  ! 2.(iii) 
  CALL write_default_params(met_global,air_global,ssoil_global,veg_global, &
        bgc_global,soil_global,canopy_global,rough_global,rad_global,logn, &
        vegparmnew,month) 
!  WRITE(45,*) 'tile, flag, sdepth(1-3),                tggsn(1-3)'
!  WRITE(45,*) 'After write_default_params: '
!  DO i = 1, mp_global
!    IF (ssoil_global%isflag(i) > 0 .AND. ANY(ssoil_global%sdepth(i,:) == 0.0)) &
!      WRITE(45,*) i, ssoil_global%isflag(i), &
!                     ssoil_global%sdepth(i,:), ssoil_global%tggsn(i,:)
!  END DO
!  WRITE(45,'(6x,6f7.2)') ssoil_global%tgg(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wb(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wbice(1672,:)

  ! initializing the appropriate snow parameters based on Mk3L values 
  smoy = 1
  CALL get_default_inits(met_global,soil_global,ssoil_global,canopy_global,logn)
!  WRITE(45,*) 'After get_default_inits: '
!  DO i = 1, mp_global
!    IF (ssoil_global%isflag(i) > 0 .AND. ANY(ssoil_global%sdepth(i,:) == 0.0)) &
!      WRITE(45,'(2i5,6f10.4)') i, ssoil_global%isflag(i), &
!                     ssoil_global%sdepth(i,:), ssoil_global%tggsn(i,:)
!  END DO
!  WRITE(45,'(6x,6f7.2)') ssoil_global%tgg(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wb(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wbice(1672,:)

  ! mapping patch to patch global array
  ALLOCATE(patch_global(mp_global)) 
  ALLOCATE(defaultLAI_global(mp_global,12))
  patch_global = patch
  ! Setting rad%latitude values to be used in the CABLE_run 
  DO i=1,mland_global 
     rad_global%latitude(landpt(i)%cstart:landpt(i)%cend) = latitude(i)
  END DO

  ! 2.(iv): load the default LAI only to open the LAI file
!  ktau = 1 ! set to one to execute the first portion of the code
          ! in get_default_lai ! commented out by BP because not used (May2010)
  met_global%moy = 1 ! set to a dummy value one 
  CALL get_default_lai 
  defaultLAI_global = defaultLAI

  !  Look for restart file (which will have all parameters and inits):
  ok = NF90_OPEN(filename%restart_in,0,ncid_rin) ! open restart file
  IF (ok /= NF90_NOERR) THEN
     WRITE(logn,*) ' Could not find restart file: ',&
          TRIM(filename%restart_in)
  ELSE ! RESTART FILE EXISTS, parameters and init will be loaded from it.
     WRITE(logn,*) ' Loading initialisations ', &
          'from restart file: ', TRIM(filename%restart_in)
     ! Load initialisations and parameters from restart file:
     CALL get_restart_data(logn,ssoil_global,canopy_global,rough_global,bgc_global, &
          & bal_global,veg_global,soil_global,rad_global,vegparmnew)     
  END IF ! if restart file exists
!  WRITE(45,*) 'After get_restart_data: '
!  DO i = 1, mp_global
!    IF (ssoil_global%isflag(i) > 0 .AND. ANY(ssoil_global%sdepth(i,:) == 0.0)) &
!      WRITE(45,'(2i5,6f10.4)') i, ssoil_global%isflag(i), &
!                     ssoil_global%sdepth(i,:), ssoil_global%tggsn(i,:)
!  END DO
!  WRITE(45,'(6x,6f7.2)') ssoil_global%tgg(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wb(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wbice(1672,:)

  ! checking for wb and wbice whether they are greater than ssat and resetting
  ! back to appropriate values
  ! Also check for very small wbice values when frozen ground thaws; to prevent
  ! floating point underflow problem (BP jul2010)
  DO tmp=1,ms
     WHERE(ssoil_global%wb(:,tmp)>soil_global%ssat)
        ssoil_global%wb(:,tmp) = soil_global%ssat
     END WHERE
     WHERE(ssoil_global%wbice(:,tmp)>0.98*soil_global%ssat)
       ssoil_global%wbice(:,tmp) = 0.98*soil_global%ssat
     END WHERE
     WHERE(ssoil_global%wbice(:,tmp)<1.0e-6)
       ssoil_global%wbice(:,tmp) = 0.0
     END WHERE
  END DO
  ! ################################
!  WRITE(45,*) 'After reset: '
!  WRITE(45,'(6x,6f7.2)') ssoil_global%tgg(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wb(1672,:)
!  WRITE(45,'(6x,6f7.3)') ssoil_global%wbice(1672,:)


  ! 2.(v): Construct derived parameters, regardless
  ! of where parameters have loaded from:
  CALL derived_parameters(soil_global,sum_flux_global,bal_global,ssoil_global, &
                           veg_global,rough_global)


  ! initialization of casaCNP (BP apr2010)
  IF (icycle >0) THEN
    PRINT *, 'mp before alloc_casavariable = ', mp
    PRINT *, '    and mland_global = ', mland_global
    CALL alloc_casavariable(casabiome,casapool_g,casaflux_g,casamet_g, &
                            casabal_g,mp)
    CALL alloc_phenvariable(phen_global,mp)
    call casa_readpoint(veg_global,soil_global,casaflux_g,casamet_g,rad_global)
    call casa_readbiome(veg_global,soil_global,casabiome, &
                        casapool_g,casaflux_g,casamet_g,phen_global)
    call casa_readphen(veg_global,casamet_g,phen_global)
    call casa_init(casapool_g,casabal_g,veg_global)
! added mvtype and mstype to define_dimensions (BP sep2010)
!    CALL alloc_casavariable(casabiome,casapool_g,casaflux_g,casamet_g, &
!                            casabal_g,mp,mvtype)
!    CALL alloc_phenvariable(phen_global,mp,mvtype)
!    call casa_readpoint(mvtype,veg_global,soil_global,casaflux_g,casamet_g, &
!                        patch_global,rad_global)
!    call casa_readbiome(mvtype,mstype,veg_global,soil_global,casabiome, &
!                        casapool_g,casaflux_g,casamet_g,phen_global)
!    call casa_readphen(mvtype,veg_global,casamet_g,phen_global)
!! for first time reading file *_1220.csv  (BP may2010)
!    call casa_init(mstype,casapool_g,casabal_g,veg_global)
!!    call casa_init(mstype,casapool_g,casabal_g)
!! end addition (BP may2010)
  END IF


  ! Check for basic inconsistencies in parameter values:
  CALL check_parameter_values(soil_global,veg_global,ssoil_global)

  ! Write per-site parameter values to log file if requested:
  IF(verbose) CALL report_parameters(logn,soil_global,veg_global,bgc_global, &
       rough_global,ssoil_global,canopy_global,vegparmnew,verbose)

  ! 3. the imsl array from Mk3L is a latxln2 (i.e., 28x128) array where each
  ! row is a vector used for an execution. For each execution set (each lat) 
  ! record the first land point and last land point in idx_start
  ! and idx_end respectively.
  ALLOCATE (idx_start(lat),idx_end(lat)) 
!  ALLOCATE (idx_start(lat), STAT = ok) 
!  ALLOCATE (idx_end(lat), STAT = ok)
  tmp = 0 ! init  
  DO j = 1,lat   ! for every execution set      
     DO i = 1, ln2 ! for each grid in an execution set
        IF(imsl(i,j) == 4) THEN  ! increase tmp when there is land found
           tmp =  tmp + 1           
        END IF
     END DO
     IF (tmp > 0) THEN        ! when tmp greater than one there is a land found
        IF(j == 1) THEN    ! check whether its the first execution set 
           idx_start(1) = 1  ! the start index is set to one 
        ELSE
           ! the start index for the subsequent rows would be the end index
           ! of the previous execution set + 1
           idx_start(j) = idx_end(j - 1) + 1  
        END IF
        idx_end(j) = tmp   ! end index is the final tmp from each execution set
     ELSE
        idx_start(j) = 0   ! set start index to 0 when tmp is zero
        idx_end(j) = 0     ! set end index to 0 when tmp is zero
     END IF
  END DO
!  WRITE(42,*) 'idx_start:'
!  WRITE(42,'(99(i5))') idx_start
!  WRITE(42,*) 'idx_end:'
!  WRITE(42,'(99(i5))') idx_end

  ! Global cstart, cend and nap (in landpt) are those collected above where mland=mland_global
  landpt_global = landpt  

!    tmp = 0
!    DO i = 1, mland_global
!    DO j = landpt_global(i)%cstart, landpt_global(i)%cend
!      tmp = tmp + 1
!      WRITE(77,'(i5,2f10.3,2i5,f10.3,2i5)') tmp, patch_global(tmp)%latitude, &
!            patch_global(tmp)%longitude, veg_global%iveg(tmp), &
!            soil_global%isoilm(tmp), patch_global(tmp)%frac, &
!            landpt_global(i)%ilon, landpt_global(i)%ilat
!    ENDDO
!    ENDDO

  ! opening the output file
  CALL open_output_file(dels,soil_global,veg_global,bgc_global,rough_global)

  DEALLOCATE(landpt) ! deallocate landpt so it can be locally reallocated for each execution step
  DEALLOCATE(patch) ! deallocate patch - it will be pointed in each execution set

END SUBROUTINE CABLE_init  ! end of initializations 


!###########################################################
!####             CABLE_run subroutine            ##########
!###########################################################
   
SUBROUTINE CABLE_run(lg,sg,rgsav,rcondx,pg,ttg,iyear,month,kdays,mins, & !INPUTS
     ln2,lat,nl,imsl,nstepsa,vmod,zmin,qcloud,qtg,rpreci,mstep,  &  ! INPUTS
     rg,scalev,tg,fg,eg,tb2,tb3,totpev,mcmax,tddd,tgg,tgf,mc,osnowd, & ! OUTPUTS
     snage,ssdnn,gflux,sgflux,tggsl,tggsn,wb,wbice,smass,ssdn3,  & ! OUTPUTS
     isflag,als,snowd,runoff,perc,wg,wg2,cabTscrn)                 ! OUTPUTS
  USE CABLE_Mk3L_module
  IMPLICIT NONE
  !######################################
  !#####      input variables      ######
  !######################################    
  INTEGER(i_d), INTENT(IN) ::  lg       ! the number of the execution step currently executing
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: sg     ! net downward shortwave surface
  REAL(r_2),DIMENSION(ln2,lat), INTENT(IN) :: rgsav   ! longwave (W/m2)
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: rcondx  ! rainfall (liquid+solid)(mm/dels)
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: pg      ! surface air pressure (mbar)
  REAL(r_2),DIMENSION(ln2,nl), INTENT(IN) :: ttg  ! surface air temperature (K)
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: vmod    ! surface wind speed (m/s): computed in hsflux.f      
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: zmin    ! height of the lowest level     
  INTEGER(i_d),INTENT(IN) :: iyear                ! local time year AD                               |
  INTEGER(i_d), INTENT(IN) :: month ! local time month of year                         |
  INTEGER(i_d), INTENT(IN) :: kdays ! local time day of year = days since 0 hr 1st Jan | included in TIMEX.f
  INTEGER(i_d), INTENT(IN) :: mins  ! minutes since 0000:01:01  00:00                  |   
  INTEGER(i_d), INTENT(IN) :: ln2       ! number of points for each step  
  INTEGER(i_d), INTENT(IN) :: nl        ! number of layers?
  INTEGER(i_d), INTENT(IN) :: lat       ! number of points in latitude 
  INTEGER(i_d), INTENT(IN), DIMENSION(ln2, lat) :: imsl ! an  array of all the points in the globe
  INTEGER(i_d), INTENT(IN) :: nstepsa               ! the execution step counter
  LOGICAL :: qcloud                                 ! used for precip_s computation
  REAL(r_2),DIMENSION(ln2,nl), INTENT(IN) :: qtg      ! water vapour mixing ratio (kg/kg)
  REAL(r_2),DIMENSION(ln2), INTENT(IN) :: rpreci      ! for precip_s computation
  INTEGER(i_d), INTENT(IN) :: mstep ! size of time step in minutes | Mk3L variable declared in TIMEX.f
  !########################################
  !#####       output variables    ########
  !########################################
  ! compared to surfa.f of Mk3L
  ! some variables are made INOUT only to make sure the non landpoints
  ! preserve their appropriate values
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: cabTscrn ! screen temperature from CABLE
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: rg      ! net longwave heating at ground, lwnet from CABLE
  REAL(r_2), DIMENSION(ln2), INTENT(OUT) :: scalev  ! scaling evaporation in mm/timestep.
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: tg      ! combined soil-canopy temperature of the surface layer
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: fg      ! sensible heat flux
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: eg      ! latent heat flux
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: tb2     ! temperature of the second layer
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: tb3     ! temperature of the lowest layer
  REAL(r_2), DIMENSION(ln2), INTENT(OUT) :: totpev  ! potential evaporation in mm/timestep.
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: mcmax   ! maximum skin reservation depth
  REAL(r_2), DIMENSION(ln2), INTENT(OUT) :: tddd    ! total wet evaporation
  ! The whole SURF1.f common block
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: tgg  ! ground temperature
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: tgf     ! vegetation covered ground temperature
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: mc
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: osnowd  ! snow depth of last time step
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: snage  ! snow age
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: ssdnn  ! average snow density
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: gflux  ! ground heat flux
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: sgflux ! snow-covered ground heat flux
  REAL(r_2), DIMENSION(ln2,ms), INTENT(INOUT) :: wb   ! soil moisture ! size of ms = 6
  REAL(r_2), DIMENSION(ln2,ms), INTENT(INOUT) :: wbice
  REAL(r_2), DIMENSION(ln2,ms), INTENT(INOUT) :: tggsl ! vegetation covered ground temperature
  REAL(r_2), DIMENSION(ln2,3), INTENT(INOUT) :: tggsn ! snow layer temperature
  REAL(r_2), DIMENSION(ln2,3), INTENT(INOUT) :: smass ! snow mass
  REAL(r_2), DIMENSION(ln2,3), INTENT(INOUT) :: ssdn3 ! snow-layer density
  INTEGER(i_d), DIMENSION(ln2), INTENT(INOUT) :: isflag
  ! end SURF1.f common block
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: als  ! albedo 

  ! compared to surfb.f of Mk3L
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: snowd   ! snow depth
  REAL(r_2), DIMENSION(ln2), INTENT(OUT) :: runoff  ! ground runoff, passed onto radin.f through surfa.f      
  REAL(r_2), DIMENSION(ln2), INTENT(OUT) :: perc ! soil moisture percolation
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: wg      ! soil moisture available for evapotranspiration at the top layer
  REAL(r_2), DIMENSION(ln2), INTENT(INOUT) :: wg2      ! soil moisture available for evapotranspiration at the lower layer

  !########################################
  !######     local variables      ########
  !########################################
  INTEGER(i_d) :: i,j,tmp               ! temporary variables for loops      
  INTEGER(i_d) :: idx_start_patch         ! temporary var to point to the start
  INTEGER(i_d) :: idx_end_patch           ! temporary var to point to the end
  REAL(r_1), DIMENSION(ms) :: zse2       ! temporary var to hold zse

  ! variables declared for each execution step
  TYPE (air_type)       :: air  ! air property variables
  TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
  TYPE (canopy_type)    :: canopy ! vegetation variables
  TYPE (met_type)       :: met  ! met input variables
  TYPE (balances_type)  :: bal  ! energy and water balance variables
  TYPE (radiation_type) :: rad  ! radiation variables
  TYPE (roughness_type) :: rough ! roughness varibles
  TYPE (soil_parameter_type) :: soil ! soil parameters
  TYPE (soil_snow_type) :: ssoil ! soil and snow variables
  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters

  ! casaCNP variables for a latitude band: (BP apr2010)
  TYPE (casa_pool)      :: casapool
  TYPE (casa_flux)      :: casaflux
  TYPE (casa_met)       :: casamet
  TYPE (casa_balance)   :: casabal
  TYPE (phen_variable)  :: phen

  INTEGER :: omp_get_thread_num
  INTEGER ::  ii, jj

  ! assigning the execution step counter to ktau, to denote which execution
  ! step is currently executing
!  ktau = (iyear-1)*365*24*60/mstep + nstepsa + 1 
!  ktau = iyear*365*24*60/mstep + nstepsa + 1 
  ktau = nstepsa + 1 

  threadID = omp_get_thread_num()

!bp Output diagnostics at several points to check on CABLE
  IF (ktau==1 .AND. lg==1) THEN
!    write(36,'(a64,a64,a70,a68,a70,a66,a47)')  &
!      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',& 
!      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
!      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
!      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
!      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
!      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
!      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
!    write(37,'(a64,a64,a70,a68,a70,a66,a47)')  &
!      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
!      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
!      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
!      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
!      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
!      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
!      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
!    write(38,'(a64,a64,a70,a68,a70,a66,a47)')  &
!      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
!      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
!      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
!      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
!      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
!      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
!      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
!    write(39,'(a64,a64,a70,a68,a70,a66,a47)')  &
!      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
!      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
!      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
!      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
!      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
!      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
!      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
    write(51,'(a64,a64,a70,a68,a70,a66,a47)')  &
      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
    write(52,'(a64,a64,a70,a68,a70,a66,a47)')  &
      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
    write(53,'(a64,a64,a70,a68,a70,a66,a47)')  &
      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
!    write(54,'(a64,a64,a70,a68,a70,a66,a47)')  &
!      'lat, lon, iveg, ktau, tk, tvair, tv, tss, fvlai1, fvlai2, vlaiw,',&
!      ' fsd, fld, extkb, extkd, qcan(:,1,1), qcan(:,1,2), qcan(:,1,3), ', &
!      'qcan(:,2,1), qcan(:,2,2), qcan(:,2,3), qssabs, albsoilsn1, albsoilsn2,',&
!      ' albsoilsn3, albedo1, albedo2, lwabv, transd, fevc, fwet, fhv, fes, ', &
!      'fhs, fns, fnv, trad, ebal, ebal_total, flws, fev, cls, fh, ga, snowd, ',&
!      'tgg1-6, tggsn1-3, wb1-6, wbice1-6, osnowd, snage, ssdnn, ssdn1-3, ', &
!      'smass1-3, ghflux, sghflux, delwc, rnof1, runof2'
  ENDIF
!  IF (lg .eq. 23) THEN 
!    WRITE(45,'(3f10.2,8i10,2f10.2,2x,a8,f10.2,i10,24f10.2,i10,6f10.2)') &
!       sg(4),rcondx(4),pg(4),iyear,month,kdays,mins, &
!       ln2,lat,nl,nstepsa,vmod(4),zmin(4),qcloud,rpreci(4),mstep, &
!       rg(4),scalev(4),tg(4),fg(4),eg(4),tb2(4),tb3(4),totpev(4), &
!       mcmax(4),tddd(4),tgg(4),tgf(4),mc(4),osnowd(4), &
!       snage(4),ssdnn(4),gflux(4),sgflux(4),tggsl(4,1), &
!       tggsn(4,1),wb(4,1),wbice(4,1),smass(4,1),ssdn3(4,1), &
!       isflag(4),als(4),snowd(4),runoff(4),wg(4),wg2(4),cabTscrn(4)
!  ENDIF

  ! 4. Set the number of land points, mland, for each execution step:
  mland = idx_end(lg) - idx_start(lg) + 1
  
  !##########################################################################
  ! below cstart and cend are reset so that they will be appopriate for
  ! referencing cable's LOCAL variables, rather than global variables.
  ! note that any references within the cbm or soilsnow code that use
  ! cstart or cend are now going to be changed
  ! #########################################################################
  
  ! Compute the cstart and cend for each landpoint in the current execution set
  ALLOCATE(landpt(mland))  ! Allocate to local mland size
  
  ! Get number of active patches from global array:
  DO i= 1,mland
     landpt(i)%nap = landpt_global(idx_start(lg)+i-1)%nap 
  END DO
  ! Set cstart and cend relative to local variables:
  DO i= 1,mland
     IF(i == 1) THEN
        landpt(i)%cstart = 1 
        landpt(i)%cend = landpt(i)%cstart +  landpt(i)%nap - 1
     ELSE
        landpt(i)%cstart =  landpt(i-1)%cend + 1  
        landpt(i)%cend =   landpt(i)%cstart +  landpt(i)%nap - 1
     END IF
  END DO
  ! Set the total # patches in this execution set to be the end patch index
  ! of the last land point in this execution set:
  mp = landpt(mland)%cend

  ! Set the start and end indicies in the global arrays for all patches in this execution set
  idx_start_patch = landpt_global(idx_start(lg))%cstart   ! the starting index of step  lg
  idx_end_patch   = landpt_global(idx_end(lg))%cend       ! the end index of step lg

  ! Point local air variables to the appropriate subset of the global array:
  air%rho    => air_global%rho(idx_start_patch : idx_end_patch)   ! dry airdensity (kg m-3)
  air%volm   => air_global%volm(idx_start_patch : idx_end_patch)  ! molar volume(m3 mol-1)
  air%rlam   => air_global%rlam(idx_start_patch : idx_end_patch)  ! latent heatfor water (j/kg)
  air%qsat   => air_global%qsat(idx_start_patch : idx_end_patch)  ! saturation specific humidity
  air%epsi   => air_global%epsi(idx_start_patch : idx_end_patch)  ! d(qsat)/dT ((kg/kg)/K)
  air%visc   => air_global%visc(idx_start_patch : idx_end_patch)  ! airkinematic viscosity (m2/s)
  air%psyc   => air_global%psyc(idx_start_patch : idx_end_patch)  ! psychrometric constant
  air%dsatdk => air_global%dsatdk(idx_start_patch : idx_end_patch)! d(es)/dT (mb/K) 
  air%cmolar => air_global%cmolar(idx_start_patch : idx_end_patch)! conv. from m/s to mol/m2/s

!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)

  ! Point local carbon pool variables to the appropriate subset of the global array:
  bgc%cplant => bgc_global%cplant(idx_start_patch : idx_end_patch,:)! plant carbon (g C/m2))
  bgc%csoil  => bgc_global%csoil(idx_start_patch : idx_end_patch,:) ! soil carbon (g C/m2)
  bgc%ratecp = bgc_global%ratecp
  bgc%ratecs = bgc_global%ratecs   

  ! Point local canopy variables to the appropriate subset of the global array:
  canopy%cansto => canopy_global%cansto(idx_start_patch : idx_end_patch) ! canopy water storage (mm)      
  canopy%cduv   => canopy_global%cduv(idx_start_patch : idx_end_patch)   ! drag coefficient for momentum
  canopy%delwc  => canopy_global%delwc(idx_start_patch : idx_end_patch)  ! change in canopy water store (mm/dels)
  canopy%dewmm  => canopy_global%dewmm(idx_start_patch : idx_end_patch)  ! dewfall (mm)
  canopy%dgdtg  => canopy_global%dgdtg(idx_start_patch : idx_end_patch)  ! derivative of gflux wrt soil temp
  canopy%fe     => canopy_global%fe(idx_start_patch : idx_end_patch)     ! total latent heat (W/m2)
  canopy%fh     => canopy_global%fh(idx_start_patch : idx_end_patch)     ! total sensible heat (W/m2)
  canopy%fpn    => canopy_global%fpn(idx_start_patch : idx_end_patch)    ! plant photosynthesis (g C s-1)
  canopy%frp    => canopy_global%frp(idx_start_patch : idx_end_patch)    ! plant respiration (g C m-2 s-1) 
  canopy%frpw   => canopy_global%frpw(idx_start_patch : idx_end_patch)   ! plant respiration (g C m-2 s-1)
  canopy%frpr   => canopy_global%frpr(idx_start_patch : idx_end_patch)   ! plant respiration (g C m-2 s-1)
  canopy%frs    => canopy_global%frs(idx_start_patch : idx_end_patch)    ! soil respiration (g C m-2 s-1)
  canopy%fnee   => canopy_global%fnee(idx_start_patch : idx_end_patch)   ! net carbon flux (g C m-2 s-1)
  canopy%frday  => canopy_global%frday(idx_start_patch : idx_end_patch)  ! daytime leaf resp
  canopy%fnv    => canopy_global%fnv(idx_start_patch : idx_end_patch)    ! net rad. avail. to canopy (W/m2)
  canopy%fev    => canopy_global%fev(idx_start_patch : idx_end_patch)    ! latent hf from canopy (W/m2)
  canopy%fevc   => canopy_global%fevc(idx_start_patch : idx_end_patch)   ! dry canopy transpiration (W/m2)
  canopy%fevw   => canopy_global%fevw(idx_start_patch : idx_end_patch)   ! lat heat fl wet canopy (W/m2)
  canopy%potev_c=> canopy_global%potev_c(idx_start_patch : idx_end_patch)! canopy potential evapotranspitation (YP & Mao, jun08)
  canopy%fhv    => canopy_global%fhv(idx_start_patch : idx_end_patch)    ! sens heatfl from canopy (W/m2)
  canopy%fhvw   => canopy_global%fhvw(idx_start_patch : idx_end_patch)   ! sens heatfl from wet canopy (W/m2)
  canopy%fns    => canopy_global%fns(idx_start_patch : idx_end_patch)    ! net rad avail to soil (W/m2)
  canopy%fes    => canopy_global%fes(idx_start_patch : idx_end_patch)    ! latent heatfl from soil (W/m2)
  canopy%fhs    => canopy_global%fhs(idx_start_patch : idx_end_patch)    ! sensible heat flux from soil
  canopy%fwet   => canopy_global%fwet(idx_start_patch : idx_end_patch)   ! fraction of canopy wet
  canopy%ga     => canopy_global%ga(idx_start_patch : idx_end_patch)     ! ground heat flux (W/m2)
  canopy%ghflux => canopy_global%ghflux(idx_start_patch : idx_end_patch) ! ground heat flux (W/m2)
  canopy%precis => canopy_global%precis(idx_start_patch : idx_end_patch) ! throughfall to soil, after snow (mm)
  canopy%qscrn  => canopy_global%qscrn(idx_start_patch : idx_end_patch)  ! specific humudity at screen height (g/g)
  canopy%rnet   => canopy_global%rnet(idx_start_patch : idx_end_patch)   ! net radiation absorbed by surface (W/m2)
  canopy%segg   => canopy_global%segg(idx_start_patch : idx_end_patch)   ! latent heatfl from soil mm (EK nov 2007)
  canopy%sghflux=> canopy_global%sghflux(idx_start_patch : idx_end_patch)! ground heat flux (W/m2)
  canopy%spill  => canopy_global%spill(idx_start_patch : idx_end_patch)  ! ground heat flux (W/m2)
  canopy%through=> canopy_global%through(idx_start_patch : idx_end_patch)! canopy throughfall (mm)
  canopy%tscrn  => canopy_global%tscrn(idx_start_patch : idx_end_patch)  ! air temperature at screen height (oC)
  canopy%tv     => canopy_global%tv(idx_start_patch : idx_end_patch)     ! vegetation temp (K)
  canopy%us     => canopy_global%us(idx_start_patch : idx_end_patch)     ! friction velocity
  canopy%uscrn  => canopy_global%uscrn(idx_start_patch : idx_end_patch)  ! wind speed at screen height (m/s)
  canopy%vlaiw  => canopy_global%vlaiw(idx_start_patch : idx_end_patch)  ! lai adjusted for snow depth for calculation of resistances
  canopy%wcint  => canopy_global%wcint(idx_start_patch : idx_end_patch)  ! canopy rainfall interception
  canopy%evapfbl=> canopy_global%evapfbl(idx_start_patch : idx_end_patch,:)
  canopy%rwater => canopy_global%rwater(idx_start_patch : idx_end_patch,:)

!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)
!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(canopy%tv)

!  IF (threadID == 1) THEN  ! delaying thread 1
!    jj = 0
!    DO ii = 1, 100000*1000
!      jj = ii + 1
!    END DO
!  END IF

  ! Point local met input variables to the appropriate subset of the global array:
  met%ca      => met_global%ca(idx_start_patch : idx_end_patch)
  met%year    => met_global%year(idx_start_patch : idx_end_patch)
  met%moy     => met_global%moy(idx_start_patch : idx_end_patch)
  met%doy     => met_global%doy(idx_start_patch : idx_end_patch)
  met%hod     => met_global%hod(idx_start_patch : idx_end_patch)
  met%fsd     => met_global%fsd(idx_start_patch : idx_end_patch)
  met%fld     => met_global%fld(idx_start_patch : idx_end_patch)
  met%precip  => met_global%precip(idx_start_patch : idx_end_patch)
  met%precip_s=> met_global%precip_s(idx_start_patch : idx_end_patch)
!  met%tc      => met_global%tc(idx_start_patch : idx_end_patch)
  met%tk      => met_global%tk(idx_start_patch : idx_end_patch)
  met%tvair   => met_global%tvair(idx_start_patch : idx_end_patch)
  met%tvrad   => met_global%tvrad(idx_start_patch : idx_end_patch)
  met%pmb     => met_global%pmb(idx_start_patch : idx_end_patch)
  met%ua      => met_global%ua(idx_start_patch : idx_end_patch)
  met%qv      => met_global%qv(idx_start_patch : idx_end_patch)
  met%qvair   => met_global%qvair(idx_start_patch : idx_end_patch)
  met%da      => met_global%da(idx_start_patch : idx_end_patch)
  met%dva     => met_global%dva(idx_start_patch : idx_end_patch)
  met%coszen  => met_global%coszen(idx_start_patch : idx_end_patch)

!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)
!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(met%fsd)

  ! Point local energy and water balance variables to the appropriate subset of the global array:
  bal%drybal     => bal_global%drybal(idx_start_patch : idx_end_patch)
  bal%ebal       => bal_global%ebal(idx_start_patch : idx_end_patch)
  bal%ebal_tot   => bal_global%ebal_tot(idx_start_patch : idx_end_patch)
  bal%evap_tot   => bal_global%evap_tot(idx_start_patch : idx_end_patch)
  bal%osnowd0    => bal_global%osnowd0(idx_start_patch : idx_end_patch)
  bal%precip_tot => bal_global%precip_tot(idx_start_patch : idx_end_patch)
  bal%rnoff_tot  => bal_global%rnoff_tot(idx_start_patch : idx_end_patch)
  bal%wbal       => bal_global%wbal(idx_start_patch : idx_end_patch)
  bal%wbal_tot   => bal_global%wbal_tot(idx_start_patch : idx_end_patch)
  bal%wbtot0     => bal_global%wbtot0(idx_start_patch : idx_end_patch)
  bal%wetbal     => bal_global%wetbal(idx_start_patch : idx_end_patch)

  ! Point local radiation variables to the appropriate subset of the global array:
  rad%albedo     => rad_global%albedo(idx_start_patch : idx_end_patch,:)
  rad%extkb      => rad_global%extkb(idx_start_patch : idx_end_patch)
  rad%extkd2     => rad_global%extkd2(idx_start_patch : idx_end_patch)
  rad%extkd      => rad_global%extkd(idx_start_patch : idx_end_patch)
  rad%flws       => rad_global%flws(idx_start_patch : idx_end_patch)
  rad%fvlai      => rad_global%fvlai(idx_start_patch : idx_end_patch,:)
  rad%gradis     => rad_global%gradis(idx_start_patch : idx_end_patch,:)
  rad%latitude   => rad_global%latitude(idx_start_patch : idx_end_patch)
  rad%lwabv      => rad_global%lwabv(idx_start_patch : idx_end_patch)
  rad%qcan       => rad_global%qcan(idx_start_patch : idx_end_patch,:,:)
  rad%qssabs     => rad_global%qssabs(idx_start_patch : idx_end_patch)
  rad%rhocdf     => rad_global%rhocdf(idx_start_patch : idx_end_patch,:)
  rad%rniso      => rad_global%rniso(idx_start_patch : idx_end_patch,:)
  rad%scalex     => rad_global%scalex(idx_start_patch : idx_end_patch,:)
  rad%transd     => rad_global%transd(idx_start_patch : idx_end_patch)
  rad%trad       => rad_global%trad(idx_start_patch : idx_end_patch)
  rad%reffdf     => rad_global%reffdf(idx_start_patch : idx_end_patch,:)
  rad%reffbm     => rad_global%reffbm(idx_start_patch : idx_end_patch,:)
  rad%extkbm     => rad_global%extkbm(idx_start_patch : idx_end_patch,:)
  rad%extkdm     => rad_global%extkdm(idx_start_patch : idx_end_patch,:)
  rad%fbeam      => rad_global%fbeam(idx_start_patch : idx_end_patch)
  rad%cexpkbm    => rad_global%cexpkbm(idx_start_patch : idx_end_patch,:)
  rad%cexpkdm    => rad_global%cexpkdm(idx_start_patch : idx_end_patch,:)

  ! Point local roughness variables to the appropriate subset of the global array:
  rough%coexp    => rough_global%coexp(idx_start_patch : idx_end_patch)
  rough%disp     => rough_global%disp(idx_start_patch : idx_end_patch)
  rough%hruff    => rough_global%hruff(idx_start_patch : idx_end_patch)
  rough%hruff_grmx => rough_global%hruff_grmx(idx_start_patch : idx_end_patch)
  rough%rt0us    => rough_global%rt0us(idx_start_patch : idx_end_patch)
  rough%rt1usa   => rough_global%rt1usa(idx_start_patch : idx_end_patch)
  rough%rt1usb   => rough_global%rt1usb(idx_start_patch : idx_end_patch)
  rough%rt1      => rough_global%rt1(idx_start_patch : idx_end_patch)
  rough%term2    => rough_global%term2(idx_start_patch : idx_end_patch)
  rough%term3    => rough_global%term3(idx_start_patch : idx_end_patch)
  rough%term5    => rough_global%term5(idx_start_patch : idx_end_patch)
  rough%term6    => rough_global%term6(idx_start_patch : idx_end_patch)
  rough%usuh     => rough_global%usuh(idx_start_patch : idx_end_patch)
  rough%za_uv    => rough_global%za_uv(idx_start_patch : idx_end_patch)
  rough%za_tq    => rough_global%za_tq(idx_start_patch : idx_end_patch)
  rough%z0m      => rough_global%z0m(idx_start_patch : idx_end_patch)
  rough%zref_uv  => rough_global%zref_uv(idx_start_patch : idx_end_patch)
  rough%zref_tq  => rough_global%zref_tq(idx_start_patch : idx_end_patch)
  rough%zruffs   => rough_global%zruffs(idx_start_patch : idx_end_patch)
  rough%z0soilsn => rough_global%z0soilsn(idx_start_patch : idx_end_patch)
  rough%z0soil   => rough_global%z0soil(idx_start_patch : idx_end_patch)

  ! Point local soil parameters to the appropriate subset of the global array:
  soil%albsoil   => soil_global%albsoil(idx_start_patch : idx_end_patch, :)
  soil%bch       => soil_global%bch(idx_start_patch : idx_end_patch)
  soil%c3        => soil_global%c3(idx_start_patch : idx_end_patch)
  soil%clay      => soil_global%clay(idx_start_patch : idx_end_patch)
  soil%cnsd      => soil_global%cnsd(idx_start_patch : idx_end_patch)
  soil%css       => soil_global%css(idx_start_patch : idx_end_patch)
  soil%hsbh      => soil_global%hsbh(idx_start_patch : idx_end_patch)
  soil%hyds      => soil_global%hyds(idx_start_patch : idx_end_patch)
  soil%i2bp3     => soil_global%i2bp3(idx_start_patch : idx_end_patch)
  soil%ibp2      => soil_global%ibp2(idx_start_patch : idx_end_patch)
  soil%isoilm    => soil_global%isoilm(idx_start_patch : idx_end_patch)
  soil%rhosoil   => soil_global%rhosoil(idx_start_patch : idx_end_patch)
  soil%rs20      => soil_global%rs20(idx_start_patch : idx_end_patch)
  soil%sand      => soil_global%sand(idx_start_patch : idx_end_patch)
  soil%sfc       => soil_global%sfc(idx_start_patch : idx_end_patch)
  soil%silt      => soil_global%silt(idx_start_patch : idx_end_patch)
  soil%ssat      => soil_global%ssat(idx_start_patch : idx_end_patch)
  soil%sucs      => soil_global%sucs(idx_start_patch : idx_end_patch)
  soil%swilt     => soil_global%swilt(idx_start_patch : idx_end_patch)
  soil%zse       =  soil_global%zse
  soil%zshh      =  soil_global%zshh

  ! Point local soil and snow variables to the appropriate subset of the global array:
  ssoil%albsoilsn => ssoil_global%albsoilsn(idx_start_patch : idx_end_patch,:)
  ssoil%cls       => ssoil_global%cls(idx_start_patch : idx_end_patch)
  ssoil%dfn_dtg   => ssoil_global%dfn_dtg(idx_start_patch : idx_end_patch)
  ssoil%dfh_dtg   => ssoil_global%dfh_dtg(idx_start_patch : idx_end_patch)
  ssoil%dfe_ddq   => ssoil_global%dfe_ddq(idx_start_patch : idx_end_patch)
  ssoil%ddq_dtg   => ssoil_global%ddq_dtg(idx_start_patch : idx_end_patch)
  ssoil%evapsn    => ssoil_global%evapsn(idx_start_patch : idx_end_patch)
  ssoil%fwtop     => ssoil_global%fwtop(idx_start_patch : idx_end_patch)
  ssoil%gammzz    => ssoil_global%gammzz(idx_start_patch : idx_end_patch,:)
  ssoil%isflag    => ssoil_global%isflag(idx_start_patch : idx_end_patch)
  ssoil%potev     => ssoil_global%potev(idx_start_patch : idx_end_patch)
  ssoil%runoff    => ssoil_global%runoff(idx_start_patch : idx_end_patch)
  ssoil%rnof1     => ssoil_global%rnof1(idx_start_patch : idx_end_patch)
  ssoil%rnof2     => ssoil_global%rnof2(idx_start_patch : idx_end_patch)
  ssoil%rtsoil    => ssoil_global%rtsoil(idx_start_patch : idx_end_patch)
  ssoil%sconds    => ssoil_global%sconds(idx_start_patch : idx_end_patch,:)
  ssoil%sdepth    => ssoil_global%sdepth(idx_start_patch : idx_end_patch,:)
  ssoil%smass     => ssoil_global%smass(idx_start_patch : idx_end_patch,:)
  ssoil%snage     => ssoil_global%snage(idx_start_patch : idx_end_patch)
  ssoil%snowd     => ssoil_global%snowd(idx_start_patch : idx_end_patch)
  ssoil%smelt     => ssoil_global%smelt(idx_start_patch : idx_end_patch)
  ssoil%osnowd    => ssoil_global%osnowd(idx_start_patch : idx_end_patch)
  ssoil%ssdn      => ssoil_global%ssdn(idx_start_patch : idx_end_patch,:)
  ssoil%ssdnn     => ssoil_global%ssdnn(idx_start_patch : idx_end_patch)
  ssoil%tgg       => ssoil_global%tgg(idx_start_patch : idx_end_patch,:)
  ssoil%tggsn     => ssoil_global%tggsn(idx_start_patch : idx_end_patch,:)
  ssoil%tss       => ssoil_global%tss(idx_start_patch : idx_end_patch)
  ssoil%wb        => ssoil_global%wb(idx_start_patch : idx_end_patch,:)
  ssoil%wbfice    => ssoil_global%wbfice(idx_start_patch : idx_end_patch,:)
  ssoil%wbice     => ssoil_global%wbice(idx_start_patch : idx_end_patch,:)
  ssoil%wblf      => ssoil_global%wblf(idx_start_patch : idx_end_patch,:)
  ssoil%wbtot     => ssoil_global%wbtot(idx_start_patch : idx_end_patch)
  ssoil%pwb_min   => ssoil_global%pwb_min(idx_start_patch : idx_end_patch)
  ssoil%wetfac    => ssoil_global%wetfac(idx_start_patch : idx_end_patch)
  ssoil%owetfac   => ssoil_global%owetfac(idx_start_patch : idx_end_patch)
!  IF (ktau==1 .AND. lg==23) WRITE(45,'(a11,150f6.3)') 'wb again = ', ssoil%wb(:,1)
!  IF (ktau==1 .AND. lg==23) WRITE(45,'(a11,150f6.3)') &
!    'wb again = ', ssoil_global%wb(idx_start_patch : idx_end_patch,1)

  ! Point local cumulative flux variables to the appropriate subset of the global array:
  sum_flux%sumpn   => sum_flux_global%sumpn(idx_start_patch : idx_end_patch)
  sum_flux%sumrp   => sum_flux_global%sumrp(idx_start_patch : idx_end_patch)
  sum_flux%sumrpw  => sum_flux_global%sumrpw(idx_start_patch : idx_end_patch)
  sum_flux%sumrpr  => sum_flux_global%sumrpr(idx_start_patch : idx_end_patch)
  sum_flux%sumrs   => sum_flux_global%sumrs(idx_start_patch : idx_end_patch)
  sum_flux%sumrd   => sum_flux_global%sumrd(idx_start_patch : idx_end_patch)
  sum_flux%dsumpn  => sum_flux_global%dsumpn(idx_start_patch : idx_end_patch)
  sum_flux%dsumrp  => sum_flux_global%dsumrp(idx_start_patch : idx_end_patch)
  sum_flux%dsumrs  => sum_flux_global%dsumrs(idx_start_patch : idx_end_patch)
  sum_flux%dsumrd  => sum_flux_global%dsumrd(idx_start_patch : idx_end_patch)
  sum_flux%sumxrp  => sum_flux_global%sumxrp(idx_start_patch : idx_end_patch)
  sum_flux%sumxrs  => sum_flux_global%sumxrs(idx_start_patch : idx_end_patch)

  ! Point local vegetation parameters to the appropriate subset of the global array:
  veg%canst1     => veg_global%canst1(idx_start_patch : idx_end_patch)
  veg%dleaf      => veg_global%dleaf(idx_start_patch : idx_end_patch)
  veg%ejmax      => veg_global%ejmax(idx_start_patch : idx_end_patch)
  veg%frac4      => veg_global%frac4(idx_start_patch : idx_end_patch)
  veg%froot      => veg_global%froot(idx_start_patch : idx_end_patch, :)
  veg%hc         => veg_global%hc(idx_start_patch : idx_end_patch)
  veg%iveg       => veg_global%iveg(idx_start_patch : idx_end_patch)
  veg%meth       => veg_global%meth(idx_start_patch : idx_end_patch)
  veg%rp20       => veg_global%rp20(idx_start_patch : idx_end_patch)
  veg%rpcoef     => veg_global%rpcoef(idx_start_patch : idx_end_patch)
  veg%shelrb     => veg_global%shelrb(idx_start_patch : idx_end_patch)
  veg%wai        => veg_global%wai(idx_start_patch : idx_end_patch)
  veg%vegcf      => veg_global%vegcf(idx_start_patch : idx_end_patch)
  veg%tminvj     => veg_global%tminvj(idx_start_patch : idx_end_patch)
  veg%tmaxvj     => veg_global%tmaxvj(idx_start_patch : idx_end_patch)
  veg%vbeta      => veg_global%vbeta(idx_start_patch : idx_end_patch)
  veg%xalbnir    => veg_global%xalbnir(idx_start_patch : idx_end_patch)
  veg%vcmax      => veg_global%vcmax(idx_start_patch : idx_end_patch)
  veg%vlai       => veg_global%vlai(idx_start_patch : idx_end_patch)
  veg%xfang      => veg_global%xfang(idx_start_patch : idx_end_patch)
  veg%extkn      => veg_global%extkn(idx_start_patch : idx_end_patch)
  veg%deciduous  => veg_global%deciduous(idx_start_patch : idx_end_patch)

  IF (icycle > 0) THEN
    ! Point local pool sizes to the appropriate subset of the global array:
    casapool%Clabile      => casapool_g%Clabile      (idx_start_patch : idx_end_patch)
    casapool%dClabiledt   => casapool_g%dClabiledt   (idx_start_patch : idx_end_patch)
    casapool%Cplant       => casapool_g%Cplant       (idx_start_patch : idx_end_patch, :)
    casapool%Nplant       => casapool_g%Nplant       (idx_start_patch : idx_end_patch, :)
    casapool%Pplant       => casapool_g%Pplant       (idx_start_patch : idx_end_patch, :)
    casapool%dCplantdt    => casapool_g%dCplantdt    (idx_start_patch : idx_end_patch, :)
    casapool%dNplantdt    => casapool_g%dNplantdt    (idx_start_patch : idx_end_patch, :)
    casapool%dPplantdt    => casapool_g%dPplantdt    (idx_start_patch : idx_end_patch, :)
    casapool%ratioNCplant => casapool_g%ratioNCplant (idx_start_patch : idx_end_patch, :)
    casapool%ratioPCplant => casapool_g%ratioPCplant (idx_start_patch : idx_end_patch, :)
    casapool%Nsoilmin     => casapool_g%Nsoilmin     (idx_start_patch : idx_end_patch)
    casapool%Psoillab     => casapool_g%Psoillab     (idx_start_patch : idx_end_patch)
    casapool%Psoilsorb    => casapool_g%Psoilsorb    (idx_start_patch : idx_end_patch)
    casapool%Psoilocc     => casapool_g%Psoilocc     (idx_start_patch : idx_end_patch)
    casapool%dNsoilmindt  => casapool_g%dNsoilmindt  (idx_start_patch : idx_end_patch)
    casapool%dPsoillabdt  => casapool_g%dPsoillabdt  (idx_start_patch : idx_end_patch)
    casapool%dPsoilsorbdt => casapool_g%dPsoilsorbdt (idx_start_patch : idx_end_patch)
    casapool%dPsoiloccdt  => casapool_g%dPsoiloccdt  (idx_start_patch : idx_end_patch)
    casapool%Clitter      => casapool_g%Clitter      (idx_start_patch : idx_end_patch, :)
    casapool%Nlitter      => casapool_g%Nlitter      (idx_start_patch : idx_end_patch, :)
    casapool%Plitter      => casapool_g%Plitter      (idx_start_patch : idx_end_patch, :)
    casapool%dClitterdt   => casapool_g%dClitterdt   (idx_start_patch : idx_end_patch, :)
    casapool%dNlitterdt   => casapool_g%dNlitterdt   (idx_start_patch : idx_end_patch, :)
    casapool%dPlitterdt   => casapool_g%dPlitterdt   (idx_start_patch : idx_end_patch, :)
    casapool%ratioNClitter=> casapool_g%ratioNClitter(idx_start_patch : idx_end_patch, :)
    casapool%ratioPClitter=> casapool_g%ratioPClitter(idx_start_patch : idx_end_patch, :)
    casapool%Csoil        => casapool_g%Csoil        (idx_start_patch : idx_end_patch, :)
    casapool%Nsoil        => casapool_g%Nsoil        (idx_start_patch : idx_end_patch, :)
    casapool%Psoil        => casapool_g%Psoil        (idx_start_patch : idx_end_patch, :)
    casapool%dCsoildt     => casapool_g%dCsoildt     (idx_start_patch : idx_end_patch, :)
    casapool%dNsoildt     => casapool_g%dNsoildt     (idx_start_patch : idx_end_patch, :)
    casapool%dPsoildt     => casapool_g%dPsoildt     (idx_start_patch : idx_end_patch, :)
    casapool%ratioNCsoil  => casapool_g%ratioNCsoil  (idx_start_patch : idx_end_patch, :)
    casapool%ratioPCsoil  => casapool_g%ratioPCsoil  (idx_start_patch : idx_end_patch, :)

    ! Point local casaflux to the appropriate subset of the global array:
    casaflux%Cgpp         => casaflux_g%Cgpp         (idx_start_patch : idx_end_patch)
    casaflux%Cnpp         => casaflux_g%Cnpp         (idx_start_patch : idx_end_patch)
    casaflux%Crp          => casaflux_g%Crp          (idx_start_patch : idx_end_patch)
    casaflux%Crgplant     => casaflux_g%Crgplant     (idx_start_patch : idx_end_patch)
    casaflux%Nminfix      => casaflux_g%Nminfix      (idx_start_patch : idx_end_patch)
    casaflux%Nminuptake   => casaflux_g%Nminuptake   (idx_start_patch : idx_end_patch)
    casaflux%Plabuptake   => casaflux_g%Plabuptake   (idx_start_patch : idx_end_patch)
    casaflux%klabile      => casaflux_g%klabile      (idx_start_patch : idx_end_patch)
    casaflux%fracClabile  => casaflux_g%fracClabile  (idx_start_patch : idx_end_patch)
    casaflux%fracCalloc   => casaflux_g%fracCalloc   (idx_start_patch : idx_end_patch, :)
    casaflux%fracNalloc   => casaflux_g%fracNalloc   (idx_start_patch : idx_end_patch, :)
    casaflux%fracPalloc   => casaflux_g%fracPalloc   (idx_start_patch : idx_end_patch, :)
    casaflux%kplant       => casaflux_g%kplant       (idx_start_patch : idx_end_patch, :)
    casaflux%Crmplant     => casaflux_g%Crmplant     (idx_start_patch : idx_end_patch, :)
    casaflux%fromPtoL     => casaflux_g%fromPtoL     (idx_start_patch : idx_end_patch, :, :)
    casaflux%Cnep         => casaflux_g%Cnep         (idx_start_patch : idx_end_patch)
    casaflux%Crsoil       => casaflux_g%Crsoil       (idx_start_patch : idx_end_patch)
    casaflux%Nmindep      => casaflux_g%Nmindep      (idx_start_patch : idx_end_patch)
    casaflux%Nminloss     => casaflux_g%Nminloss     (idx_start_patch : idx_end_patch)
    casaflux%Nminleach    => casaflux_g%Nminleach    (idx_start_patch : idx_end_patch)
    casaflux%Nupland      => casaflux_g%Nupland      (idx_start_patch : idx_end_patch)
    casaflux%Nlittermin   => casaflux_g%Nlittermin   (idx_start_patch : idx_end_patch)
    casaflux%Nsmin        => casaflux_g%Nsmin        (idx_start_patch : idx_end_patch)
    casaflux%Nsimm        => casaflux_g%Nsimm        (idx_start_patch : idx_end_patch)
    casaflux%Nsnet        => casaflux_g%Nsnet        (idx_start_patch : idx_end_patch)
    casaflux%fNminloss    => casaflux_g%fNminloss    (idx_start_patch : idx_end_patch)
    casaflux%fNminleach   => casaflux_g%fNminleach   (idx_start_patch : idx_end_patch)
    casaflux%Pdep         => casaflux_g%Pdep         (idx_start_patch : idx_end_patch)
    casaflux%Pwea         => casaflux_g%Pwea         (idx_start_patch : idx_end_patch)
    casaflux%Pleach       => casaflux_g%Pleach       (idx_start_patch : idx_end_patch)
    casaflux%Ploss        => casaflux_g%Ploss        (idx_start_patch : idx_end_patch)
    casaflux%Pupland      => casaflux_g%Pupland      (idx_start_patch : idx_end_patch)
    casaflux%Plittermin   => casaflux_g%Plittermin   (idx_start_patch : idx_end_patch)
    casaflux%Psmin        => casaflux_g%Psmin        (idx_start_patch : idx_end_patch)
    casaflux%Psimm        => casaflux_g%Psimm        (idx_start_patch : idx_end_patch)
    casaflux%Psnet        => casaflux_g%Psnet        (idx_start_patch : idx_end_patch)
    casaflux%fPleach      => casaflux_g%fPleach      (idx_start_patch : idx_end_patch)
    casaflux%kplab        => casaflux_g%kplab        (idx_start_patch : idx_end_patch)
    casaflux%kpsorb       => casaflux_g%kpsorb       (idx_start_patch : idx_end_patch)
    casaflux%kpocc        => casaflux_g%kpocc        (idx_start_patch : idx_end_patch)
    casaflux%kmlabP       => casaflux_g%kmlabP       (idx_start_patch : idx_end_patch)
    casaflux%Psorbmax     => casaflux_g%Psorbmax     (idx_start_patch : idx_end_patch)
    casaflux%klitter      => casaflux_g%klitter      (idx_start_patch : idx_end_patch, :)
    casaflux%ksoil        => casaflux_g%ksoil        (idx_start_patch : idx_end_patch, :)
    casaflux%fromLtoS     => casaflux_g%fromLtoS     (idx_start_patch : idx_end_patch, :, :)
    casaflux%fromStoS     => casaflux_g%fromStoS     (idx_start_patch : idx_end_patch, :, :)
    casaflux%fromLtoCO2   => casaflux_g%fromLtoCO2   (idx_start_patch : idx_end_patch, :)
    casaflux%fromStoCO2   => casaflux_g%fromStoCO2   (idx_start_patch : idx_end_patch, :)
    casaflux%FluxCtolitter=> casaflux_g%FluxCtolitter(idx_start_patch : idx_end_patch, :)
    casaflux%FluxNtolitter=> casaflux_g%FluxNtolitter(idx_start_patch : idx_end_patch, :)
    casaflux%FluxPtolitter=> casaflux_g%FluxPtolitter(idx_start_patch : idx_end_patch, :)
    casaflux%FluxCtosoil  => casaflux_g%FluxCtosoil  (idx_start_patch : idx_end_patch, :)
    casaflux%FluxNtosoil  => casaflux_g%FluxNtosoil  (idx_start_patch : idx_end_patch, :)
    casaflux%FluxPtosoil  => casaflux_g%FluxPtosoil  (idx_start_patch : idx_end_patch, :)
    casaflux%FluxCtoco2   => casaflux_g%FluxCtoco2   (idx_start_patch : idx_end_patch)

    ! Point local casamet to the appropriate subset of the global array:
    casamet%glai     => casamet_g%glai    (idx_start_patch : idx_end_patch)
    casamet%lnonwood => casamet_g%lnonwood(idx_start_patch : idx_end_patch)
    casamet%Tairk    => casamet_g%Tairk   (idx_start_patch : idx_end_patch)
    casamet%precip   => casamet_g%precip  (idx_start_patch : idx_end_patch)
    casamet%tsoilavg => casamet_g%tsoilavg(idx_start_patch : idx_end_patch)
    casamet%moistavg => casamet_g%moistavg(idx_start_patch : idx_end_patch)
    casamet%btran    => casamet_g%btran   (idx_start_patch : idx_end_patch)
    casamet%Tsoil    => casamet_g%Tsoil   (idx_start_patch : idx_end_patch, :)
    casamet%moist    => casamet_g%moist   (idx_start_patch : idx_end_patch, :)
    casamet%iveg2    => casamet_g%iveg2   (idx_start_patch : idx_end_patch)
    casamet%ijgcm    => casamet_g%ijgcm   (idx_start_patch : idx_end_patch)
    casamet%isorder  => casamet_g%isorder (idx_start_patch : idx_end_patch)
    casamet%lat      => casamet_g%lat     (idx_start_patch : idx_end_patch)
    casamet%lon      => casamet_g%lon     (idx_start_patch : idx_end_patch)
    casamet%areacell => casamet_g%areacell(idx_start_patch : idx_end_patch)

    ! Point local casabal to the appropriate subset of the global array:
    casabal%FCgppyear    => casabal_g%FCgppyear    (idx_start_patch : idx_end_patch)
    casabal%FCnppyear    => casabal_g%FCnppyear    (idx_start_patch : idx_end_patch)
    casabal%FCrpyear     => casabal_g%FCrpyear     (idx_start_patch : idx_end_patch)
    casabal%FCrsyear     => casabal_g%FCrsyear     (idx_start_patch : idx_end_patch)
    casabal%FCneeyear    => casabal_g%FCneeyear    (idx_start_patch : idx_end_patch)
    casabal%FNdepyear    => casabal_g%FNdepyear    (idx_start_patch : idx_end_patch)
    casabal%FNfixyear    => casabal_g%FNfixyear    (idx_start_patch : idx_end_patch)
    casabal%FNsnetyear   => casabal_g%FNsnetyear   (idx_start_patch : idx_end_patch)
    casabal%FNupyear     => casabal_g%FNupyear     (idx_start_patch : idx_end_patch)
    casabal%FNleachyear  => casabal_g%FNleachyear  (idx_start_patch : idx_end_patch)
    casabal%FNlossyear   => casabal_g%FNlossyear   (idx_start_patch : idx_end_patch)
    casabal%FPweayear    => casabal_g%FPweayear    (idx_start_patch : idx_end_patch)
    casabal%FPdustyear   => casabal_g%FPdustyear   (idx_start_patch : idx_end_patch)
    casabal%FPsnetyear   => casabal_g%FPsnetyear   (idx_start_patch : idx_end_patch)
    casabal%FPupyear     => casabal_g%FPupyear     (idx_start_patch : idx_end_patch)
    casabal%FPleachyear  => casabal_g%FPleachyear  (idx_start_patch : idx_end_patch)
    casabal%FPlossyear   => casabal_g%FPlossyear   (idx_start_patch : idx_end_patch)
    casabal%glaimon      => casabal_g%glaimon      (idx_start_patch : idx_end_patch, :)
    casabal%glaimonx     => casabal_g%glaimonx     (idx_start_patch : idx_end_patch, :)
    casabal%cplantlast   => casabal_g%cplantlast   (idx_start_patch : idx_end_patch, :)
    casabal%nplantlast   => casabal_g%nplantlast   (idx_start_patch : idx_end_patch, :)
    casabal%pplantlast   => casabal_g%pplantlast   (idx_start_patch : idx_end_patch, :)
    casabal%clitterlast  => casabal_g%clitterlast  (idx_start_patch : idx_end_patch, :)
    casabal%nlitterlast  => casabal_g%nlitterlast  (idx_start_patch : idx_end_patch, :)
    casabal%plitterlast  => casabal_g%plitterlast  (idx_start_patch : idx_end_patch, :)
    casabal%csoillast    => casabal_g%csoillast    (idx_start_patch : idx_end_patch, :)
    casabal%nsoillast    => casabal_g%nsoillast    (idx_start_patch : idx_end_patch, :)
    casabal%psoillast    => casabal_g%psoillast    (idx_start_patch : idx_end_patch, :)
    casabal%nsoilminlast => casabal_g%nsoilminlast (idx_start_patch : idx_end_patch)
    casabal%psoillablast => casabal_g%psoillablast (idx_start_patch : idx_end_patch)
    casabal%psoilsorblast=> casabal_g%psoilsorblast(idx_start_patch : idx_end_patch)
    casabal%psoilocclast => casabal_g%psoilocclast (idx_start_patch : idx_end_patch)
    casabal%cbalance     => casabal_g%cbalance     (idx_start_patch : idx_end_patch)
    casabal%nbalance     => casabal_g%nbalance     (idx_start_patch : idx_end_patch)
    casabal%pbalance     => casabal_g%pbalance     (idx_start_patch : idx_end_patch)
    casabal%sumcbal      => casabal_g%sumcbal      (idx_start_patch : idx_end_patch)
    casabal%sumnbal      => casabal_g%sumnbal      (idx_start_patch : idx_end_patch)
    casabal%sumpbal      => casabal_g%sumpbal      (idx_start_patch : idx_end_patch)
    casabal%clabilelast  => casabal_g%clabilelast  (idx_start_patch : idx_end_patch)

    phen%phase    => phen_global%phase   (idx_start_patch : idx_end_patch)
    phen%doyphase => phen_global%doyphase(idx_start_patch : idx_end_patch, :)
    phen%TKshed   => phen_global%TKshed    ! This variable is for veg types
  END IF  ! icycle > 0

  ! patch array mapped
  patch  => patch_global(idx_start_patch : idx_end_patch)

! removed this because it has not been declared threadprivate and got mixed up
! by other threads in parallel runs (BP jun2010)
!  ! LAI array mapped
!  defaultLAI => defaultLAI_global(idx_start_patch : idx_end_patch,:)

!!!###   WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)
!!!###   WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(defaultLAI,1)
!!!###   IF (threadID == 4) THEN  ! delaying thread 1
!!!###     jj = 0
!!!###     DO ii = 1, 100000*1000
!!!###       jj = ii + 1
!!!###     END DO
!!!###   END IF

!!!### !$OMP BARRIER

  ! 6. 
  ! met data values are mapped from Mk3L (met%qvair and  met%da and met%dva
  ! are left untouched)
  ! * Note that Mk3L is sending values in 1xln2 dimension (i.e., MAP FROM
  !   1xln2 -> 1xmp)
  ! * Check for landpoints in the execution set for each points (of ln2 = 128)
  !   map the land point values into the pointer array inorder
  ! * Each parameter value from Mk3L will be mapped into elements equals to
  !   the max veg patches for a certain point in CABLE
  ! * convert the 64-bit Mk3L variables to 32-bit CABLE variables  
  ! * call get_defualt_LAI for each landpoint in the Mk3L set
  tmp = 1  ! init landpt counter
  DO i = 1, ln2  ! loop for each point in an execution set 
     IF (imsl(i,lg) .eq. 4) THEN   ! check for land point for that particular set (i.e., lg)
        met%fsd(landpt(tmp)%cstart:landpt(tmp)%cend)    = REAL(sg(i),r_1)/(1-SUM(     &
             & (rad%albedo(landpt(tmp)%cstart:landpt(tmp)%cend,1) & 
             & +rad%albedo(landpt(tmp)%cstart:landpt(tmp)%cend,2))*0.5* &
             & patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac))  ! surface downward shortwave
        met%fld(landpt(tmp)%cstart:landpt(tmp)%cend)    =  REAL(-rgsav(i,lg),r_1)  ! longwave   ! defined in MASIV3.f 
        met%precip(landpt(tmp)%cstart:landpt(tmp)%cend) =  MAX(REAL(rcondx(i),r_1),0.0) ! precipitation in mms of water per step
        met%pmb(landpt(tmp)%cstart:landpt(tmp)%cend)    =  REAL(pg(i),r_1)     ! surface pressure
        met%tk(landpt(tmp)%cstart:landpt(tmp)%cend)     =  REAL(ttg(i,1),r_1)    ! surface air temp. in K
!        met%tc(landpt(tmp)%cstart:landpt(tmp)%cend)     =  REAL(ttg(i,1),r_1)-273.15 ! surface air temp. in C
        met%ua(landpt(tmp)%cstart:landpt(tmp)%cend)     =  max(REAL(vmod(i),r_1),1.0) ! surface wind speed (m/s)
        met%qv(landpt(tmp)%cstart:landpt(tmp)%cend)     =   REAL(qtg(i,1),r_1)       ! specific humidity in the air lowest layer
        met%doy(landpt(tmp)%cstart:landpt(tmp)%cend)    =  kdays               ! local time day of year = days since 0 hr 1st Jan  
        IF(qcloud) THEN
           met%precip_s(landpt(tmp)%cstart:landpt(tmp)%cend) = REAL(rpreci(i),r_1)
        ELSE 
           met%precip_s(landpt(tmp)%cstart:landpt(tmp)%cend) = 0.0  
        END IF
        rough%za_uv(landpt(tmp)%cstart:landpt(tmp)%cend) = REAL(zmin(i),r_1)
        rough%za_tq(landpt(tmp)%cstart:landpt(tmp)%cend) = REAL(zmin(i),r_1)
        ! compute local hour of day from mins 
        met%hod(landpt(tmp)%cstart:landpt(tmp)%cend)    =  mod(mins, 1440)/60 + (longitude(idx_start(lg)+tmp-1)/180.0)*12.0
        met%year(landpt(tmp)%cstart:landpt(tmp)%cend)   =  iyear               ! local time year AD
        met%moy(landpt(tmp)%cstart:landpt(tmp)%cend)    =  month               ! local time month of year
        met%ca(landpt(tmp)%cstart:landpt(tmp)%cend)     =  REAL(rrvco2_global,r_1)  ! CO2 concentration 
        met%tvair(landpt(tmp)%cstart:landpt(tmp)%cend)  =  met%tk(landpt(tmp)%cstart:landpt(tmp)%cend)         ! it will be updated in CBM 
        met%tvrad(landpt(tmp)%cstart:landpt(tmp)%cend)  =  met%tk(landpt(tmp)%cstart:landpt(tmp)%cend)         ! it will be updated in CBM
!        !! These albedo values are from original Mk3L monthly file (albedo.nc). Added by BP nov 2009
!        IF (als(i) > 0.7) THEN
!          ! snow grid
!          soil%albsoil(landpt(tmp)%cstart:landpt(tmp)%cend,1)  = 0.9
!          soil%albsoil(landpt(tmp)%cstart:landpt(tmp)%cend,2)  = 0.6
!        ELSE
!          soil%albsoil(landpt(tmp)%cstart:landpt(tmp)%cend,1)  = (2.0/3.0) * REAL(als(i),r_1) ! monthly-varying albedo
!          soil%albsoil(landpt(tmp)%cstart:landpt(tmp)%cend,2)  = (4.0/3.0) * REAL(als(i),r_1) 
!        END IF
!        soil%albsoil(landpt(tmp)%cstart:landpt(tmp)%cend,3) = 0.0

        tmp = tmp + 1                         ! increment the land counter
     END IF
  END DO

!  !! These albedo values are from original Mk3L monthly file (albedo.nc). Added by BP nov 2009
!  ssoil%albsoilsn = soil%albsoil     ! initialization
!  IF (lg == 23) THEN
!    WRITE(46,*) 'albsoil Before'
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,1)
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,2)
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,3)
!    WRITE(46,*) 'albsoilsn Before'
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,1)
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,2)
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,3)
!  END IF

  ! Reset local time after longitude adjustment from GMT if negative hour of day
  WHERE(met%hod < 0.0)
     met%hod = met%hod + 24.0
  END WHERE

  ! computed from CABLE values
  met%coszen =  MAX(sinbet(met%doy,rad%latitude, met%hod),1.e-5)

  ! Set LAI of each patch to the present month
  veg%vlai(:) = defaultLAI_global(idx_start_patch : idx_end_patch, &
                                  met%moy(landpt(1)%cstart))
! removed this because it has not been declared threadprivate and got mixed up
! by other threads in parallel runs (BP jun2010)
!  veg%vlai(:) = defaultLAI(:,met%moy(landpt(1)%cstart)) 
!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)
!  WRITE(42,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(defaultLAI,1)

  ! Make sure LAI > 0  
  WHERE(veg%vlai<0.0001)    ! (BP jul2010)
    veg%vlai=0.0001
  END WHERE
! modified because this mess up the veg threshold of 0.01 LAI (BP jul2010)
!  WHERE(veg%vlai<0.05)
!     veg%vlai=0.05
!  END WHERE
  ! Further make sure ice grids do not have veg (BP jul2010)
  WHERE(veg%iveg==15)       ! ice for IGBP=15 
    veg%vlai=MIN(0.008,veg%vlai)
  END WHERE

!!bp  ! Make evergreen needle leaf forest at least 3.0 for LAI (quick fix by BP nov2009)
  WHERE(veg%iveg==1 .AND. veg%vlai<3.0)  veg%vlai = 3.0
!!bp  ! the following LAI changes are for testing, should delete afterwards (BP)
!!bp  veg_global%vlai(690) = 0.001
!!bp  veg_global%vlai(718) = 0.001
!!bp  veg_global%vlai(1586) = 0.001
!!bp  veg_global%vlai(1672) = 0.001

  ! 7. calling cable 
  !    ktau = nstepsa + 1, kstart = 2 variable not necessary.
  !    cable local arrays are passed
  !    nvegst and mstype are loaded using load parameters.
!  IF (lg == 23) THEN
!    i = 1039
!    j = 1672
!    write(39,'(a7,f7.2,i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      'Before ', longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy%delwc(j), ssoil%rnof1(j), ssoil%rnof2(j)
!  ENDIF
!  IF (lg == 18) THEN
!    i = 852
!    j = 1306
!    write(54,'(a7,f7.2,i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      'Before ', longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy%delwc(j), ssoil%rnof1(j), ssoil%rnof2(j)
!  ENDIF
!  IF (ktau >= 11904) write(46,*) 'before  : (lg,ktau,mp = )', lg, ktau, mp
!  IF (ktau == 11905 .AND. lg == 18) THEN
!    WRITE(47,'(6i9)') threadID,lg,idx_start_patch,idx_end_patch,(idx_end_patch-idx_start_patch+1),size(air%rho)
!    WRITE(47,*) 'nap', landpt(:)%nap
!    WRITE(47,*) 'tk(lg=18):', met%tk
!    WRITE(47,*) 'patch 1306:', patch_global(1306)%latitude, patch_global(1306)%longitude, patch_global(1306)%frac, met_global%tk(1306)
!  END IF
!  IF (ktau == 1 .AND. lg == 1) THEN
!    WRITE(47,*) 'lat 850-855:', latitude(850:855)
!    WRITE(47,*) 'lon 850-855:', longitude(850:855)
!  END IF
  DO ii = 1, mp
    IF (met%tk(ii) > 343.0) WRITE(47,*) 'extra high global tk:', ktau, lg, ii, patch(ii)%frac, patch(ii)%latitude, patch(ii)%longitude, met%tk(ii), canopy%vlaiw(ii)
  END DO

  CALL cbm(ktau,kstart,kend,dels,air,bgc,canopy,met, &
       bal,rad,rough,soil,ssoil,sum_flux,veg)
!  CALL cbm(ktau,kstart,kend,dels,air,bgc,canopy,met, &
!       bal,rad,rough,soil,ssoil,sum_flux,veg,mvtype,mstype)  
  IF (icycle > 0) CALL bgcdriver(ktau,kstart,kend,dels,met,ssoil,canopy, &
       veg,soil,casabiome,casapool,casaflux,casamet,casabal,phen)
  CALL sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
       soil, ssoil, sum_flux, veg, met, casaflux)
!  CALL sumcflux(ktau, kstart, kend, dels, mvtype, mstype, bgc, canopy,  &
!       soil, ssoil, sum_flux, veg, met, casaflux)
!  IF (lg == 18) THEN
!    i = 852
!    j = 1306
!    write(54,'(a7,f7.2,i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      'After  ', longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy%delwc(j), ssoil%rnof1(j), ssoil%rnof2(j)
!  ENDIF
!
!  IF (lg == 23) THEN
!    WRITE(46,*) 'albsoil After'
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,1)
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,2)
!    WRITE(46,'(i5,200f5.2)') ktau, soil%albsoil(:,3)
!    WRITE(46,*) 'albsoilsn After'
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,1)
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,2)
!    WRITE(46,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,3)
!  END IF

!   print *,"Testing AFTER CBM"

  ! 8.  
  ! Compute output variables from CABLE to Mk3L
  ! map the values from the CABLE array to the Mk3L global size array (i.e.,
  ! 1xmp*max_vegpatches -> 1xln2)
  ! the average values using the patchfrac for each variable from CABLE for
  ! each landpoint is mapped to Mk3L variable 
  !  - 32 bit values from CABLE are mapped to 64 bit values for Mk3L
  scalev = 0.0
  totpev = 0.0 
  tddd = 0.0 
  runoff = 0.0
  perc = 0.0
  wb = 0.0 
  tmp = 1 !  init counter for landpt
  DO i = 1, ln2  ! loop for each point in an execution set
     IF (imsl(i,lg) .eq. 4) THEN   ! check for land point for that particular set (i.e., lg)
        cabTscrn(i)=REAL(SUM(canopy%tscrn(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2) + 273.16
        ! variables similar to surfa  
        rg(i)=-REAL(SUM(met%fld(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac)-  &
             sboltz*emleaf* &
             SUM(canopy%tv(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac)**4*           &
             (1-SUM(rad%transd(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac))           &
             -SUM(rad%flws(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac)*               &
             SUM(rad%transd(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)         ! net longwave heating at ground
        totpev(i)=REAL((SUM(ssoil%potev(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac)+           &
             SUM(canopy%potev_c(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac))*dels/  &
             SUM(air%rlam(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2) ! potential evaporation in mm/timestep.
        scalev(i)=totpev(i)                      ! scaling evaporation
        tg(i)     =  REAL(SUM(rad%trad(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! combined soil-canopy temperature of the surface layer
        fg(i)     =  REAL(SUM(canopy%fh(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)    ! sensible heat flux
        eg(i)     =  REAL(SUM(canopy%fe(landpt(tmp)%cstart:landpt(tmp)%cend)*patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)    ! latent heat flux                
        tb2(i)    =  REAL(SUM(ssoil%tgg(landpt(tmp)%cstart:landpt(tmp)%cend,2)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)  ! temperature of the second layer 
        tb3(i)    =  REAL(SUM(ssoil%tgg(landpt(tmp)%cstart:landpt(tmp)%cend,6)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)  ! temperature of the lowest layer      
        mcmax(i)  =  0.2  ! this is set to a constant and not been used anywhere in Mk3L: maximum skin reservoir depth
        tddd(i)   =  REAL(SUM(canopy%wcint(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)    ! canopy interception
        tgg(i) =  REAL(SUM(ssoil%tgg(landpt(tmp)%cstart:landpt(tmp)%cend,1)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2) ! not really used in Mk3L
        tgf(i) =  REAL(SUM(canopy%tv(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)       ! vegetation covered ground temperature
        mc(i)  =  REAL(SUM(canopy%cansto(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)       ! canopy water storage (mm)
        osnowd(i) = REAL(SUM(ssoil%osnowd(landpt(tmp)%cstart:landpt(tmp)%cend)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! snow depth from previous time step
        snage(i)  = REAL(SUM(ssoil%snage(landpt(tmp)%cstart:landpt(tmp)%cend)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! snow age
        ssdnn(i)  = REAL(SUM(ssoil%ssdnn(landpt(tmp)%cstart:landpt(tmp)%cend)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! average snow density
        gflux(i)  = REAL(SUM(canopy%ghflux(landpt(tmp)%cstart:landpt(tmp)%cend)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! ground heat flux (W/m2)
        sgflux(i)  = REAL(SUM(canopy%sghflux(landpt(tmp)%cstart:landpt(tmp)%cend)* &
                patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! ground heat flux (W/m2)
        IF (SUM(ssoil%isflag(landpt(tmp)%cstart:landpt(tmp)%cend)) > 0) THEN
          isflag(i) = 1  ! have snow
        ELSE
          isflag(i) = 0  ! no snow
        END IF
        DO j=1,ms
          wb(i,j) = REAL(SUM(ssoil%wb(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)
          wbice(i,j) = REAL(SUM(ssoil%wbice(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)
          tggsl(i,j) =  REAL(SUM(ssoil%tgg(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! ground temperature at the previous time step     
        END DO
        DO j=1,3
          tggsn(i,j) =  REAL(SUM(ssoil%tggsn(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! snow layer temperature
          smass(i,j) =  REAL(SUM(ssoil%smass(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! snow mass
          ssdn3(i,j) =  REAL(SUM(ssoil%ssdn(landpt(tmp)%cstart:landpt(tmp)%cend,j)* &
               patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)     ! snow density
        END DO
        ! variables similar to surfb
        snowd(i) = REAL(SUM(ssoil%snowd(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)  ! snow depth
        runoff(i)= REAL(SUM(ssoil%runoff(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2) ! grnd runoff
        perc(i) = REAL(SUM(ssoil%rnof2(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2) ! percolation
        wg(i)     =  REAL(MAX(0.,wb(i,1)-SUM(soil%swilt(landpt(tmp)%cstart:landpt(tmp)%cend)* &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac))) ! soil moisture available for evaporation at the top layer
        ! computation for wg2: volume weighted soil moisture
        zse2 = soil%zse
        zse2 = zse2/SUM(zse2) ! Normalization
        ! the average values of soil%wb are used for wg2 computation 
        ! wg2(i)    =  SUM(ssoil%wb(j,:)*zse2)  ! 
        wg2(i) = SUM(wb(i,:)*zse2)
        als(i) = REAL(SUM((rad%albedo(landpt(tmp)%cstart:landpt(tmp)%cend,1)+ &
             rad%albedo(landpt(tmp)%cstart:landpt(tmp)%cend,2))*0.5*      &
             patch(landpt(tmp)%cstart:landpt(tmp)%cend)%frac),r_2)      ! albedo 
        tmp = tmp + 1                    ! landpoint counter incremented
     END IF
  END DO
!bp Output diagnostics at several points to check on CABLE
!  write(42,*) 'ktau, lg = ', ktau, lg
!  write(42,*) 'landpt(:)%nap, cstart, cend ='
!  write(42,'(99(i6))') landpt(:)%nap
!  write(42,'(99(i6))') landpt(:)%cstart
!  write(42,'(99(i6))') landpt(:)%cend
!  write(42,*) 'idx_start(lg), idx_end(lg), mp, idx_start_patch, idx_end_patch:'
!  write(42,*) idx_start(lg), idx_end(lg), mp, idx_start_patch, idx_end_patch
!   write(42,'(a9,130(e10.2))') 'totpev = ', totpev
!   write(42,'(a9,130(e10.2))') 'canopyPev', canopy%potev_c
!   write(42,'(a9,130(e10.2))') 'ssoilPev ', ssoil%potev
!  write(42,*) 'latitude: (lg,ktau = )', lg, ktau
!  write(42,'(130f7.1)') latitude(idx_start(lg) : idx_end(lg))
!  write(42,*) 'longitude: (lg,ktau = )', lg, ktau
!  write(42,'(130f7.1)') longitude(idx_start(lg) : idx_end(lg))
!  write(43,*) 'latitude: (lg,ktau = )', latitude(idx_start(lg)), ktau
!  write(43,'(a12,130i7)') '# patches:  ', landpt(:)%nap
!  write(43,'(a12,130f7.1)') 'longitudes: ', longitude(idx_start(lg):idx_end(lg))
!  write(44,*) 'latitude: (lg,ktau = )', latitude(idx_start(lg)), ktau
!  write(44,'(a12,130i7)') '# patches:  ', landpt(:)%nap
!  write(44,'(a12,130f7.1)') 'longitudes: ', longitude(idx_start(lg):idx_end(lg))
!  IF (ktau>=11905) write(46,*) 'latitude: (lg,ktau = )', lg,ktau,latitude(idx_start(lg))
!  write(47,*) 'latitude: (lg,ktau = )', latitude(idx_start(lg)), ktau
!  write(48,*) 'latitude: (lg,ktau = )', latitude(idx_start(lg)), ktau
!  IF (lg == 9) THEN
!    i = 526
!    j = 690
!    write(36,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
!    i = 540
!    j = 718
!    write(37,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
!    write(41,*) 'longitude', ktau
!    write(41,'(100f6.1)') longitude(idx_start(lg) : idx_end(lg))
!    write(41,*) 'nap'
!    write(41,'(100i6)') landpt(:)%nap
!    write(41,*) 'veg%iveg', ktau
!    write(41,'(100(i2,4x))') veg%iveg
!    write(41,*) 'veg%vlai', ktau
!    write(41,'(100f6.3)') veg%vlai
!    write(41,*) 'canopy%vlaiw', ktau
!    write(41,'(100f6.3)') canopy%vlaiw
!    write(41,*) 'sum(rad%fvlai,2)', ktau
!    write(41,'(100f6.3)') sum(rad%fvlai,2)
!  ENDIF
!  IF (lg == 21) THEN
!    i = 988
!    j = 1586
!    write(38,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
!!    write(38,*) 'soil type = ', soil_global%isoilm(j)
!  ENDIF
!  IF (lg == 23) THEN
!    i = 1039
!    j = 1672
!    write(39,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
!      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
!      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
!      (rad_global%fvlai(j,tmp),tmp=1,2), &
!      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
!      rad_global%extkb(j), rad_global%extkd(j), &
!      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
!      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
!      (rad_global%albedo(j,tmp),tmp=1,2), &
!      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
!      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
!      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
!      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
!      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
!      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
!      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
!      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
!      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
!      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
!      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
!      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
!!!    write(39,*) 'soil type = ', soil_global%isoilm(j)
!  ENDIF
  IF (lg == 22) THEN
    i = 1033
    j = 1627   ! Brazilian woody savannas
    write(51,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
      (rad_global%fvlai(j,tmp),tmp=1,2), &
      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
      rad_global%extkb(j), rad_global%extkd(j), &
      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
      (rad_global%albedo(j,tmp),tmp=1,2), &
      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
    j = 1629   ! Brazilian savannas
    write(52,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
      (rad_global%fvlai(j,tmp),tmp=1,2), &
      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
      rad_global%extkb(j), rad_global%extkd(j), &
      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
      (rad_global%albedo(j,tmp),tmp=1,2), &
      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
  ENDIF
  IF (lg == 28) THEN  ! Amazon evergreen broadleaf forest
    i = 1217
    j = 1985
    write(53,'(2(f7.2),i3,i8,7(f8.2),3(f12.2),24(f10.2),f10.1,43(f8.2))') &
      latitude(i), longitude(i), veg_global%iveg(j), ktau, met_global%tk(j), &
      met_global%tvair(j), canopy_global%tv(j), ssoil_global%tss(j), &
      (rad_global%fvlai(j,tmp),tmp=1,2), &
      canopy_global%vlaiw(j), met_global%fsd(j), met_global%fld(j), &
      rad_global%extkb(j), rad_global%extkd(j), &
      (rad_global%qcan(j,1,tmp),tmp=1,3), (rad_global%qcan(j,2,tmp),tmp=1,3), &
      rad_global%qssabs(j), (ssoil_global%albsoilsn(j,tmp),tmp=1,3), &
      (rad_global%albedo(j,tmp),tmp=1,2), &
      rad_global%lwabv(j), rad_global%transd(j), canopy_global%fevc(j), &
      canopy_global%fwet(j), canopy_global%fhv(j), canopy_global%fes(j), &
      canopy_global%fhs(j), canopy_global%fns(j), canopy_global%fnv(j), &
      rad_global%trad(j), bal_global%ebal(j), bal_global%ebal_tot(j), &
      rad_global%flws(j), canopy_global%fev(j), ssoil_global%cls(j), &
      canopy_global%fh(j), canopy_global%ga(j), ssoil_global%snowd(j), &
      (ssoil_global%tgg(j,tmp),tmp=1,6), (ssoil_global%tggsn(j,tmp),tmp=1,3), &
      (ssoil_global%wb(j,tmp),tmp=1,6), (ssoil_global%wbice(j,tmp),tmp=1,6), &
      ssoil_global%osnowd(j), ssoil_global%snage(j), ssoil_global%ssdnn(j), &
      (ssoil_global%ssdn(j,tmp),tmp=1,3), (ssoil_global%smass(j,tmp),tmp=1,3), &
      canopy_global%ghflux(j), canopy_global%sghflux(j), met_global%precip(j), &
      met_global%precip_s(j), canopy_global%delwc(j), ssoil_global%rnof1(j), ssoil_global%rnof2(j)
  ENDIF
  ! Deallocate landpt
   DEALLOCATE(landpt)

END SUBROUTINE CABLE_run  ! end of CABLE step

!###########################################################
!####      CABLE_write_output subroutine          ##########
!###########################################################
SUBROUTINE CABLE_write_output(nstepsa,mstep)
  USE CABLE_Mk3L_module
  INTEGER(i_d), INTENT(IN) :: nstepsa               ! the execution step counter
  INTEGER(i_d), INTENT(IN) :: mstep                 ! step size

  ! set the mp to mp_global so that all the local variables in the
  ! output routine size will match with global arrays
  mp = mp_global
  mland = mland_global
  ALLOCATE(landpt(mland_global))
  landpt = landpt_global
  patch => patch_global
  ktau = nstepsa + 1  ! set ktau to values local to current run
  kend = INT(365 * 24 * 60/mstep)
  CALL write_output(ktau,dels,met_global,canopy_global,ssoil_global, &
       rad_global,bal_global,air_global,soil_global,veg_global)  
  IF(icycle>0 .AND. ktau==kend) CALL casa_poolout(ktau,veg_global,soil_global, &
       casabiome,casapool_g,casaflux_g,casamet_g,casabal_g,phen_global)

  DEALLOCATE(landpt)

END SUBROUTINE CABLE_write_output ! end of CABLE_write_output

!###########################################################
!####             CABLE_close subroutine          ##########
!########################################################### 
  SUBROUTINE CABLE_close()
    USE CABLE_Mk3L_module
    ! 9. write global restart file
    ! Write restart file if requested
    ! necessary variables have to be reset to their global values
    mland = mland_global
    mp = mp_global  ! needed for close_output_file
    ALLOCATE(landpt(mland_global))  ! define the size to the local mland for each step
    landpt = landpt_global
    patch => patch_global ! isnt it the restart file only for landpts????
    IF(output%restart) CALL create_restart(  &
         logn,ktau,dels,soil_global,veg_global,ssoil_global, &
         canopy_global,rough_global,rad_global,bgc_global,bal_global)
!    IF(output%restart) CALL create_restart(  &
!         logn,ktau,dels,soil_global,veg_global,ssoil_global, &
!         canopy_global,rough_global,rad_global,bgc_global,bal_global,mvtype,mstype)
    
    ! close output file
    ! the global arrays are deallocated in this routine 
    CALL close_output_file(bal_global, air_global, bgc_global, canopy_global, &
             met_global, rad_global, rough_global, soil_global, ssoil_global, &
             sum_flux_global, veg_global)

    ! Close log file
    CLOSE(logn)
  END SUBROUTINE CABLE_close ! end of CABLE close
