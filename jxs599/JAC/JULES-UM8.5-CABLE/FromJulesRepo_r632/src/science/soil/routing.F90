#if !defined(UM_JULES)
!###############################################################################
!###############################################################################
! Code for runoff routing.
!###############################################################################
!###############################################################################

! subroutine route_drive
! Driver routine for runoff routing.

  SUBROUTINE route_drive( land_pts,sub_surf_roff,surf_roff,roffAccumLand )

  USE route_mod, ONLY :  &
!  imported scalar parameters
     routeTypeTRIP  &
!  imported scalars with intent(in)
    ,routeTimestep,routeType  &
!  imported scalars with intent(inout)
    ,routeCount,routeStep

  USE timeConst, ONLY : &
     iSecInDay

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Array arguments with intent(in)
  INTEGER, INTENT(in) :: land_pts              !  number of land points

! Array arguments with intent(in)
  REAL, INTENT(in) :: sub_surf_roff(land_pts)  !  Sub-surface runoff (kg m-2 s-1)
  REAL, INTENT(in) :: surf_roff(land_pts)      !  Surface runoff (kg m-2 s-1)

! Array arguments with intent(inout)
  REAL, INTENT(inout) :: roffAccumLand(land_pts) !  average runoff (production)
!                                    rate between calls to routing (kg m-2 s-1)

!-------------------------------------------------------------------------------

  IF ( routeCount == 0 ) THEN

    SELECT CASE ( routeType )
      CASE ( routeTypeTRIP )
!       Reset accumulated runoff to zero.
        roffAccumLand(:) = 0.0
      CASE default
        WRITE(*,*)'ERROR: route_drive: do not recognise routeType=',TRIM(routeType)
        STOP
    END SELECT

  ENDIF

! Increment counters.
  routeStep = routeStep + 1
  routeCount = routeCount + 1

!-------------------------------------------------------------------------------
! Accumulate runoff.
!-------------------------------------------------------------------------------
  SELECT CASE ( routeType )

    CASE ( routeTypeTRIP )

!     Accumulate the total runoff.
      roffAccumLand(:) = roffAccumLand(:) + ( surf_roff(:) + sub_surf_roff(:) )  &
                           / REAL(routeTimeStep)

    CASE default
      WRITE(*,*)'ERROR: route_drive: do not recognise routeType=',TRIM(routeType)
      STOP

  END SELECT

!-------------------------------------------------------------------------------
! Decide if routing model is to be called on this timestep and call.
!-------------------------------------------------------------------------------
  IF ( routeStep >= routeTimeStep ) THEN

!-------------------------------------------------------------------------------
!   Routing is to be called.
!-------------------------------------------------------------------------------

    SELECT CASE ( routeType )

      CASE ( routeTypeTRIP )
!-------------------------------------------------------------------------------
!       Correct average runoff rate if did not get expected number of increments.
!       This deals with awkward cases - e.g. when spin-up interupts.
!-------------------------------------------------------------------------------
        IF ( routeCount /= routeTimeStep ) THEN
          IF ( routeCount > 0 ) THEN
            roffAccumLand(:) = roffAccumLand(:) * REAL(routeTimeStep)  &
                       / REAL(routeCount)
          ELSE
!           At present I'm not sure if we'd ever get here! But just in case...
            WRITE(*,*) 'ERROR: route_drive div by routeCount=0'
            STOP
          ENDIF
        ENDIF

!-------------------------------------------------------------------------------
!       Call TRIP routing driver.
!-------------------------------------------------------------------------------
        CALL route_drive_trip( roffAccumLand )

      CASE default
!-------------------------------------------------------------------------------
!       Default case for routeType.
!-------------------------------------------------------------------------------

        WRITE(*,*)'ERROR: route_drive: do not recognise routeType=',TRIM(routeType)
        STOP

    END SELECT

!-------------------------------------------------------------------------------
!   Reset counters after a call to routing.
!-------------------------------------------------------------------------------
    routeStep = 0
    routeCount = 0

!  ELSE

!   Routing is not done on this timestep.
!   Set routing diagnostics to missing data. Prognostics will retain current value.
!    select case ( routeType )
!      case ( routeTypeTRIP )
!        if ( outPosRFlow > 0 ) outValExtra(outPosRFlow:outPosRFlow+npRoute-1) = undefOut
!        if ( outPosRRun > 0 ) outValExtra(outPosRRun:outPosRRun+npRoute-1) = undefOut
!      case default
!        write(*,*)'ERROR: route_drive: do not recognise routeType=',trim(routeType)
!        stop
!    end select

  ENDIF   !  routeStep

  END SUBROUTINE route_drive
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
! subroutine route_drive_trip
! Driver routine for runoff routing by the TRIP (linear) model.

  SUBROUTINE route_drive_trip( roffAccumLand )

!-------------------------------------------------------------------------------
! Modules used:

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts,nx=>row_length,ny=>rows  &
!  imported arrays with intent(in)
    ,land_index

  USE earth_utils, ONLY :  &
!  imported procdures
     earthArea

  USE grid_utils, ONLY : &
!  imported procedures
     getXYpos

  USE offline_diag, ONLY :  &
!  imported scalars with intent(in)
     useRflowDiag,useRrunDiag  &
!  imported arrays with intent(inout)
    ,rflowDiag,rrunDiag

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     npRoute,nxRoute,nyRoute,routeDlat,routeDlon,routeLat1,routeLon1  &
    ,routeRegLatLon,routeRegrid   &
!  imported arrays with intent(in)
    ,routeIndex  &
!  imported arrays with intent(inout)
    ,routeNext

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     nxGrid,nyGrid,regDlat,regDlon,regLat1,regLon1  &
!  imported arrays with intent(in)
    ,latitude,longitude
!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Array arguments with intent(in)

  REAL, INTENT(in) :: roffAccumLand(land_pts)  !  average rate of runoff since
!                                                   last routing call (kg m-2 s-1)

! Local scalar variables.

  INTEGER :: i      !  work
  INTEGER :: ip     !  loop counter
  INTEGER :: ix     !  work
  INTEGER :: iy     !  work
  INTEGER :: maxLin ! an array size

! Local array variables.

  INTEGER :: x(npRoute)      !  x index on the routing grid of active routing points
  INTEGER :: xNext(npRoute)  !  x location of the next downstream point
!                                 <0 indicates next point is off grid
  INTEGER :: y(npRoute)      !  y index on the routing grid of active routing points
  INTEGER :: yNext(npRoute)  !  y location of the next downstream point
!                                 If xNext(i)<0, yNext(i) gives the direction of the
!                                 next point, in terms of index in flowDirSet.

  REAL :: outFlow(nxRoute,nyRoute)     !  rate of channel flow leaving gridbox (kg s-1)
  REAL :: routeArea(nxRoute,nyRoute)   !  area of each gridbox on routing grid (m2)
  REAL :: routeLat(nxRoute,nyRoute)    !  latitude of points on routing grid (degrees)
  REAL :: routeLength(nxRoute,nyRoute) !  distance to next downstream gridpoint (m)
!                                           This is only calculated at "active" routing points.
  REAL :: routeLon(nxRoute,nyRoute)    !  longitude of points on routing grid (degrees)
  REAL :: runoff(nxRoute,nyRoute)      !  average rate of runoff, on routing grid (kg s-1)
!                                           During work stages, units may be kg m-2 s-1.
  REAL :: roffAccumGrid(nxGrid,nyGrid) !  average rate of runoff since last routing
!                                             call (kg m-2 s-1), on a grid that
!              is a superset of (can be equal to) the model grid

  LOGICAL :: routeMask(nxRoute,nyRoute) !  TRUE at land points on routing grid,
!                                          else FALSE

!-------------------------------------------------------------------------------
! Curent code can only deal with "regular" lat/lon grids.
!-------------------------------------------------------------------------------
  IF ( .NOT. routeRegLatLon ) THEN
    WRITE(*,*)'ERROR: route_drive_trip: code only exists for regular lat/lon grids.'
    WRITE(*,*)'Sorry. Maybe you''d like to write some new code...?'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Calculate location on routing grid of each land point, and find the point
! immediately downstream. Also set a mask to show land points on the routing
! grid.
!-------------------------------------------------------------------------------

! Initialise.
  xNext(:) = 0
  yNext(:) = 0
  routeMask(:,:) = .FALSE.

  DO ip=1,npRoute

!   Get location of this point.
    CALL getXYPos( routeIndex(ip),nxRoute,nyRoute,x(ip),y(ip) )

!   Set mask to TRUE at this land point.
    routeMask(x(ip),y(ip)) = .TRUE.

!   Get location of next downstream point.
    IF ( routeNext(ip) > 0 ) THEN
!     Downstream point is on the grid.
      CALL getXYPos( routeNext(ip),nxRoute,nyRoute,xNext(ip),yNext(ip) )
    ELSEIF ( routeNext(ip) == -9 ) THEN
!     There is no defined flow direction, and so no downstream point.
      xNext(ip) = -9
      yNext(ip) = -9
    ELSE
!     Downstream point is off edge of grid.
!     routeNext has been set to -1*index of direction in flowDirSet.
      xNext(ip) = -1 * routeNext(ip)
      yNext(ip) = xNext(ip)
    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
! Calculate latitude and longitude of each gridpoint on routing grid.
!-------------------------------------------------------------------------------
  DO iy=1,nyRoute
    DO ix=1,nxRoute
      routeLat(ix,iy) = routeLat1 + (iy-1)*routeDlat
      routeLon(ix,iy) = routeLon1 + (ix-1)*routeDlon
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Calculate distance between grid points.
!-------------------------------------------------------------------------------
  CALL getRouteLen( x,xNext,y,yNext,routeLat,routeLon,routeLength )

!-------------------------------------------------------------------------------
! Calculate area of each routing gridbox.
!-------------------------------------------------------------------------------
  DO iy=1,nyRoute
    DO ix=1,nxRoute
      routeArea(ix,iy) = earthArea(  &
                routeLat(ix,iy)-0.5*routeDlat,routeLat(ix,iy)+0.5*routeDlat  &
               ,routeLon(ix,iy)-0.5*routeDlon,routeLon(ix,iy)+0.5*routeDlon ) * 1.0e6
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Put accumulated runoff onto a grid.
! If the "main" model grid is 2-D, this target grid is the 2-D grid.
! If the "main" grid is a vector (in offline applications of JULES this is
! possible if points from a larger grid have been compressed - e.g. land points
! selected from a larger grid.), the target grid is the larger grid, across
! which the points are to be scattered.
!
! Arguably this would be better done at a higher level, since it's unlikely that
! the requirement to grid the runoff will be unique to TRIP. But it's here for
! now. An advantage of being here is that we only grid the fields we need for
! this particular model (if done higher, we would have to test which fields to
! grid for which model - not impossible!).
!
! Note that the regridded runoff may be non-zero at TRIP sea points.
!
!-------------------------------------------------------------------------------

  roffAccumGrid(:,:) = 0.0

  IF ( nxGrid /= nx .OR. nyGrid /= ny ) THEN

!   The "main" model grid is not a regular lat/lon grid, but is a subset of
!   such a regular grid.
!   We cannot use land_index to scatter the points across this grid. Instead,
!   use lat/lon to locate in the larger grid of shape (nxGrid,nyGrid).
!   Note that this code will not currently work with grids that are just not
!   a subset of a "regular" grid - e.g. rotated grid - since then nxGrid and
!   nyGrid have not been given useful values....

    DO i=1,land_pts
      ip = land_index(i)
      ix = NINT( (longitude(ip,1) - regLon1) / regDlon ) + 1
      iy = NINT( (latitude(ip,1) - regLat1) / regDlat ) + 1
      roffAccumGrid(ix,iy) = roffAccumLand(ip)
    ENDDO

  ELSE

!   The "main" model grid is a regular, lat/lon grid.
!   land_index can be used to scatter the points across this grid.
    DO ip=1,land_pts
      CALL getXYPos( land_index(ip),nx,ny,ix,iy )
      roffAccumGrid(ix,iy) = roffAccumLand(ip)
    ENDDO

  ENDIF

!-------------------------------------------------------------------------------
! Get runoff from land grid onto routing grid, and convert from flux density
! (kg m-2 s-1) to flux (kg s-1).
!-------------------------------------------------------------------------------
  runoff(:,:) = 0.0

  IF ( .NOT. routeRegrid ) THEN

!-------------------------------------------------------------------------------
!   "Main" and routing grids are identical (including points in same order).
!   No need to regrid, simply convert units.
!-------------------------------------------------------------------------------
    runoff(:,:) = roffAccumGrid(:,:) * routeArea(:,:)

  ELSE

!-------------------------------------------------------------------------------
!   Regrid from land to routing grid (and convert to flux).
!-------------------------------------------------------------------------------
!   Calculate an array size used in regridding.
    maxLin = ( nx + nxRoute )*( ny + nyRoute )

    CALL route_regrid( maxLin,roffAccumGrid,runoff )

!   Units of runoff are currently those of roffAccumGrid *(kg m-2 s-1).
!   Convert from flux density to flux.
    runoff(:,:) = runoff(:,:) * routeArea(:,:)

  ENDIF

!-------------------------------------------------------------------------------
! Call the routing routine.
!-------------------------------------------------------------------------------
  CALL route_trip( routeLength,runoff,outFlow )

!-------------------------------------------------------------------------------
! Add runoff on TRIP sea points to the river flow variable.
! Regridding may have resulted in some land runoff appearing in TRIP sea
! gridboxes, so we need to account for this water. Even without regridding,
! a land point that is sea in TRIP will not have been added to flow.
!-------------------------------------------------------------------------------
  WHERE( .NOT. routeMask(:,:) ) outFlow(:,:) = runoff(:,:)

!  do iy=1,nyRoute
!    do ix=1,nxRoute
!      if ( .NOT. routeMask(ix,iy) ) outFlow(ix,iy) = runoff(ix,iy)
!    enddo
!  enddo


!-------------------------------------------------------------------------------
! Save any extra diagnostics that are requested for offline application.
!-------------------------------------------------------------------------------
  IF ( useRflowDiag ) rFlowDiag(:,:) = outFlow(:,:)
  IF ( useRrunDiag ) rRunDiag(:,:) = runoff(:,:)

  END SUBROUTINE route_drive_trip

!###############################################################################
!###############################################################################
!###############################################################################
! subroutine getRouteLen
! Driver routine that calls procedures that calculates distance between gridpoints.

  SUBROUTINE getRouteLen( x,xNext,y,yNext,routeLat,routeLon,length )

  USE earth_utils, ONLY :  &
!  imported procedures
     giveLength

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     npRoute,nxRoute,nyRoute,routeDlat,routeDlon  &
!  imported arrays with intent(in)
    ,flowDirDelta,flowDirSet

  USE time_loc, ONLY :  &
!  imported arrays with intent(in)
     latitude,longitude
!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Array arguemnts with intent(in)
!-------------------------------------------------------------------------------

  INTEGER, INTENT(in) ::  &
    x(npRoute)      &!  x index on the routing grid of active routing points
   ,xNext(npRoute)  &!  x location of the next downstream point
!                         <0 indicates next point is off grid and gives the
!                         direction of the next point, in terms of index in
!                         flowDirSet.
   ,y(npRoute)      &!  y index on the routing grid of active routing points
   ,yNext(npRoute)   !  y location of the next downstream point
!                         <0 interpreted as for xNext.

  REAL, INTENT(in) ::  &
    routeLat(nxRoute,nyRoute)  &!  latitude of points on routing grid (degrees)
   ,routeLon(nxRoute,nyRoute)   !  longitude of points on routing grid (degrees)

!-------------------------------------------------------------------------------
! Array arguments with intent(out)
!-------------------------------------------------------------------------------

  REAL, INTENT(out) ::  &
    length(nxRoute,nyRoute)  !  distance between each gridpoint and downstream gridpoint (m)

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------

  INTEGER ::  &
    ip,ix,iy      &!  work
   ,jx,jy          !  work

  REAL ::  &
    dx,dy          !  work

!-------------------------------------------------------------------------------

! Only need to calculate length at land points on routing grid.

  DO ip=1,npRoute

!-------------------------------------------------------------------------------
!   Get coords of this point and any point immediately downstream.
!-------------------------------------------------------------------------------
    ix = x(ip)
    iy = y(ip)
    jx = xNext(ip)
    jy = yNext(ip)

    IF ( jx > 0 ) THEN

!-------------------------------------------------------------------------------
!     There is an outflow direction, and the downstream point is on the grid.
!-------------------------------------------------------------------------------
      length(ix,iy) = giveLength( routeLat(ix,iy),routeLon(ix,iy)  &
                               ,routeLat(jx,jy),routeLon(jx,jy) ) * 1000.0

    ELSEIF ( jx == -9 ) THEN

!-------------------------------------------------------------------------------
!       Outflow to sea, or no outflow (pit or depression).
!       This is relatively common, hence treatment separately from flow across grid edge.
!       Use distance to an adjacent point at same latitude.
!-------------------------------------------------------------------------------
        length(ix,iy) = giveLength( routeLat(ix,iy),routeLon(ix,iy)  &
                                 ,routeLat(ix,iy),routeLon(ix,iy)+routeDlon ) * 1000.0

    ELSE

!-------------------------------------------------------------------------------
!       Outflow across the edge of the grid. We would still like to calculate a length,
!       so that, as far as possible, a regional implementation gives the same results
!       as a global implementation, even for these edge points.
!       NB This works for outflow from grid, but if a regional run is missing an
!       inflow from outside the grid, answers will generally be different in global
!       and regional implementations.
!       The lat/lon of the downstream point
!       is not available in these cases, so we ASSUME a "regular" lat/lon grid.
!       jx is -1 * index in flowDirSet.
!-------------------------------------------------------------------------------
!       Get number of gridboxes in each direction to the downstream location.
        dx = REAL( flowDirDelta(ABS(jx),1) )
        dy = REAL( flowDirDelta(ABS(jx),2) )
        length(ix,iy) = giveLength( routeLat(ix,iy),routeLon(ix,iy)  &
                            ,routeLat(ix,iy)+dy*routeDlat  &
                            ,routeLon(ix,iy)+dx*routeDlon ) * 1000.0

    ENDIF

  ENDDO

  END SUBROUTINE getRouteLen
!###############################################################################
!###############################################################################
!###############################################################################
! subroutine route_trip
! Calculate river outflow and update channel storage for TRIP routing model.
! Based on UM routine outflow1.
! At present there is no regridding and therefore no need to account for water
! that may have been allocated to e.g. a sea point during regridding.

  SUBROUTINE route_trip( routeLength,runoff,outFlow )

  USE grid_utils, ONLY : &
!  imported procedures
     getXYpos

  USE prognostics, ONLY :  &
!  imported arrays with intent(inout)
     routeStore

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     npRoute,nxRoute,nyRoute,routeMeander,routeSpeed,routeTimeStep  &
!  imported arrays with intent(in)
    ,routeIndex,routeNext,routeOrder

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timeStep
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    i,ip,ix,iy,jx,jy       !  work/loop counters

  REAL, INTENT(in) ::  &!  in arrays
    routeLength(nxRoute,nyRoute)  &!  distance between gridpoints (m)
   ,runoff(nxRoute,nyRoute)        !  rate of runoff generation in each gridbox (kg s-1)

  REAL, INTENT(out) ::  &!  out arrays
    outFlow(nxRoute,nyRoute)   !  rate of channel flow leaving gridbox (kg s-1)

  REAL ::  &!  local scalars
    coeff     &!  coefficient in the routing model (s-1)
   ,dt        &!  timestep of routing model (s)
   ,storeOld   !  channel storage (kg)

  REAL ::  &!  local arrays
    inflow(nxRoute,nyRoute)  !  rate of channel flow entering gridbox (kg s-1)
!-------------------------------------------------------------------------------

  dt = REAL(routeTimeStep) * timeStep

! Initialise inflow with runoff generated over each gridbox.
  inflow(:,:) = runoff(:,:)

! Initialise outlow (for clarity at non-routing points).
  outFlow(:,:) = 0.0

! Loop over active routing points.
  DO i=1,npRoute

!   Get index (location in routing vector) of the point to consider.
    ip = routeOrder(i)

!   Get location of this point in routing grid.
    CALL getXyPos( routeIndex(ip),nxRoute,nyRoute,ix,iy )

!   Calculate the coefficient "c" of the model.
!   c=u/(d*r), where u is effective flow speed, d is distance between gridpoints,
!   and r is meander ratio.
    coeff = routeSpeed / ( routeLength(ix,iy) * routeMeander )

!   Save value of channel storage at start of timestep.
    storeOld = routeStore(ix,iy)

!   Calculate channel storage at end of timestep.
!   Eqn.4 of Oki et al, 1999, J.Met.Soc.Japan, 77, 235-255.
    routeStore(ix,iy) = storeOld * EXP(-(coeff*dt))  &
                   + (1.0 - EXP(-(coeff*dt))) * inflow(ix,iy) / coeff

!   Calculate outflow as inflow minus change in storage.
    outflow(ix,iy) = inflow(ix,iy) + (storeOld-routeStore(ix,iy)) / dt

!   Add outflow to inflow of next downstream point.
    IF ( routeNext(ip) > 0 ) THEN
!     Get location in grid of next downstream point.
      CALL getXYPos( routeNext(ip),nxRoute,nyRoute,jx,jy )
      inflow(jx,jy) = inflow(jx,jy) + outflow(ix,iy)
    ENDIF

  ENDDO   !  points

  END SUBROUTINE route_trip

!###############################################################################
!###############################################################################
!###############################################################################
! subroutine route_regrid
! Handles regridding of runoff from "main" to routing grids.

  SUBROUTINE route_regrid( maxLin,roffAccumGrid,runoff )
!-------------------------------------------------------------------------------

! Description:
!   Regrids runoff from a source grid to a target grid (the routing grid).
!   Both grids must be regular in latitude and longitude.
!
!-------------------------------------------------------------------------------
! Modules used:

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     nxRoute,nyRoute,routeDlat,routeDlon,routeLat1,routeLon1

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     nxGrid,nyGrid,regDlat,regDlon,regLat1,regLon1

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar arguments with intent(in)

  INTEGER, INTENT(in) :: maxLin  !  a vector length

! Array arguments with intent(in)

  REAL, intent(in) :: roffAccumGrid(nxGrid,nyGrid) !  runoff rate on model grid (kg m-2 s-1)

! Array arguments with intent(out)

  REAL, intent(out) :: runoff(nxRoute,nyRoute)  !  runoff rate on routing grid (kg m-2 s-1)

!-------------------------------------------------------------------------------
! Local scalar variables.

  INTEGER :: adjust  !  mode (option number) of area averaging
  INTEGER :: icode   !  exit code from subroutines
  INTEGER :: ix      !  loop counter
  INTEGER :: iy      !  loop counter
  INTEGER :: maxL    !  number of entries in output from pre_areaver

  LOGICAL :: cyclic_srce  !  TRUE source (model) grid is cyclic in x
  LOGICAL :: cyclic_targ  !  TRUE target (model) grid is cyclic in x
  LOGICAL :: invert_srce  !  TRUE source (model) grid runs N to S
!                            FALSE source grid runs S to N
  LOGICAL :: spherical  !  TRUE  coordinates are lat/lon on a sphere
!                          FALSE Cartesian axes.
  LOGICAL :: want       !  Value of masks at locations for which data are required

  CHARACTER(len=80) :: cmessage  !  error message from subroutines

!-------------------------------------------------------------------------------
! Local array variables.

  INTEGER :: count_targ(nxRoute,nyRoute)  !  number of model gridboxes that
!                                            contribute to each routing gridbox
  INTEGER :: base_targ(nxRoute,nyRoute)   !  base (starting) index for each
!                   routing gridbox (i.e. location of first element in lists)
  INTEGER :: index_srce(maxLin)          !  list of source (model) gridboxes
!                   that contribute to each routing gridbox

  REAL :: sourceLat(nyGrid+1)  !  latitudes of edges of source (model) gridboxes.
!         First value is the S edge of first gridbox, all other values are the
!         N edge of each gridbox.
  REAL :: sourceLon(nxGrid+1)  !  longitudes of edges of source (model) gridboxes
!         First value is the W edge of first gridbox, all other values are the
!         E edge of each gridbox.
  REAL :: targetLat(nyRoute+1)  !  latitudes of edges of target (mrouting) gridboxes.
!         First value is the S edge of first gridbox, all other values are the
!         N edge of each gridbox.
  REAL :: targetLon(nxRoute+1)  !  longitudes of edges of target (routing) gridboxes
!         First value is the W edge of first gridbox, all other values are the
!         E edge of each gridbox.

  REAL :: adjust_targ(nxRoute,nyRoute)  !  adjustment factors (not used)
  REAL :: weight(maxLin)     !   lists of weights for each source (model) gridbox

  LOGICAL :: mask_srce(nxGrid,nyGrid)
  LOGICAL :: mask_targ(nxRoute,nyRoute)
!-------------------------------------------------------------------------------

! Set values.
  spherical = .TRUE.    !  Calculations are for lat/lon coordinates.
  adjust = 0            !  "normal" area averaging
  maxL = maxlIn
  invert_srce = .FALSE. !  model grid runs S to N

! Set coordinates of edges of model gridboxes.
  DO ix=1,nxGrid+1
    sourceLon(ix) = regLon1 + (REAL(ix)-0.5) * regDlon
  ENDDO
  DO iy=1,nyGrid+1
    sourceLat(iy) = regLat1 + (REAL(iy)-0.5) * regDlat
  ENDDO

! Set coordinates of edges of routing gridboxes.
  DO ix=1,nxRoute+1
    targetLon(ix) = routeLon1 + (REAL(ix)-0.5) * routeDlon
  ENDDO
  DO iy=1,nyRoute+1
    targetLat(iy) = routeLat1 + (REAL(iy)-0.5) * routeDlat
  ENDDO

! Decide if grids are cyclic in longitude.
  cyclic_srce = .FALSE.
  cyclic_targ = .FALSE.
  IF ( REAL(nxGrid)*regDlon > 359.9 ) cyclic_srce = .TRUE.
  IF ( REAL(nxRoute)*routeDlon > 359.9 ) cyclic_targ = .TRUE.

! Set masks to indicate that all points in both grids are to be used.
  want = .TRUE.
  mask_srce(:,:) = want
  mask_targ(:,:) = want

! Call setup routing for averaging.
  CALL PRE_AREAVER( nxGrid,sourceLon,nyGrid,sourceLat  &
                   ,cyclic_srce,nxGrid,want,mask_srce  &
                   ,nxRoute,targetLon,nyRoute,targetLat  &
                   ,cyclic_targ,spherical  &
                   ,maxL,count_targ,base_targ,index_srce,weight,icode,cmessage )

! Call averaging routine.
  CALL DO_AREAVER( nxGrid,nyGrid,nxGrid    &
            ,invert_srce,roffAccumGrid,nxRoute,nyRoute,count_targ &
            ,base_targ,nxRoute,want,mask_targ,index_srce,weight,adjust &
            ,runoff,adjust_targ,icode,cmessage )

  END SUBROUTINE route_regrid
!###############################################################################
!###############################################################################
!###############################################################################
#endif
