#if defined(OASIS4)
MODULE OASIS4_grid32_mod
! Description: This module contains interface routines for OASIS4
!              PRISM library calls and associated variables. 
!              The module must be compiled without the UM's normal 
!              promotion of all real and integer variables to 64 bits.
!
!              PRISM calls refer to the prebuilt OASIS4 library
!              and associated include files. The PRISM modules 
!              and library must be specified in the compile and
!              load options. FCM must be prevented from 
!              performing dependency checks on PRISM.!              
!
!              Link stage will also require linking to an existing
!              XML library and potentially netcdf library depending
!              on the options used to build the OASIS4 PRISM libs. 
!              (All done, typically via a compile override file).
!
! History:
! UM Version    Date    Description
! ----------  --------  ----------------------------------
!    6.4      Jan 2007  Original Code. R. Barnes, R. Hill
!                                      A. Treshansky
!=================================================================

      USE PRISM
      USE OASIS4_atm_data_mod

      IMPLICIT NONE

#include "c_kinds.h"

      Type(PRISM_Time_Struct) :: date
      Type(PRISM_Time_Struct) :: window(2)

      ! Allocate space for several different grid types
      ! Arbitrarily 6, but hope to only need 3 eventually.
      INTEGER (kind=integer32), PARAMETER :: n_grids=6

      INTEGER (kind=integer32)  :: grid_id(n_grids)

      INTEGER (kind=integer32), PARAMETER :: t_grid = 1
      INTEGER (kind=integer32), PARAMETER :: u_grid = 2
      INTEGER (kind=integer32), PARAMETER :: v_grid = 3

      ! These values (t_grid, u_grid, v_grid) must match the values
      ! in fldtype.h (fld_type_p, fld_type_u, fld_type_v) or there will
      ! be problems; the fldtype values are used in c_grid_64.F
      INTEGER (kind=integer32), PARAMETER :: n_fields = 20      

      INTEGER (kind=integer32)  :: n_invar
      INTEGER (kind=integer32)  :: n_outvar
      INTEGER (kind=integer32)  :: um_start_date(6)

      INTEGER (kind=integer32), PARAMETER :: n_tranvar=40 ! max no. of 
                                                          ! transients

      INTEGER (kind=integer32)  :: var_id(n_tranvar)  ! Incoming/outgoing
      LOGICAL :: no_neighbour_N

CONTAINS
      subroutine GRID32_OASIS4(comp_64,MODEL_BASIS_TIME,ierr64    &
     &, valid_sh64,actual_sh64,nGridDims64,grid_type64,nCorners64 &
     &, t_longitudes,t_latitudes,t_verticals,v_off64              &
     &, nbr_subdomains64, offset_array64, extent_array64)

! Description: This routine calls the OASIS4 PRISM library routines
!              to define grid names, IDs, grid-point locations,
!              corner grid-box locations, masks (if any) and 
!              transient variables.
!
! History:
! UM Version    Date    Description
! ----------  --------  ----------------------------------
!    6.4      Jan 2007  Original Code. R. Hill, R Barnes,
!                                              A. Treshansky
!=================================================================

      implicit none

      integer (kind=integer64) :: comp_64
      integer (kind=integer32) :: comp_id
      character*32 :: grid_name
      integer (kind=integer64) :: ierr64
      integer (kind=integer32) :: ierror
      integer (kind=integer32) :: info

      real (kind=real64) :: v_off64

      integer (kind=integer64) :: MODEL_BASIS_TIME(6)

      logical, parameter :: verbose = .true.

  ! this module defines an Arakawa C grid:
  !   - this is an instance of PRISM_reglonlatvrt
  !   - 3D grid
  !   - points are created for the first time here
  !   - each grid cell is a cube w/ 8 vertices
  !   - C grid has t, u, v points

      integer (kind=integer64) :: grid_type64
      integer (kind=integer32) :: grid_type
      integer (kind=integer64) :: nGridDims64
      integer (kind=integer32) :: nGridDims
      logical, parameter :: new_points = .true.
      integer (kind=integer64) :: nCorners64
      integer (kind=integer32) :: nCorners
      integer (kind=integer64) :: t_points_64,u_points_64,v_points_64
      integer (kind=integer32) :: t_points_id,u_points_id,v_points_id

      integer (kind=integer64) :: valid_sh64(2,nGridDims64)
      integer (kind=integer32) :: valid_shape(2,nGridDims64)
      integer (kind=integer64) :: actual_sh64(2,nGridDims64)
      integer (kind=integer32) :: actual_shape(2,nGridDims64)

      ! V grid has 1 fewer row
      integer (kind=integer32) :: valid_shape_v(2,nGridDims64)
      integer (kind=integer32) :: actual_shape_v(2,nGridDims64)  


      integer (kind=integer64) :: nbr_subdomains64
      integer (kind=integer64) ::                                       &
     & offset_array64(n_fields,nbr_subdomains64,nGridDims64)
      integer (kind=integer64) ::                                       &
     & extent_array64(n_fields,nbr_subdomains64,nGridDims64)
      integer(kind=integer64) :: out_unit_64


      integer (kind=integer32) :: var_nodims(2)
      integer (kind=integer32) :: blank_mask_id
      integer (kind=integer32) :: mask_id(3)  ! One each for T U V

      ! local variables...
      integer (kind=integer32) :: i,j,k,ip1

      ! Incoming longitudes
      real (kind=real64), dimension(actual_sh64(1,1):actual_sh64(2,1))  &
     & :: t_longitudes, u_longitudes, v_longitudes

      ! Incoming latitudes 
      real (kind=real64), dimension(actual_sh64(1,2):actual_sh64(2,2))  &
     & :: t_latitudes, u_latitudes
     
      ! Incoming vertical positions 
      real (kind=real64), dimension(actual_sh64(1,3):actual_sh64(2,3))  &
     & :: t_verticals, u_verticals, v_verticals

      ! locally defined corners on c_grid 
      real(kind=real64),dimension(actual_sh64(1,1):actual_sh64(2,1),2)  &
     & :: corner_longitudes, corner_longitudes_u

      ! locally defined corners on c_grid
      real(kind=real64),dimension(actual_sh64(1,2):actual_sh64(2,2),2)  &
     & :: corner_latitudes

      ! Latitudes for v pooints
      real(kind=real64), DIMENSION(:), ALLOCATABLE   &
     & :: v_latitudes

      ! locally calculated lat corners v grid.
      real(kind=real64),dimension(:,:),ALLOCATABLE  &
     & :: corner_latitudes_v

     ! local vertical corners suitable for p, u and v grids
      real(kind=real64),dimension(actual_sh64(1,3):actual_sh64(2,3),2)  &
     & :: corner_verticals

      real(kind=real64) :: u_off  ! U grid EW staggering
      real(kind=real64) :: v_off  ! V grid NS staggering

      ! local (32-bit) versions of arguments to PRISM_def_partition
      integer(kind=integer32)  :: nbr_subdomains
      integer(kind=integer32), dimension(:,:,:), allocatable ::         &
     &                            offset_array
      integer(kind=integer32), dimension(:,:,:), allocatable ::         &
     &                            extent_array

!
! define the grid and its relatives
!
      INTEGER :: iat42
      INTEGER :: jat42
      INTEGER :: kat42

      LOGICAL,  DIMENSION(:,:,:), ALLOCATABLE :: t_mask
      LOGICAL,  DIMENSION(:,:,:), ALLOCATABLE :: u_mask
      LOGICAL,  DIMENSION(:,:,:), ALLOCATABLE :: v_mask

      INTEGER :: PRISM_method_id, PRISM_mask_id
      INTEGER :: PRISM_points_t
      INTEGER :: PRISM_points_u
      INTEGER :: PRISM_points_ur
      INTEGER :: PRISM_points_us
      INTEGER :: PRISM_points_v


      LOGICAL :: ll_points, ll_mask
      INTEGER :: II,JJ

      ! Set up our atmosphere grid dimensions
      iat42 =  actual_sh64(2,1)
      jat42 =  actual_sh64(2,2)
      kat42 =  1

      ! Move valid and actual shape arrays to local copies.
      valid_shape(:,:) = valid_sh64(:,:) 
      actual_shape(:,:) = actual_sh64(:,:)

      ! Set valid shape for v grid
      valid_shape_v(:,:) = valid_shape(:,:)
      actual_shape_v(:,:) = actual_shape(:,:)

      ! If this PE is at the northern domain limit
      ! subtract 1 row from the v grid dimension.
      ! What this effectively means is that we're saying
      ! if a PE owns a T point, then it must also
      ! own the v point to the north of the t point 
      ! (and also the u point to the east).
      IF (no_neighbour_N) THEN
         valid_shape_v(2,2) = valid_shape_v(2,2)-1
         actual_shape_v(2,2) = actual_shape_v(2,2)-1
      end if 

! 2. Define the gaussian grid and its relatives

! 2.1.2. The  lon, lat, z and corner values are passed into this routine

! 2.1.3. Allocate space for the masks
      allocate(t_mask(iat42,jat42,kat42))
      allocate(u_mask(iat42,jat42,kat42))
      allocate(v_mask(iat42,valid_shape_v(2,2),kat42))


      ! Set up masks on each set of grid points
      DO J = valid_shape(1,2),valid_shape(2,2)
         DO I = valid_shape(1,1),valid_shape(2,1)
            t_mask(i,j,1) = (OASIS4_TMASK(I,J,1))
            u_mask(i,j,1) = (OASIS4_UMASK(I,J,1))
         end do
      end do

      DO J = valid_shape_v(1,2),valid_shape_v(2,2)
         DO I = valid_shape_v(1,1),valid_shape_v(2,1)
            v_mask(i,j,1) = (OASIS4_VMASK(I,J,1))
         end do
      end do



! 2.2. Define the grid and its relatives with the PRISM interface

      ! Move component ID to local copy
      comp_id = comp_64
      OASIS4_comp_id = comp_id


      ! 1st we do the "main" grid
      grid_name = "c_grid"
      call PRISM_def_grid (grid_id(t_grid),grid_name,comp_id,    &
     &                      valid_shape,PRISM_reglonlatvrt,ierror)

      ! The U and T grids are effectively one and the same.

      ! A bug in OASIS4 prevents us from receieving things 
      ! ONTO U points - define a separate grid for U points
      grid_name = "c_gridu"
      call PRISM_def_grid (grid_id(u_grid),grid_name,comp_id,    &
     &                      valid_shape,PRISM_reglonlatvrt,ierror)
       
      ! The V and T grids are different.

      ! Define the V grid. This is the same shape and size
      ! as the T grid in the EW dimension but 1 row fewer in the NS
      ! dimension.
      grid_name = "c_gridv"
      call PRISM_def_grid (grid_id(v_grid),grid_name,comp_id,    &
     &                      valid_shape_v,PRISM_reglonlatvrt,ierror)


! 2.2.2. Define the positions of the points on each grid

      ! Calculate the U grid EW stagger
      u_off = (t_longitudes(2) - t_longitudes(1))/2.0

      ! Set up the V grid NS stagger - which we already have
      ! available from elsewhere.       
      v_off = v_off64

      ! Define the T points
      ll_points = .true.  ! New points being defined
      call PRISM_set_points(PRISM_points_t, "AT_points",    &
     & grid_id(t_grid),                                      &
     & valid_shape,                                          &
     & t_longitudes(valid_shape(1,1):valid_shape(2,1)),      & 
     & t_latitudes(valid_shape(1,2):valid_shape(2,2)),       &
     & t_verticals(valid_shape(1,3):valid_shape(2,3)),       &
     & ll_points, ierror )

      ! Set up longitudes for U and V points    
      DO I = valid_shape(1,1),valid_shape(2,1)
         u_longitudes(i) = t_longitudes(i) + u_off 
         v_longitudes(i) = t_longitudes(i) 
      end do
        
      ! U grid latitudes are the same as T grid latitudes
      DO J = valid_shape(1,2),valid_shape(2,2) 
         u_latitudes(J) = t_latitudes(J)
      end do

      ! For the V points, we have 1 row fewer. Allocate space
      ! for latitudes and latitude corners      
      allocate(v_latitudes(valid_shape_v(2,2)))
      allocate(corner_latitudes_v(valid_shape_v(2,2),2))
 
      DO J = valid_shape_v(1,2),valid_shape_v(2,2)
         ! We restrict our corner latitudes to +/- 90 degrees
         ! may need to develop this to cope with LAMS? Maybe
         ! that doesnt matter in LAMS. Nobody ever couples LAMs
         ! but "coupling" has a slightly different meaning
         ! in a FLUME context.
         v_latitudes(J) = t_latitudes(J) + v_off

         corner_latitudes_v(j,1) = v_latitudes(j) - v_off
         corner_latitudes_v(j,2) = v_latitudes(j) + v_off
      end do

      ! Set verticals for U and V points
      DO K = valid_shape(1,3),valid_shape(2,3)
         u_verticals(K) = t_verticals(K)
         v_verticals(K) = t_verticals(K)
      end do
        
      ! Define the U points
      ll_points = .true.
      ! Our U grid is the same size as the T grid 
      call PRISM_set_points(PRISM_points_us, "AU_points",   &
     & grid_id(t_grid),                                      &
     & valid_shape,                                          &
     & u_longitudes(valid_shape(1,1):valid_shape(2,1)),      &
     & u_latitudes(valid_shape(1,2):valid_shape(2,2)),       &
     & u_verticals(valid_shape(1,3):valid_shape(2,3)),       &
     & ll_points, ierror )


      ! Define the UR points
      ll_points = .true.
      ! Our UR grid is the same size as the T grid
      call PRISM_set_points(PRISM_points_ur, "AUR_points",  &
     & grid_id(u_grid),                                      &
     & valid_shape,                                          &
     & u_longitudes(valid_shape(1,1):valid_shape(2,1)),      &
     & u_latitudes(valid_shape(1,2):valid_shape(2,2)),       &
     & u_verticals(valid_shape(1,3):valid_shape(2,3)),       &
     & ll_points, ierror )


      ! Our V grid is 1 less row than the T grid.
      ll_points = .true.
      call PRISM_set_points(PRISM_points_v, "AV_points",    &
     & grid_id(v_grid),                                      &
     & valid_shape_v,                                        &
     & v_longitudes(valid_shape_v(1,1):valid_shape_v(2,1)),  &
     & v_latitudes(valid_shape_v(1,2):valid_shape_v(2,2)),   &
     & v_verticals(valid_shape_v(1,3):valid_shape_v(2,3)),   &
     & ll_points, ierror )


! 2.2.3. Define the land-sea mask on t-grid
      ll_mask = .true.
      call PRISM_set_mask(mask_id(t_grid),grid_id(t_grid), &
     &                     valid_shape,t_mask(:,:,:),       &
     &                     ll_mask, ierror )

! Mask on u points
      ll_mask = .true.
      call PRISM_set_mask(mask_id(u_grid),grid_id(u_grid), &
     &                     valid_shape,u_mask(:,:,:),       &
     &                     ll_mask, ierror )

! Mask on v points
      ll_mask = .true.
      call PRISM_set_mask(mask_id(v_grid),grid_id(v_grid), &
     &                     valid_shape_v,v_mask(:,:,:),     &
     &                     ll_mask, ierror )
    

! 2.2.4. Define the corners

      ! The RH longitude corners are effectively at the same 
      ! position as the u points, though they needn't be.
      do i = valid_shape(1,1), valid_shape(2,1)
         corner_longitudes(i,1) = t_longitudes(i)-u_off
         corner_longitudes(i,2) = t_longitudes(i)+u_off
 
         corner_longitudes_u(i,1) =corner_longitudes(i,1)+u_off
         corner_longitudes_u(i,2) =corner_longitudes(i,2)+u_off
      end do

      ! Set up our corner latitudes - restricting things 
      ! to the real world!
      do j = valid_shape(1,2), valid_shape(2,2)
         corner_latitudes(j,1) = MAX((t_latitudes(j)-v_off),-90.0)
         corner_latitudes(j,2) = MIN((t_latitudes(j)+v_off),+90.0)
      end do

      ! Set the vertical corners and grid positions
      ! it must touch the ground to couple!
      do k=valid_shape(1,3), valid_shape(2,3)
         t_verticals(k) = 0.0  !
         ! with a delta-z of 1, working out the corners is trivial...
         corner_verticals(k,1) = (k-1)*1.0
         corner_verticals(k,2) = (k)*1.0
      end do


      ! Set the corners for points on the c_grid (T and U points)
      call PRISM_set_corners( grid_id(t_grid), 8, valid_shape,      &
     &  corner_longitudes, corner_latitudes, corner_verticals,       &
     &  ierror )
 
      ! Set the corners for points on the c_gridu (UR points)
      call PRISM_set_corners( grid_id(u_grid), 8, valid_shape,      &
     &  corner_longitudes_u, corner_latitudes, corner_verticals,       &
     &  ierror )

      ! Define the V grid corners. Only the lats
      ! differ from the T grid (and we have 1 less row)
      call PRISM_set_corners( grid_id(v_grid), 8, valid_shape_v,    &
     &  corner_longitudes, corner_latitudes_v, corner_verticals,    &
     &  ierror )


      ! Define the MPP partition for OASIS4. i.e. this is the 
      ! MPP decompostion as the atmosphere currently knows
      ! it, defined in a way that OASIS4 can understand it.
      nbr_subdomains = nbr_subdomains64

      allocate(offset_array(size(offset_array64,1),                     &
     & size(offset_array64,2),size(offset_array64,3)))

      allocate(extent_array(size(extent_array64,1),                     &
     & size(extent_array64,2),size(extent_array64,3)))

      offset_array = offset_array64
      extent_array = extent_array64

      ! Set up the MPP partition on the t grid
      call PRISM_def_partition(grid_id(t_grid), nbr_subdomains   &
     &, offset_array(t_grid,:,:), extent_array(t_grid,:,:), ierror )
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_def_partition (t_grid)",ierror
      else
        write(6,*) "PRISM_def_partition (t_grid) called correctly"
      end if

      ! Set up the MPP partition on the u grid
      call PRISM_def_partition(grid_id(u_grid), nbr_subdomains   &
     &, offset_array(u_grid,:,:), extent_array(u_grid,:,:), ierror )
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_def_partition (u_grid)",ierror
      end if

      ! Call def partition for the v points
      call PRISM_def_partition(grid_id(v_grid), nbr_subdomains   &
     &, offset_array(v_grid,:,:), extent_array(v_grid,:,:), ierror )
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_def_partition (v_grid)",ierror
      end if


      var_nodims(1) = nGridDims64
      var_nodims(2) = 0

      ! Set up a blank mask
      blank_mask_id = PRISM_UNDEFINED
      !var_type = PRISM_double_precision

      ! Define output ustress - u points on a c_grid
      call PRISM_def_var ( var_id(vind_taux),                    &
     &  "taux", grid_id(t_grid),                                 &
     &  PRISM_points_us, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )

      ! Define output vstress - v points on a c_gridv
      call PRISM_def_var ( var_id(vind_tauy),                    &
     &  "tauy", grid_id(v_grid),                                 &
     &  PRISM_points_v, blank_mask_id, var_nodims, valid_shape_v, &
     &  PRISM_Double_Precision, ierror )

      ! Define output wme
      call PRISM_def_var ( var_id(vind_wme),                     &
     &  "wme", grid_id(t_grid),                                  &
     &  PRISM_points_t, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )

      ! Define output pme
      call PRISM_def_var ( var_id(vind_pme),                     &
     &  "pme", grid_id(t_grid),                                  &
     &  PRISM_points_t, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )

      ! Define output runoff
      call PRISM_def_var ( var_id(vind_runoff),                  &
     &  "runoff", grid_id(t_grid),                               &
     &  PRISM_points_t, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )

      ! Define output pen solar radn (blue band) 
      call PRISM_def_var ( var_id(vind_pen_solar),               &
     &  "pen_solar", grid_id(t_grid),                            &
     &  PRISM_points_t, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )


      ! Define output field for total heat flux.
       call PRISM_def_var ( var_id(vind_tot_atm_hflux),          &
     &  "heat_flux", grid_id(t_grid),                            &
     &  PRISM_points_t, blank_mask_id, var_nodims, valid_shape,  &
     &  PRISM_Double_Precision, ierror )

      ! Define incoming fields
      ! Define input field for ocean sst.
       call PRISM_def_var ( var_id(vind_ocn_sst),                &
     &  "ocn_sst", grid_id(t_grid),                              &
     &  PRISM_points_t, mask_id(t_grid), var_nodims, valid_shape,&
     &  PRISM_Double_Precision, ierror )

      ! Define input field for ocean freeze (AICE as far as the UM
      ! is concerned.
       call PRISM_def_var ( var_id(vind_ocn_freeze),             &
     &  "ocn_freeze", grid_id(t_grid),                           &
     &  PRISM_points_t, mask_id(t_grid), var_nodims, valid_shape,&
     &  PRISM_Double_Precision, ierror )

      ! Define input field for ocean u. u points on c_grid
       call PRISM_def_var ( var_id(vind_ocn_u),                   &
     &  "sunocean", grid_id(u_grid),                              &
     &  PRISM_points_ur, mask_id(u_grid), var_nodims, valid_shape,&
     &  PRISM_Double_Precision, ierror )

      ! Define input field for ocean v. v points on c_gridv
       call PRISM_def_var ( var_id(vind_ocn_v),                   &
     &  "svnocean", grid_id(v_grid),                              &
     &  PRISM_points_v,mask_id(v_grid), var_nodims, valid_shape_v,&
     &  PRISM_Double_Precision, ierror )


      ! Finish the PRISM definition phase and perform inter component
      ! integrity checking. 99 times out of 100 this is where 
      ! a failure will occur but don't expect any helpful
      ! messages!      
      call PRISM_enddef(ierror)
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_enddef"
      end if
      ierr64 = ierror

      ! Get hold of the job start date/time as defined in the SCC XML file
      date = PRISM_jobstart_date

      WRITE(6,*) "UM MODEL BASIS",MODEL_BASIS_TIME 

      ! Cross check the UM start date against the PRISM start date
      ! While we can't cross check with other component models directly,
      ! if all our component models cross check their own start date with
      ! the PRISM start date, we effectively have a system-wide cross 
      ! check. We probably ought to have a run length cross check too.    
      if (date%year.NE.MODEL_BASIS_TIME(1).or.                    &
     &    date%month.NE.MODEL_BASIS_TIME(2).or.                   &
     &    date%day.NE.MODEL_BASIS_TIME(3).or.                     &
     &    date%hour.NE.MODEL_BASIS_TIME(4).or.                    &
     &    date%minute.NE.MODEL_BASIS_TIME(5).or.                  &
     &    date%second.NE.MODEL_BASIS_TIME(6)) THEN 
          WRITE(6,*) "*********************************************"
          WRITE(6,*) "********    W A R N I N G   *****************"    
          WRITE(6,*) "UM Start date does not match PRISM Start date"
          WRITE(6,*) "While not strictly necessary for things to"
          WRITE(6,*) "work, it may create confusion if you go ahead"
          WRITE(6,*) " "
          WRITE(6,*) "The best way to fix this is to edit your" 
          WRITE(6,*) "SCC XML file - but remember, all your component"
          WRITE(6,*) "models must have the same start date so you may"
          WRITE(6,*) "need to adjust your UM job!" 
          WRITE(6,*) "        UM  PRISM"
          WRITE(6,*) "Year",MODEL_BASIS_TIME(1),date%year
          WRITE(6,*) "Month",MODEL_BASIS_TIME(2),date%month
          WRITE(6,*) "Day",MODEL_BASIS_TIME(3),date%day
          WRITE(6,*) "Hour",MODEL_BASIS_TIME(4),date%hour
          WRITE(6,*) "Minute",MODEL_BASIS_TIME(5),date%minute
          WRITE(6,*) "Second",MODEL_BASIS_TIME(6),date%second
      end if

      PRISM_timestep = -99.0 ! set invalid start timestep to act as flag

      write(6,*)' Atmos PRISM start Date is ',date

      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_put (heat_flux)",ierror,info
      end if
      ierr64 = ierror

      end subroutine grid32_OASIS4

!======================================================================
      subroutine set_oasis_window32(ierr64)
! Description: Given a date, this routine calculates the appropriate
!              time window either side of the new date
!              for use in put/get operations.
!======================================================================

      implicit none
      integer (kind=integer64) :: ierr64
      integer (kind=integer32) :: ierror

      window(1) = date
      window(2) = date
      ! So our window is +/- half our timestep

      call PRISM_calc_newdate(window(1),(prism_timestep/(-2.0)),ierror)
      call PRISM_calc_newdate(window(2),(prism_timestep/(2.0)),ierror)

      ierr64=ierror

      end subroutine set_oasis_window32

!======================================================================
      subroutine OASIS4_advance_date32(ierr64)
!
! Description: This routine advances the PRISM date and calculates
!              the appropriate time window either side of the new date
!              for use in put/get operations.
!
! History:
! UM Version    Date    Description
! ----------  --------  ----------------------------------
!    6.4      Jan 2007  Original Code. R. Hill
!======================================================================

      implicit none
      integer (kind=integer64) :: ierr64
      integer (kind=integer32) :: ierror

      ! Advance the date/time as used by PRISM
       
      call PRISM_calc_newdate ( date, PRISM_timestep, ierror )

      ierr64=ierror

      call set_oasis_window32(ierr64)

      end subroutine OASIS4_advance_date32

!=======================================================================
      subroutine PUT32_OASIS4(data_64,rowl_64,rows_64,ierr64)
!
! Description: PUT data to OASIS4 in 64 bit form (PRISM double
!              precision).
!
! History:
! UM Version        Date        Description
! ----------     -----------    ----------------------------------
!    6.4          Jan 2007      Original Code. R. Hill
!=======================================================================
      implicit none

      integer (kind=integer64) :: rowl_64
      integer (kind=integer64) :: rows_64

      integer (kind=integer64) :: ierr64
      integer (kind=integer32) :: ierror
      integer (kind=integer32) :: info
      logical, parameter :: verbose = .true.

      real (kind=real64),dimension(1:rowl_64,1:rows_64,1) :: data_64

      ! Hard coded for heat flux test
!      write(6,*)'atmos put window ',window

!      write(6,*) "VAR_IND at PUT is",var_ind,var_id(var_ind)


!      call PRISM_put_inquire(var_id(var_ind), date, window, &
!     &                 info, ierror )

!      write(6,*) "WILL WE PUT data_64?",info, ierror
             
!      write(6,*) "Data size for put is" , rowl_64,rows_64


      info =0
      ierror=0


      ! PUT the appropriate bit of the array
      call PRISM_put(var_id(var_ind),date,window,data_64,info,ierror)
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_put(data_64)",ierror,info,    &
     & var_id,date,window,data_64(rowl_64/2,rows_64/2,1)
      end if
      ierr64 = ierror

      end subroutine PUT32_OASIS4
!=======================================================================
      subroutine GET32_OASIS4(data_64,rowl_64,rows_64,o4info64,ierr64)
    
!
! Description: Get date from PRISM in 64 bit form (double
!              precision as far as PRISM is concerned, but 
!              standard precision for the UM).
! History:
! UM Version        Date        Description
! ----------     -----------    ----------------------------------
!    6.4          Jan 2007      Original Code. R. Hill
!=======================================================================
      implicit none

      integer (kind=integer64) :: rowl_64
      integer (kind=integer64) :: rows_64

      integer (kind=integer64) :: ierr64
      integer (kind=integer32) :: ierror
      integer (kind=integer64) :: o4info64
      integer (kind=integer32) :: info

      logical, parameter :: verbose = .true.
      real (kind=real64),dimension(1:rowl_64,1:rows_64,1) :: data_64

!      write(6,*)'atmos put window ',window

      info=0
      ierror=0

      call PRISM_get(var_id(var_ind),date,window,data_64,info,ierror)
      if (verbose .and. (ierror /= 0)) then
        write(6,*) "error in PRISM_get(data_64)",ierror,info,    &
     & var_id,date,window,data_64(rowl_64/2,rows_64/2,1)
      end if

      o4info64 = info 
      ierr64 = ierror

      end subroutine get32_OASIS4

end MODULE OASIS4_grid32_mod
#endif
