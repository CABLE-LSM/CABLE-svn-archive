#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      subroutine OASIS4_GRID64(comp_id,MODEL_BASIS_TIME,ierror)

! Description: This subroutine contains the intermediate call 
!              between the main UM code and the OASIS4 interface
!              code needed in order to hide the 32 bit dependencies of 
!              PRISM from the 64 bit dependencies of the UM.
!
! This routine is compiled -ew 64bit and computes the  basic values
! needed for OASIS4 grid definition.  However the OASIS4 software
! needs 32bit integers, so these calls are wrapped inside 
! subroutine grid32_OASIS4, compiled -dw with explicit kinds,
! where integers are converted from 64-32bit.
!
! History:
! UM Version        Date        Description
! ----------     -----------    ----------------------------------
!    6.4          Jan 2007      Original Code. R. Hill, R Barnes
!                                              A. Treshansky
!=================================================================

! DEPENDS ON: OASIS4_grid32_mod
      USE OASIS4_grid32_mod
! DEPENDS ON: OASIS4_atm_data_mod
      USE OASIS4_atm_data_mod

      implicit none

      integer         :: comp_id
      integer         :: ierror
      integer         :: MODEL_BASIS_TIME(6)

#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "version.h"
#include "model.h"


      logical, parameter :: verbose = .true.

! this module defines an Arakawa C grid:
!   - this is an instance of PRISM_reglonlatvrt
!   - 3D grid
!   - points are created for the first time here
!   - each grid cell is a cube w/ 8 vertices
!   - C grid has t, u, v points

      integer            :: grid_type
      integer, parameter :: nGridDims = 3 
      logical, parameter :: new_points = .true.
      integer, parameter :: nCorners = 8 
      integer :: t_points_id, u_points_id, v_points_id

      integer, save :: valid_shape(2,nGridDims)
      integer, save :: actual_shape(2,nGridDims)

      ! arguments needed for PRISM_def_partition, passed to grid32_OASIS4
      integer :: nbr_subdomains
      integer, dimension(:,:,:), allocatable :: offset_array
      integer, dimension(:,:,:), allocatable :: extent_array

      ! local variables...
      integer :: i,j,k

      real :: bottom_lat
      real :: v_off64
      real, dimension(:), allocatable :: t_longitudes
      real, dimension(:), allocatable :: t_latitudes
      real, dimension(:), allocatable :: t_verticals

      ! valid_shape = actual_shape since I'm not using halos, right?
      valid_shape(1,1) = 1
      valid_shape(1,2) = 1
      valid_shape(1,3) = 1

      ! Set the I extent on this PE (no halos for valid shape)
      valid_shape(2,1) = blsize(1,fld_type_p)


      ! The following is a short term modification to get MPP runs 
      ! to bit compare on up to 5 PEs in HadGAM1 resolution models.

      ! Set the J extent on this PE (no halos)
      ! Give 2 rows to each PE except the control PE
      IF (O4PE.EQ.O4CNTLPE) THEN
         ! Give most of the rows to the top pe 
         valid_shape(2,2) = glsize(2,fld_type_p)-(O4CNTLPE*2)
      ELSE
         valid_shape(2,2) =2 
      END IF
 
      ! Set the K extent on this PE (we dont decompose over K
      ! and, anyway, we only care about the bottom level)
      valid_shape(2,3) = 1 

      write(6,*)'grid64 - valid_shape ',valid_shape

      ! Now set up the actual shape. These can contain halos
      ! if required and therefore must be >= valid shape equivalents
      ! but for now we keep things the same. 
      actual_shape(1,1) = 1
      actual_shape(1,2) = 1
      actual_shape(1,3) = 1
      actual_shape(2,1) = valid_shape(2,1)     ! I (E-W) 
      actual_shape(2,2) = valid_shape(2,2)     ! J (S-N)
      actual_shape(2,3) = 1                    ! K (Vertical) 

      write(6,*)'grid64 - actual_shape ',actual_shape

      allocate(t_longitudes(actual_shape(1,1):actual_shape(2,1)))
      allocate(t_latitudes(actual_shape(1,2):actual_shape(2,2)))
      allocate(t_verticals(actual_shape(1,3):actual_shape(2,3)))


! work out where (in degrees lon, degrees lat, metres vrt) t points 
! are located
      do i = actual_shape(1,1),actual_shape(2,1)
        t_longitudes(i) = h_a_firstlong + (datastart(1)+i-2)*h_a_ewspace
      END DO


      ! Set the start latitudes for each PE 2 rows apart
      bottom_lat = -90.0 + (O4PE*2*h_a_nsspace) 
      do j = actual_shape(1,2),actual_shape(2,2)
         t_latitudes(j) = bottom_lat + (j-1)*h_a_nsspace
      END DO

      v_off64 = h_a_nsspace/2.0


!RSRH Set to have grid start at 0 in the vertical  - more intutive
! probably not important functionally
      do k = actual_shape(1,3), actual_shape(2,3)
        t_verticals(k) = 0.0  ! 
      END DO

      nbr_subdomains = 1 ! this is hard-wired into the UM?
      if (nGridDims .ne. Ndim_max) then
        write(6,*) "Error in OASIS4_grid64: "                           &
     &, "the grid that PRISM thinks it's on "                           &
     &, "is not the grid the UM thinks it's on."
      end if

      allocate(offset_array(Nfld_max,nbr_subdomains,nGridDims))
      allocate(extent_array(Nfld_max,nbr_subdomains,nGridDims))

      do i=1,nGridDims,2
         !! offset_array value is decreased by 1
         !! because PRISM _adds_ it to the starting location
         !! but UM uses it to specify the _exact_ starting location?
              offset_array(fld_type_p,nbr_subdomains,i) =               &
     & datastart_f(i,fld_type_p) - 1
              offset_array(fld_type_u,nbr_subdomains,i) =               &
     & datastart_f(i,fld_type_u) - 1
              offset_array(fld_type_v,nbr_subdomains,i) =               &
     & datastart_f(i,fld_type_v) - 1
              extent_array(fld_type_p,nbr_subdomains,i) =               &
     & blsize(i,fld_type_p)
              extent_array(fld_type_u,nbr_subdomains,i) =               &
     & blsize(i,fld_type_u)
              extent_array(fld_type_v,nbr_subdomains,i) =               &
     & blsize(i,fld_type_v)
      END DO

      ! Set special values for NS direction assuming 1xN 
      ! decomposition.
      i=2
         !! offset_array value is decreased by 1
         !! because PRISM _adds_ it to the starting location
         !! but UM uses it to specify the _exact_ starting location?
              offset_array(fld_type_p,nbr_subdomains,i) =               &
     & MYPE*2
              offset_array(fld_type_u,nbr_subdomains,i) =               &
     & MYPE*2 
              offset_array(fld_type_v,nbr_subdomains,i) =               &
     & MYPE*2

      IF (MYPE.LT.O4CNTLPE) THEN
              ! Each non O4CNTLPE cpu just gets 1 row of data 
              extent_array(fld_type_p,nbr_subdomains,i) =               &
     & 2
              extent_array(fld_type_u,nbr_subdomains,i) =               &
     & 2
              extent_array(fld_type_v,nbr_subdomains,i) =               &
     & 2
       ELSE
              extent_array(fld_type_p,nbr_subdomains,i) =               &
     & glsize(i,fld_type_p)-(O4CNTLPE*2)
              extent_array(fld_type_u,nbr_subdomains,i) =               &
     & glsize(i,fld_type_u)-(O4CNTLPE*2)
              extent_array(fld_type_v,nbr_subdomains,i) =               &
     & glsize(i,fld_type_v)-(O4CNTLPE*2)

       END IF

      ! Detect whether this PE is at the northern limit of our
      ! global domain.
      no_neighbour_N = at_extremity(Pnorth)

      ! Call the routine where our values will be used in the PRISM
      ! initialisation, grid definition and variable definitions. 
! DEPENDS ON: OASIS4_grid32_mod
      call grid32_OASIS4(comp_id,MODEL_BASIS_TIME,ierror                & 
     &, valid_shape,actual_shape,nGridDims,grid_type,nCorners           &
     &, t_longitudes,t_latitudes,t_verticals,v_off64                    &
     &, nbr_subdomains, offset_array, extent_array                      &
     &)

      return
      end subroutine OASIS4_grid64
#endif
