!jhan:Althoughthis isonly calling subrs and not using data - still needs revision to only call once and use L_tile_pts from here
module cbl_masks_mod
  public L_tile_pts,   mask_active_tile_pts
  public fveg_mask
  public fsunlit_mask
  public fsunlit_veg_mask
  public veg_mask
  public sunlit_mask
  public sunlit_veg_mask
  
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable :: L_tile_pts(:,:)
  LOGICAL, allocatable :: veg_mask(:), sunlit_mask(:), sunlit_veg_mask(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_active_tile_mask_cbl()

USE cable_types_mod,          ONLY: l_tile_pts_types => l_tile_pts
USE ancil_info,               ONLY: frac_surft, land_pts
USE jules_surface_types_mod,  ONLY: ntype

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

INTEGER :: i, j

! Determine active tiles map
IF ( .NOT. ALLOCATED(l_tile_pts)) ALLOCATE( l_tile_pts(land_pts, ntype) )

l_tile_pts(:,:) = .FALSE.

DO j = 1, ntype
  DO i = 1, land_pts
    IF ( frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO 

l_tile_pts_types = l_tile_pts

RETURN

END SUBROUTINE init_active_tile_mask_cbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fveg_mask( mp, lai_thresh, reducedLAIdue2snow )

implicit none
integer ::  mp
real :: lai_thresh
real :: reducedLAIdue2snow(mp)
 
IF ( .NOT. ALLOCATED(veg_mask)) ALLOCATE( veg_mask(mp) )

  veg_mask = reducedLAIdue2snow > lai_thresh 

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( mp, coszen_tols, coszen )

implicit none
integer ::  mp
real :: coszen_tols
real :: coszen(mp)   

  IF ( .NOT. ALLOCATED(sunlit_mask)) ALLOCATE( sunlit_mask(mp) )

  ! Define sunlit mask:
  sunlit_mask =  coszen(:) > coszen_tols
   
End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( mp )

implicit none
integer ::  mp

  IF ( .NOT. ALLOCATED(sunlit_veg_mask)) ALLOCATE( sunlit_veg_mask(mp) )

  ! Define sunlit AND vegetation mask:
  sunlit_veg_mask = veg_mask .AND.  sunlit_mask

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
