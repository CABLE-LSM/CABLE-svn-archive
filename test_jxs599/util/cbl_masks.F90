!jhan:Althoughthis isonly calling subrs and not using data - still needs revision to only call once and use L_tile_pts from here
module cbl_masks_mod
  public L_tile_pts
  public init_active_tile_mask_cbl
  public fveg_mask
  public fsunlit_mask
  public fsunlit_veg_mask
  public veg_mask
  public sunlit_mask
  public sunlit_veg_mask
!H! remove SAVE attr later 
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable,SAVE :: L_tile_pts(:,:)
  LOGICAL, allocatable, SAVE :: veg_mask(:), sunlit_mask(:), sunlit_veg_mask(:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )

implicit none
LOGICAL, allocatable :: veg_mask(:)
integer ::  mp
real :: lai_thresh
real :: reducedLAIdue2snow(mp)
 
IF ( .NOT. ALLOCATED(veg_mask)) ALLOCATE( veg_mask(mp) )

  veg_mask = reducedLAIdue2snow > lai_thresh 

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )

implicit none
LOGICAL, allocatable :: sunlit_mask(:)
integer ::  mp
real :: coszen_tols
real :: coszen(mp)   

  IF ( .NOT. ALLOCATED(sunlit_mask)) ALLOCATE( sunlit_mask(mp) )

  ! Define sunlit mask:
  sunlit_mask =  coszen(:) > coszen_tols
   
End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( sunlit_veg_mask, mp )

implicit none
LOGICAL, allocatable :: sunlit_veg_mask(:)
integer ::  mp

  IF ( .NOT. ALLOCATED(sunlit_veg_mask)) ALLOCATE( sunlit_veg_mask(mp) )

  ! Define sunlit AND vegetation mask:
  sunlit_veg_mask = veg_mask .AND.  sunlit_mask

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
