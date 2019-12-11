module cbl_masks_mod
  public L_tile_pts,   mask_active_tile_pts
  public fveg_mask
  public fsunlit_mask
  public fsunlit_veg_mask
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable, save :: L_tile_pts(:,:)
  LOGICAL, allocatable, save :: mask_active_tile_pts(:,:)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )

implicit none
integer ::  mp
real :: lai_thresh
logical :: veg_mask(mp)
real :: reducedLAIdue2snow(mp)
 
  veg_mask = reducedLAIdue2snow > lai_thresh 

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )

implicit none
integer ::  mp
real :: coszen_tols
logical :: sunlit_mask(mp)
real :: coszen(mp)   

  ! Define sunlit mask:
  sunlit_mask =  coszen(:) > coszen_tols
   
End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( sunlit_veg_mask, mp, veg_mask, sunlit_mask)

implicit none
integer ::  mp
logical :: veg_mask(mp)
logical :: sunlit_mask(mp)
logical :: sunlit_veg_mask(mp)

   ! Define sunlit AND vegetation mask:
   sunlit_veg_mask = veg_mask .AND.  sunlit_mask

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
