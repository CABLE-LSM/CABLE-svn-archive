module cbl_masks_mod
  
  !mask TRUE where tile fraction is greater than zero
  LOGICAL, allocatable, save :: L_tile_pts(:,:)
  LOGICAL, allocatable, save :: mask_active_tile_pts(:,:)

  LOGICAL, allocatable, save :: veg_mask(:) 
  LOGICAL, allocatable, save :: sunlit_mask(:) 
  LOGICAL, allocatable, save :: sunlit_veg_mask(:) 

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
subroutine fmask_active_tile_pts( land_pts, ntiles, tile_frac )
implicit none
integer :: land_pts, ntiles
real :: tile_frac(land_pts,ntiles)      ! surface type fraction 
integer :: i,j 

!---def. vector length for cable(mp) & logical l_tile_pts
!--- IF the tile is "active"
if( .NOT. allocated(L_tile_pts) ) allocate ( L_tile_pts(land_pts,ntiles) ) 
L_TILE_PTS = .FALSE.

DO i=1,land_pts
  DO j=1,ntiles
      
    IF( TILE_FRAC(i,j) .GT. 0.0 ) THEN 
      L_TILE_PTS(i,j) = .TRUE.
    ENDIF
   
  ENDDO
ENDDO

mask_active_tile_pts = L_TILE_PTS

End subroutine fmask_active_tile_pts 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fveg_mask( mp, reducedLAIdue2snow )

Use cable_other_constants_mod, ONLY : lai_thresh

implicit none
integer ::  mp
real :: reducedLAIdue2snow(mp)

  if( .NOT. allocated(veg_mask) ) &
    allocate ( veg_mask(mp) ) 

  ! Define vegetated mask:
  veg_mask = reducedLAIdue2snow > lai_thresh 

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( mp, nrb, coszen_tols, coszen )
!subroutine fsunlit_mask( mp, nrb, sw_down )
Use cable_other_constants_mod, ONLY : rad_thresh

implicit none
integer ::  mp
integer ::  nrb
real :: coszen_tols
real :: coszen(mp)   

  if( .NOT. allocated(sunlit_mask) ) &
    allocate ( sunlit_mask(mp) ) 

  ! Define sunlit mask:
  sunlit_mask =  coszen(:) > coszen_tols
   
End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( mp, veg_mask_IN, sunlit_mask_IN)

implicit none
integer ::  mp
logical :: veg_mask_IN(mp)
logical :: sunlit_mask_IN(mp)

  if( .NOT. allocated(sunlit_veg_mask) ) &
    allocate ( sunlit_veg_mask(mp) ) 

   ! Define sunlit AND vegetation mask:
   sunlit_veg_mask = veg_mask_IN .AND.  sunlit_mask_IN

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
