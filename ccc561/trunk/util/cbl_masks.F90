!jhan:Althoughthis isonly calling subrs and not using data - still needs revision to only call once and use L_tile_pts from here
!******************************************************************************
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
!******************************************************************************
MODULE cbl_masks_mod

  !-----------------------------------------------------------------------------
  ! Description:
  !   Computes various masks for CABLE:
  !     - veg_mask: TRUE where LAI > LAI threshold value
  !     - sunlit_mask: TRUE where the cosine of zenith angle is higher than a 
  !                    specified tolerance.
  !     - sunlit_veg_mask: TRUE where veg_mask and sunlit_mask are TRUE.
  !
  ! This MODULE is USEd in:
  !     cable_land_albedo_mod_cbl.F90 (JULES),
  !     cable_cbm.F90 (ESM1.5),
  !     cable_rad_driver.F90 (ESM1.5),
  !     cbl_model_driver_offline.F90 (CABLE)
  ! 
  ! This MODULE contains 3 public Subroutine:
  !     fveg_mask,
  !     fsunlit_mask,
  !     fsunlit_veg_mask
  !
  ! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
  ! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
  !-----------------------------------------------------------------------------
  
  public L_tile_pts
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
! Description:
!   Computes mask where LAI is higher than a given threshold.

implicit none
LOGICAL, allocatable :: veg_mask(:)
integer ::  mp
real :: lai_thresh
real :: reducedLAIdue2snow(mp)
integer :: i
 
IF ( .NOT. ALLOCATED(veg_mask)) ALLOCATE( veg_mask(mp) )

veg_mask = .FALSE.
! Define vegetation mask:
do i=1, mp  
  if( reducedLAIdue2snow(i) > lai_thresh ) veg_mask(i) = .true.
end do

End subroutine fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )
! Description:
!   Computes mask where the sun is above the horizon within a given tolerance.

implicit none
LOGICAL, allocatable :: sunlit_mask(:)
integer ::  mp
real :: coszen_tols
real :: coszen(mp)   
integer :: i

IF ( .NOT. ALLOCATED(sunlit_mask)) ALLOCATE( sunlit_mask(mp) )

sunlit_mask = .FALSE.
! Define sunlit mask:
do i=1, mp  
  if( coszen(i) > coszen_tols ) sunlit_mask(i) = .true.  
end do

End subroutine fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fsunlit_veg_mask( sunlit_veg_mask, mp )
! Description:
!   Computes mask where the sun is above the horizon within a given tolerance
!   and the LAI is above a given threshold. This is the union of the
!   masks from fsunlit_mask() and fveg_mask().

implicit none
LOGICAL, allocatable :: sunlit_veg_mask(:)
integer ::  mp
integer :: i

IF ( .NOT. ALLOCATED(sunlit_veg_mask) ) ALLOCATE( sunlit_veg_mask(mp) )

sunlit_veg_mask = .FALSE.
! Define sunlit AND vegetation mask:
do i=1, mp  
  if( veg_mask(i) .AND.  sunlit_mask(i) ) sunlit_veg_mask(i) = .true.
end do

End subroutine fsunlit_veg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


End module cbl_masks_mod
