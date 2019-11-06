module cable_wide_mod 

!CLAI_thresh

public

!already passed at rad call 
real, allocatable, save :: SnowDepth(:)             !Formerly: ssnow%snowd 
real, allocatable, save :: SnowDensity(:)           !Formerly: ssnow%ssdnn 

real, allocatable, save :: SnowODepth(:)             !Formerly: ssnow%osnowd
 
!computed from UM HT(LAI)_PFT passed in explicit call - need at rad call
real, allocatable, save :: LAI_pft_cbl(:)           !Formerly: ~veg%vlai
real, allocatable, save :: HGT_pft_cbl(:)           !Formerly:  ~veg%hc 

!can compute from  z0surf_min, HGT_pft_cbl, SnowDepth, SnowDensity )
real, allocatable, save :: HeightAboveSnow(:)       !Formerly: rough%hruff

! LAI_pft_cbl, HGT_pft_cbl, HeightAboveSnow
!real, allocatable, save :: reducedLAIdue2snow(:)  !Formerly: canopy%vlaiw 

contains

subroutine allocate_cable_wide( mp )

implicit none
integer :: mp

!if(.NOT. allocated(reducedLAIdue2snow)) allocate( reducedLAIdue2snow(mp) )
if(.NOT. allocated(HeightAboveSnow)) allocate( HeightAboveSnow(mp) )
if(.NOT. allocated(SnowDepth)) allocate( SnowDepth(mp) )
if(.NOT. allocated(SnowODepth)) allocate( SnowODepth(mp) )
if(.NOT. allocated(SnowDensity)) allocate( SnowDensity(mp) )
if(.NOT. allocated(LAI_pft_cbl)) allocate( LAI_pft_cbl(mp) )
if(.NOT. allocated(HGT_pft_cbl)) allocate( HGT_pft_cbl(mp) )

End subroutine allocate_cable_wide

End module cable_wide_mod 
