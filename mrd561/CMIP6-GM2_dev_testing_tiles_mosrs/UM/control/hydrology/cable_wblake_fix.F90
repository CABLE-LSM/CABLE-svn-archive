module cable_wblake_mod
  implicit none

  real  :: wblake_ratio               &! ratio of wblake/subroff

  real, allocatable, save ::                                                         &
    WBLAKE(:), TOT_WBLAKE(:), TOT_SUBRUN(:)

  real, allocatable, save ::  WB_LAKE(:,:) ! 

contains
!      Real                                                              &
!     &  sub_surf_roff(land_points)                                      &
!                                    ! sub-surface runoff
!     &, snomlt_sub_htf(land_points)&! subsurface snowmelt heat flux
!     &, wblake_ratio               &! ratio of wblake/subroff
!     &, WBLAKE(land_points)        &! 
     &, SUBROFF(land_points)       &! 
!     &, SUBROFF0(land_points)       &! 
!     &, SUBROFF1(land_points)       &! 
     &, TOT_WBLAKE(land_points)    &! 
     &, TOT_SUBRUN(land_points)   &! 
     &, WB_LAKE(land_points,ntiles) ! 

subroutine cable_wblake_fix_alloc( land_points, ntiles )
  integer :: land_points, ntiles
  logical, save :: first_call==.true.

  if(first_call) then
    allocate ( WBLAKE(land_points),           &! 
               TOT_WBLAKE(land_points),       &! 
               TOT_SUBRUN(land_points),       &! 
               WB_LAKE(land_points,ntiles) )   ! 
  endif

  wblake  = 0.0
  subroff = 0.0

  first_call = .FALSE.               
End subroutine cable_wblake_fix_alloc

End module cable_wblake_mod
