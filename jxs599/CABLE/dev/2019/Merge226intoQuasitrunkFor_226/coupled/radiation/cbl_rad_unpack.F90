module cable_rad_unpack_mod
  
contains
 
subroutine cable_rad_unpack( &
#                            include "cbl_rad_unpack_args.inc"
                           )
!data  
USE cbl_allocate_types_mod, ONLY : rad

implicit none

!___ re-decl input args
#                            include "cbl_rad_unpack_decs.inc"

!___ local vars
INTEGER :: i,J,K,L,N
REAL :: miss = 0.0
 
! std template args 
character(len=*), parameter :: subr_name = "cable_rad_unpack"
  
! only for land points, at present do not have a method for treating 
! mixed land/sea or land/seaice points as yet.
alb_surft(:,:,1) = UNPACK(EffSurfRefl_beam(:,1),L_TILE_PTS, miss)
alb_surft(:,:,3) = UNPACK(EffSurfRefl_beam(:,2),L_TILE_PTS, miss)
alb_surft(:,:,2) = UNPACK(EffSurfRefl_dif(:,1),L_TILE_PTS, miss)
alb_surft(:,:,4) = UNPACK(EffSurfRefl_dif(:,2),L_TILE_PTS, miss)

LAND_ALBEDO=0

DO N=1,nsurft
  DO K=1,TILE_PTS(N)
    L = TILE_INDEX(K,N)
    J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
    I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
    
    ! direct beam visible
    LAND_ALBEDO(I,J,1) = LAND_ALBEDO(I,J,1) +              &
                               TILE_FRAC(L,N)*ALB_surft(L,N,1)

    ! diffuse beam visible
    LAND_ALBEDO(I,J,2) = LAND_ALBEDO(I,J,2) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,2)

    ! direct beam nearinfrared 
    LAND_ALBEDO(I,J,3) = LAND_ALBEDO(I,J,3) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,3)

    ! diffuse beam nearinfrared
    LAND_ALBEDO(I,J,4) = LAND_ALBEDO(I,J,4) +               &
                               TILE_FRAC(L,N)*ALB_surft(L,N,4)
  ENDDO
ENDDO

return
     
End subroutine cable_rad_unpack
End module cable_rad_unpack_mod
