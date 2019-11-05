module cable_rad_unpack_mod
  
contains
 
subroutine cable_rad_unpack( &
#                            include "cbl_rad_unpack_args.inc"
                           )
!data  
USE cbl_allocate_types_mod, ONLY : rad

  !diag 
  USE cable_fprint_module, ONLY : cable_fprintf
  USE cable_Pyfprint_module, ONLY : cable_Pyfprintf
  USE cable_fFile_module, ONLY : fprintf_dir_root, fprintf_dir, L_cable_fprint,&
                                 L_cable_Pyfprint, unique_subdir
implicit none

!___ re-decl input args
#                            include "cbl_rad_unpack_decs.inc"

!___ local vars
INTEGER :: i,J,K,L,N
REAL :: miss = 0.0
 
! std template args 
character(len=*), parameter :: subr_name = "cable_rad_unpack"
real :: Sumreffbm(mp) 
real :: Sumreffdf(mp) 
  
# include "../../../util/cable/diag/cable_fprint.txt"
! only for land points, at present do not have a method for treating 
! mixed land/sea or land/seaice points as yet.
alb_surft(:,:,:) = 0.0
alb_surft(:,:,1) = UNPACK(EffSurfRefl_beam(:,1),L_TILE_PTS, miss)
alb_surft(:,:,3) = UNPACK(EffSurfRefl_beam(:,2),L_TILE_PTS, miss)
alb_surft(:,:,2) = UNPACK(EffSurfRefl_dif(:,1),L_TILE_PTS, miss)
alb_surft(:,:,4) = UNPACK(EffSurfRefl_dif(:,2),L_TILE_PTS, miss)

do i=1,land_pts
  do j=1,nsurft
              
  if( alb_surft(i,j,1)> 1.) then
    print *, 'alb > 1',alb_surft(i,j,1)
    STOP
  ENDIF
 
  ENDDO
ENDDO

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

Sumreffbm = 0.0
Sumreffdf = 0.0
do i=1, mp
Sumreffbm(i) = Sumreffbm(i) + ( EffSurfRefl_beam(i,1)+EffSurfRefl_beam(i,2) )
Sumreffdf(i) = Sumreffdf(i) + ( EffSurfRefl_dif(i,1)+EffSurfRefl_dif(i,2) )
enddo

fprintf_dir="/home/599/jxs599/"
vname='HEffSRefl_beam'; dimx=mp 
call cable_Pyfprintf( cDiag1, vname, Sumreffbm, dimx, .true.)
vname='HEffSRefl_dif'; dimx=mp 
call cable_Pyfprintf( cDiag2, vname, Sumreffdf, dimx, .true.)

return
     
End subroutine cable_rad_unpack
End module cable_rad_unpack_mod
