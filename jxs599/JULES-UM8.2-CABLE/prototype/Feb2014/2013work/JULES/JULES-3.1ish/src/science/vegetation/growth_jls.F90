#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine GROWTH -------------------------------------------------
!
! Purpose : Increments leaf, root and wood carbon.
!
! -------------------------------------------------------------------
SUBROUTINE growth (land_pts,trif_pts,trif_index                   &
,                  n,dpcg_dlai,forw,GAMMA,pc_g                    &
,                  leaf,root,wood)

USE descent
USE pftparm
USE trif

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 land_pts                                                         &
                            ! IN Total number of land points.
,trif_pts                                                         &
                            ! IN Number of points on which
!                                 !    TRIFFID may operate
,trif_index(land_pts)                                             &
                            ! IN Indices of land points on
!                                 !    which TRIFFID may operate
,n                                                                &
                            ! IN Vegetation type.
,l,t                        ! WORK Loop counters

REAL                                                              &
 dpcg_dlai(land_pts)                                              &
                            ! IN Rate of change of PC_G with
!                                 !    leaf area index
!                                 !    (kg C/m2/360days/LAI).
,forw                                                             &
                            ! IN Forward timestep weighting.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,pc_g(land_pts)                                                   &
                            ! IN Net carbon flux available
!                                 !    for growth (kg C/m2/360days).
,leaf(land_pts)                                                   &
                            ! INOUT Leaf biomass (kg C/m2).
,root(land_pts)                                                   &
                            ! INOUT Root biomass (kg C/m2).
,wood(land_pts)             ! INOUT Woody biomass (kg C/m2).

REAL                                                              &
 denom                                                            &
                            ! WORK Denominator of update
!                                 !      equation.
,dleaf,droot,dwood                                                &
                            ! WORK Increments to leaf, root
!                                 !      and woody biomass (kg C/m2).
,dl_dw                                                            &
                            ! WORK Rate of change of leaf
!                                 !      carbon with wood carbon.
,dlai_dw                                                          &
                            ! WORK Rate of change of leaf area
!                                 !      index with wood carbon
!                                 !      (LAI m2/kg C).
,dr_dw                                                            &
                            ! WORK Rate of change of root
!                                 !      carbon with wood carbon.
,numer                                                            &
                            ! WORK Numerator of the update
!                                 !      equation.
,wood_max                                                         &
                            ! WORK Maximum wood carbon (kg C/m2).
,wood_min                   ! WORK Minimum wood carbon (kg C/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('GROWTH',zhook_in,zhook_handle)

DO t=1,trif_pts
  l=trif_index(t)

!----------------------------------------------------------------------
! Calculate the increment to the wood carbon
!----------------------------------------------------------------------
  dl_dw = leaf(l)/(b_wl(n)*wood(l))
  dr_dw = dl_dw
  dlai_dw = dl_dw/sigl(n)

  numer = pc_g(l)
  denom = (1+dl_dw+dr_dw)*GAMMA-forw*dlai_dw*dpcg_dlai(l)
  denom = MAX(denom,denom_min)

  dwood = numer/denom

!----------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!----------------------------------------------------------------------
  wood_min = a_wl(n)*lai_min(n)**b_wl(n)
  wood_max = a_wl(n)*lai_max(n)**b_wl(n)
  dwood = MAX((wood_min-wood(l)),dwood)
  dwood = MIN((wood_max-wood(l)),dwood)

!----------------------------------------------------------------------
! Diagnose the increments to leaf and root carbon
!----------------------------------------------------------------------
  dleaf = sigl(n)*((wood(l)+dwood)/a_wl(n))**(1.0/b_wl(n))        &
         -leaf(l)
  droot = dleaf

!----------------------------------------------------------------------
! Update carbon contents
!----------------------------------------------------------------------
  leaf(l) = leaf(l)+dleaf
  root(l) = root(l)+droot
  wood(l) = wood(l)+dwood

END DO

IF (lhook) CALL dr_hook('GROWTH',zhook_out,zhook_handle)
RETURN
END SUBROUTINE growth
#endif
