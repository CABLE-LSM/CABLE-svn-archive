#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine VEGCARB ------------------------------------------------
!
! Purpose : Updates carbon contents of the vegetation.
!
! -------------------------------------------------------------------
 SUBROUTINE vegcarb (land_pts,trif_pts,trif_index                 &
,                    n,forw,GAMMA,g_leaf,npp,resp_w               &
,                    leaf,root,wood,dcveg,pc_s)

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
                            ! IN Plant functional type.
,l,t                        ! WORK Loop counters

REAL                                                              &
 forw                                                             &
                            ! IN Forward timestep weighting.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,g_leaf(land_pts)                                                 &
                            ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
,npp(land_pts)                                                    &
                            ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
,resp_w(land_pts)                                                 &
                            ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
,leaf(land_pts)                                                   &
                            ! INOUT Leaf biomass (kg C/m2).
,root(land_pts)                                                   &
                            ! INOUT Root biomass (kg C/m2).
,wood(land_pts)                                                   &
                            ! INOUT Woody biomass (kg C/m2).
,dcveg(land_pts)                                                  &
                            ! OUT Change in vegetation carbon
!                                 !     during the timestep
!                                 !     (kg C/m2/timestep).
,pc_s(land_pts)                                                   &
                            ! OUT Net carbon flux available
!                                 !     for spreading (kg C/m2/360days).
,dfpar_dlai                                                       &
                            ! WORK Rate of change of FPAR
!                                 !      with leaf area index.
,dlai                                                             &
                            ! WORK Increment to the leaf area
!                                 !      index.
,dlamg_dlai,dlit_dlai                                             &
                            ! WORK Required for calculation
!                                 !      of the equilibrium increments.
,dnpp_dlai(land_pts)                                              &
                            ! WORK Rate of change of NPP
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,dpc_dlai(land_pts)                                               &
                            ! WORK Rate of change of PC
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,dpcg_dlai(land_pts)                                              &
                            ! WORK Rate of change of PC_G
!                                 !      with leaf area index
!                                 !      (kg C/m2/360days/LAI).
,drespw_dlai                                                      &
                            ! WORK Rate of change of RESP_W
!                                 !      with leaf area index
,fpar                                                             &
                            ! WORK PAR interception factor.
,lai(land_pts)                                                    &
                            ! WORK Leaf area index.
,lambda_g                                                         &
                            ! WORK Fraction of NPP available
!                                 !      for spreading.
,lit_c_l(land_pts)                                                &
                            ! WORK Local rate of Carbon Litter
!                                 !      production (kg C/m2/360days).
,pc(land_pts)                                                     &
                            ! WORK Net carbon flux available
!                                 !      to vegetation (kg C/m2/360days)
,pc_g(land_pts)             ! WORK Net carbon flux available
!                                 !      for growth (kg C/m2/360days).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('VEGCARB',zhook_in,zhook_handle)

DO t=1,trif_pts
  l=trif_index(t)

  lai(l) = leaf(l)/sigl(n)
!----------------------------------------------------------------------
! Calculate the local production rate for carbon litter
!----------------------------------------------------------------------
  lit_c_l(l) = g_leaf(l)*leaf(l)+g_root(n)*root(l)                &
               + g_wood(n)*wood(l)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the vegetation
!----------------------------------------------------------------------
  pc(l) = npp(l) - lit_c_l(l)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
  dlit_dlai = (g_leaf(l)*leaf(l)+g_root(n)*root(l))/lai(l)        &
            + b_wl(n)*g_wood(n)*wood(l)/lai(l)

  fpar = (1 - EXP(-kpar(n)*lai(l)))/kpar(n)
  dfpar_dlai = EXP(-kpar(n)*lai(l))

  dnpp_dlai(l) = npp(l)*dfpar_dlai/fpar                           &
               + (1-r_grow(n))*resp_w(l)                          &
               *(dfpar_dlai/fpar-b_wl(n)/lai(l))

  lambda_g = 1 - (lai(l) - lai_min(n))                            &
                /(lai_max(n) - lai_min(n))
  dlamg_dlai = -1.0/(lai_max(n) - lai_min(n))

  pc_g(l) = lambda_g * npp(l) - lit_c_l(l)
  dpcg_dlai(l) = lambda_g*dnpp_dlai(l)                            &
               + dlamg_dlai*npp(l)                                &
               - dlit_dlai
  dpc_dlai(l) = dnpp_dlai(l) - dlit_dlai

END DO

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  dcveg(l) = leaf(l)+root(l)+wood(l)
END DO

! DEPENDS ON: growth
CALL growth (land_pts,trif_pts,trif_index                         &
,            n,dpcg_dlai,forw,GAMMA,pc_g,leaf,root,wood)

DO t=1,trif_pts
  l=trif_index(t)
  dcveg(l) = leaf(l)+root(l)+wood(l)-dcveg(l)
END DO

!----------------------------------------------------------------------
! Diagnose the carbon available for spreading and apply implicit
! corrections to the driving fluxes.
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  dlai = leaf(l)/sigl(n) - lai(l)
  pc_s(l) = pc(l) + forw*dpc_dlai(l)*dlai - dcveg(l)*GAMMA

  fpar = (1 - EXP(-kpar(n)*lai(l)))/kpar(n)
  dfpar_dlai = EXP(-kpar(n)*lai(l))
  drespw_dlai = resp_w(l)*b_wl(n)/lai(l)

  npp(l) = npp(l) + forw*dnpp_dlai(l)*dlai
  resp_w(l) = resp_w(l) + forw*drespw_dlai*dlai
END DO

IF (lhook) CALL dr_hook('VEGCARB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vegcarb
#endif
