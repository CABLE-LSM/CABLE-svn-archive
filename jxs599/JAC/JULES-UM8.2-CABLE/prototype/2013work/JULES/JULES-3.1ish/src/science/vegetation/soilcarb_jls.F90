#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine SOILCARB -----------------------------------------------
!
! Purpose : Updates carbon contents of the soil.
!
! -------------------------------------------------------------------
SUBROUTINE soilcarb (land_pts,trif_pts,trif_index                 &
,                    lit_c,frac,frac_agric,resp_frac              &
,                    forw,GAMMA,lit_c_t,resp_s,cs)

USE nstypes

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
,l,t                        ! WORK Loop counters

REAL                                                              &
 forw                                                             &
                            ! IN Forward timestep weighting.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,lit_c_t(land_pts)                                                &
                            ! IN Total carbon litter
!                                 !    (kg C/m2/360days).
,lit_c(land_pts,npft)                                             &
                            ! IN Carbon Litter (kg C/m2/360days).
,resp_s(land_pts,5)                                               &
                            ! INOUT Soil respiration
!                                 !    (kg C/m2/360days).
,cs(land_pts,4)                                                   &
                            ! INOUT Soil carbon (kg C/m2).
!                                 !    The 4 soil C pools are
!                                 !    DPM, RPM, biomass and humus.
,dcs(land_pts,4)                                                  &
                            ! WORK Increment to the soil carbon
!                                 !      (kg C/m2).
,dpc_dcs(land_pts,4)                                              &
                            ! WORK Rate of change of PC with
!                                 !      soil carbon (/360days).
,pc(land_pts,4)                                                   &
                            ! WORK Net carbon accumulation in
!                                 !      the soil (kg C/m2/360days).
,frac(land_pts,ntype)                                             &
                            ! INOUT Fractional cover of each
,frac_agric(land_pts)                                             &
                            ! INOUT Fractional cover of each
,dpm_ratio(land_pts)                                              &
                            ! WORK ratio of litter input into DPM
,resp_frac(land_pts)        ! respired fraction of RESP_S

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOILCARB',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Note: The 4 soil carbon pools are  1 decomposable plant material,
! 2 resistant plant material, 3 biomass, 4 humus.
!----------------------------------------------------------------------
! sum total respiration
DO l=1,land_pts
  resp_s(l,5) = resp_s(l,1) + resp_s(l,2) +                       &
                resp_s(l,3) + resp_s(l,4)
END DO

! calculate DPM:RPM ratio of input litter Carbon
! DEPENDS ON: dpm_rpm
CALL dpm_rpm(land_pts,trif_pts,trif_index,frac,                   &
             frac_agric,lit_c,dpm_ratio)

DO t=1,trif_pts
  l=trif_index(t)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the soil
!----------------------------------------------------------------------
  pc(l,1) = dpm_ratio(l)*lit_c_t(l) - resp_s(l,1)
  pc(l,2) = (1-dpm_ratio(l))*lit_c_t(l) - resp_s(l,2)
  pc(l,3) = 0.46*resp_frac(l)*resp_s(l,5) - resp_s(l,3)
  pc(l,4) = 0.54*resp_frac(l)*resp_s(l,5) - resp_s(l,4)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
  dpc_dcs(l,1) = resp_s(l,1)/cs(l,1)
  dpc_dcs(l,2) = resp_s(l,2)/cs(l,2)
  dpc_dcs(l,3) = resp_s(l,3)/cs(l,3)
  dpc_dcs(l,4) = resp_s(l,4)/cs(l,4)

!----------------------------------------------------------------------
! Save current value of soil carbon
!----------------------------------------------------------------------
  dcs(l,1) = cs(l,1)
  dcs(l,2) = cs(l,2)
  dcs(l,3) = cs(l,3)
  dcs(l,4) = cs(l,4)

END DO


!----------------------------------------------------------------------
! Update soil carbon
!----------------------------------------------------------------------
! DEPENDS ON: decay
CALL decay (land_pts,trif_pts,trif_index                          &
,           dpc_dcs,forw,GAMMA,pc,cs)

!----------------------------------------------------------------------
! Apply implicit correction to the soil respiration rate.
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)

  dcs(l,1) = cs(l,1) - dcs(l,1)
  dcs(l,2) = cs(l,2) - dcs(l,2)
  dcs(l,3) = cs(l,3) - dcs(l,3)
  dcs(l,4) = cs(l,4) - dcs(l,4)

  resp_s(l,1) = resp_s(l,1) + forw*dpc_dcs(l,1)*dcs(l,1)
  resp_s(l,2) = resp_s(l,2) + forw*dpc_dcs(l,2)*dcs(l,2)
  resp_s(l,3) = resp_s(l,3) + forw*dpc_dcs(l,3)*dcs(l,3)
  resp_s(l,4) = resp_s(l,4) + forw*dpc_dcs(l,4)*dcs(l,4)

END DO

! sum total respiration
DO l=1,land_pts
  resp_s(l,5) = resp_s(l,1) + resp_s(l,2) +                       &
                resp_s(l,3) + resp_s(l,4)
END DO

IF (lhook) CALL dr_hook('SOILCARB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE soilcarb
#endif
