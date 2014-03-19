#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine TRIFFID ------------------------------------------------
!
!                     Top-down
!                     Representation of
!                     Interactive
!                     Foliage and
!                     Flora
!                     Including
!                     Dynamics
!
! Purpose : Simulates changes in vegetation structure, areal
!           coverage and the carbon contents of vegetation and soil.
!           can be used to advance these variables dynamically
!           (GAMMA=1/TIMESTEP) or to iterate towards  equilibrium
!           (GAMMA --> 0.0, FORW=1.0).
!
! -------------------------------------------------------------------
SUBROUTINE triffid (land_pts,trif_pts,trif_index,forw,GAMMA       &
,                   frac_vs,frac_agric,g_leaf,npp,resp_s          &
,                   resp_w,cs,frac,ht,lai,resp_frac               &
,                   c_veg,cv,lit_c,lit_c_t)

USE nstypes
USE pftparm
USE trif

USE switches, ONLY: l_veg_compete

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
,l,n,t                      ! WORK Loop counters

REAL                                                              &
 forw                                                             &
                            ! IN Forward timestep weighting.
,frac_vs(land_pts)                                                &
                            ! IN Total fraction of gridbox
!                                 !    covered by veg or soil.
,GAMMA                                                            &
                            ! IN Inverse timestep (/360days).
,frac_agric(land_pts)                                             &
                            ! IN Fraction of agriculture.
,g_leaf(land_pts,npft)                                            &
                            ! IN Turnover rate for leaf and
!                                 !    fine root biomass (/360days).
,npp(land_pts,npft)                                               &
                            ! INOUT Net primary productivity
!                                 !       (kg C/m2/360days).
,resp_s(land_pts,5)                                               &
                            ! INOUT Soil respiration
!                                 !       (kg C/m2/360days).
,resp_w(land_pts,npft)                                            &
                            ! INOUT Wood maintenance respiration
!                                 !       (kg C/m2/360days).
,cs(land_pts,4)                                                   &
                            ! INOUT Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                             &
                            ! INOUT Fractional cover of each
!                                 !       Functional Type.
,ht(land_pts,npft)                                                &
                            ! INOUT Vegetation height (m).
,lai(land_pts,npft)                                               &
                            ! INOUT Leaf area index.
,c_veg(land_pts,npft)                                             &
                            ! OUT Total carbon content of
!                                 !     the vegetation (kg C/m2).
,cv(land_pts)                                                     &
                            ! OUT Gridbox mean vegetation
!                                 !     carbon (kg C/m2).
,lit_c(land_pts,npft)                                             &
                            ! OUT Carbon Litter (kg C/m2/360days).
,lit_c_t(land_pts)          ! OUT Gridbox mean carbon litter
!                                 !     (kg C/m2/360days).

REAL                                                              &
 dcveg(land_pts,npft)                                             &
                            ! WORK Change in vegetation carbon
!                                 !      during the timestep
!                                 !      (kg C/m2/timestep).
,dfrac(land_pts,npft)                                             &
                            ! WORK Change in areal fraction
!                                 !      during the timestep
!                                 !      (/timestep).
,frac_flux                                                        &
                            ! WORK PFT fraction to be used
!                                 !      in the calculation of
!                                 !      the gridbox mean fluxes.
,resp_frac(land_pts)                                              &
                            ! WORK  respired fraction of RESP_S
,lai_bal(land_pts,npft)                                           &
                            ! WORK Leaf area index in balanced
!                                 !      growth state.
,leaf(land_pts,npft)                                              &
                            ! WORK Leaf biomass (kg C/m2).
,pc_s(land_pts,npft)                                              &
                            ! WORK Net carbon flux available
!                                 !      for spreading
!                                 !      (kg C/m2/yr).
,phen(land_pts,npft)                                              &
                            ! WORK Phenological state.
,root(land_pts,npft)                                              &
                            ! WORK Root biomass (kg C/m2).
,wood(land_pts,npft)        ! WORK Woody biomass (kg C/m2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Loop through Functional Types
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('TRIFFID',zhook_in,zhook_handle)

DO n=1,npft

!----------------------------------------------------------------------
! Loop through TRIFFID points
!----------------------------------------------------------------------
  DO t=1,trif_pts
    l=trif_index(t)

!----------------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the associated leaf,
! wood, root and total vegetation carbon
!----------------------------------------------------------------------
    lai_bal(l,n) = (a_ws(n)*eta_sl(n)*ht(l,n)                     &
              /a_wl(n))**(1.0/(b_wl(n)-1))
    leaf(l,n) = sigl(n)*lai_bal(l,n)
    root(l,n) = leaf(l,n)
    wood(l,n) = a_wl(n)*(lai_bal(l,n)**b_wl(n))
    c_veg(l,n) = leaf(l,n) + root(l,n) + wood(l,n)
!----------------------------------------------------------------------
! Diagnose the phenological state
!----------------------------------------------------------------------
    phen(l,n) = lai(l,n)/lai_bal(l,n)

  END DO

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
! DEPENDS ON: vegcarb
  CALL vegcarb (land_pts,trif_pts,trif_index,n,forw               &
,               GAMMA,g_leaf(:,n),npp(:,n),resp_w(:,n)            &
,               leaf(:,n),root(:,n),wood(:,n)                     &
,               dcveg(:,n),pc_s(:,n))

END DO

!-----------------------------------------------------------------------
! Diagnose the new value of Canopy Height, Leaf Area Index and Total
! Vegetation Carbon
!-----------------------------------------------------------------------
DO n=1,npft

  DO t=1,trif_pts
    l=trif_index(t)

    ht(l,n) = wood(l,n) / (a_ws(n) * eta_sl(n))                   &
            * (a_wl(n)/wood(l,n))**(1.0/b_wl(n))
    lai_bal(l,n) = leaf(l,n) / sigl(n)
    lai(l,n) = phen(l,n) * lai_bal(l,n)
    c_veg(l,n) = leaf(l,n) + root(l,n) + wood(l,n)

  END DO

END DO

!----------------------------------------------------------------------
! Update the areal coverage of each functional type
!----------------------------------------------------------------------
IF ( l_veg_compete ) THEN
! DEPENDS ON: lotka
  CALL lotka (land_pts,trif_pts,trif_index                          &
,             c_veg,forw,frac_vs,frac_agric,GAMMA,lai_bal,pc_s      &
,             frac,dfrac)
ELSE
  dfrac(:,:) = 0.0
END IF

!----------------------------------------------------------------------
! Diagnose the litterfall from the carbon balance of each vegetation
! type
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)

  lit_c_t(l) = 0.0

  DO n=1,npft
    frac_flux=frac(l,n)-(1.0-forw)*dfrac(l,n)
    lit_c(l,n) = npp(l,n)-GAMMA/frac_flux*(c_veg(l,n)*frac(l,n)   &
               -(c_veg(l,n)-dcveg(l,n))*(frac(l,n)-dfrac(l,n)))
    lit_c_t(l) = lit_c_t(l)+frac_flux*lit_c(l,n)
  END DO
END DO

!----------------------------------------------------------------------
! Call SOIL_C to update the soil carbon content
!----------------------------------------------------------------------
! DEPENDS ON: soilcarb
CALL soilcarb (land_pts,trif_pts,trif_index                       &
,              lit_c,frac,frac_agric,resp_frac                    &
,              forw,GAMMA,lit_c_t,resp_s,cs)

!----------------------------------------------------------------------
! Diagnose the gridbox mean vegetation carbon
!----------------------------------------------------------------------
DO t=1,trif_pts
  l=trif_index(t)
  cv(l) = 0.0
  DO n=1,npft
    cv(l) = cv(l)+frac(l,n)*c_veg(l,n)
  END DO
END DO

IF (lhook) CALL dr_hook('TRIFFID',zhook_out,zhook_handle)
RETURN
END SUBROUTINE triffid
#endif
