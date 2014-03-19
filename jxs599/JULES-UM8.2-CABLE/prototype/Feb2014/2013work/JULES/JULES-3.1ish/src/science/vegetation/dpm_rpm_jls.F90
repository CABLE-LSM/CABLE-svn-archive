#if defined(L19_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine DPM_RPM ------------------------------------------------
!
! Purpose : Calculates the DPM_RPM ratio of litter input for use
!            in the RothC soil carbon sub-model
!
! -------------------------------------------------------------------

SUBROUTINE dpm_rpm(land_pts,trif_pts,trif_index,frac,             &
                     frac_agric,lit_c,dpm_ratio)

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
 frac(land_pts,ntype)                                             &
                          ! IN  Fractional cover of each
,frac_agric(land_pts)                                             &
                          ! IN  Agricultural (disturbed) frac
,lit_c(land_pts,npft)                                             &
                          ! IN  PFT carbon litter
!                                 !    (kg C/m2/360days).
,dpm_ratio(land_pts)      ! OUT DPM ratio of litter input

REAL                                                              &
 grass(land_pts)                                                  &
                          ! WORK  fractional grass cover
,crop(land_pts)                                                   &
                          ! WORK  fractional crop cover
,ltree(land_pts)                                                  &
                          ! WORK  litter input from trees
,lshrub(land_pts)                                                 &
                          ! WORK  litter input from shrubs
,lgrass_crop(land_pts)                                            &
                          ! WORK  litter input from grass+crop
,lgrass(land_pts)                                                 &
                          ! WORK  litter input from grass
,lcrop(land_pts)                                                  &
                          ! WORK  litter input from crops
,dpm(land_pts)            ! WORK  litter input to DPM

REAL                                                              &
 rt                                                               &
                            ! DPM:RPM ratio for trees
,rs                                                               &
                            ! DPM:RPM ratio for shrubs
,rg                                                               &
                            ! DPM:RPM ratio for grass
,rc                         ! DPM:RPM ratio for crops

PARAMETER(rt=0.25, rs=0.33, rg=0.67, rc=1.44)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('DPM_RPM',zhook_in,zhook_handle)

! calculate grass/crop fractions from C3, C4 cover, and FRAC_AGRIC
DO t=1,trif_pts
  l=trif_index(t)

  grass(l) = frac(l,3) + frac(l,4)
  IF (grass(l)  >   frac_agric(l)) THEN
    crop(l)  = frac_agric(l)
    grass(l) = grass(l) - frac_agric(l)
  ELSE
    crop(l)  = grass(l)
    grass(l) = 0.0
  END IF

END DO

! calculate different litter-type amounts
DO t=1,trif_pts
  l=trif_index(t)

  ltree(l)       = lit_c(l,1)*frac(l,1) + lit_c(l,2)*frac(l,2)
  lshrub(l)      = lit_c(l,5)*frac(l,5)
  lgrass_crop(l) = lit_c(l,3)*frac(l,3) + lit_c(l,4)*frac(l,4)
  lgrass(l)      = lgrass_crop(l)*grass(l) / (grass(l)+crop(l))
  lcrop(l)       = lgrass_crop(l)*crop(l) / (grass(l)+crop(l))

END DO

! calculate DPM litter input and hence ratio
DO t=1,trif_pts
  l=trif_index(t)

  dpm(l) = ltree(l)  * rt/(1+rt) + lshrub(l)*rs/(1+rs) +          &
           lgrass(l) * rg/(1+rg) +  lcrop(l)*rc/(1+rc)
  dpm_ratio(l) = dpm(l) /                                         &
               (ltree(l) + lshrub(l) + lgrass(l) + lcrop(l))

END DO


IF (lhook) CALL dr_hook('DPM_RPM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dpm_rpm
#endif
