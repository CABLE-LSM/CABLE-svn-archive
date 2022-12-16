!******************************************************************************
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) can be found
! at https://github.com/CABLE-LSM/CABLE/blob/main/
!
!******************************************************************************

MODULE hruff_eff_LAI_mod_cbl

!-----------------------------------------------------------------------------
! Description:
!   Computes the height above the surface given that there is snow present
!   and the effective LAI of the canopy given the effect of snow
!
! This MODULE is USEd in:
!     cable_land_albedo_mod_cbl.F90 (JULES)
!
! This MODULE contains 2 public Subroutines:
!     HgtAboveSnow,
!     LAI_eff
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC :: HgtAboveSnow
PUBLIC :: LAI_eff

CONTAINS

!variable formerly known as rough%hruff
SUBROUTINE HgtAboveSnow( HeightAboveSnow, mp, z0surf_min, HGT_pft,             &
                         SnowDepth, SnowDensity )
! Description:
!   Computes height above surface given that therre is snow present

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mp         ! CABLE VECTOR LENGTH

REAL, INTENT(OUT) :: HeightAboveSnow(mp) ! result to return

REAL, INTENT(IN) :: z0surf_min      ! min. surface roughness
REAL, INTENT(IN) :: HGT_pft(mp)     ! canopy height
REAL, INTENT(IN) :: SnowDepth(mp)   ! snow depth (liquid water )
REAL, INTENT(IN) :: SnowDensity(mp) ! snow density

!local_vars:
REAL, PARAMETER  :: fmin = 10.0 ! [meters]?
REAL, PARAMETER  :: SnowDensity_min = 100.0 ! min. snow density

REAL :: SnowDensity_eff(mp)         ! Effective snow density range restricted
REAL :: HgtAboveSnow_min(mp)        ! min. canopy height
REAL :: HgtAboveSnow_comp(mp)       ! computed canopy height above snow

! restrict Effective snow density to be .GE. "a" minimum
SnowDensity_eff= MAX( SnowDensity_min, SnowDensity )

! min. allowed canopy height (fixed @ 10* min. surface roughness)
HgtAboveSnow_min =  fmin * z0surf_min

!Canopy Hgt above snow given computed snow depth & PFT height
HgtAboveSnow_comp =   HGT_pft - ( 1.2 * SnowDepth / SnowDensity_eff )

! Finally Set Effective canopy height above snow level
HeightAboveSnow = MAX( HgtAboveSnow_min,  HgtAboveSnow_comp )

RETURN

END SUBROUTINE HgtAboveSnow

!variable formerly known as canopy%vlaiw
SUBROUTINE LAI_eff( mp, lai_pft, Hgt_PFT, HgtAboveSnow,                        &
                    reducedLAIdue2snow )
! Description:
!   Computes Effective LAI of exposed canopy given effect of snow present

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mp           ! CABLE VECTOR LENGTH
! return result - considered LAI seen given snow coverage
REAL, INTENT(OUT) :: reducedLAIdue2snow(mp)
REAL, INTENT(IN) :: lai_pft(mp)       ! LAI
REAL, INTENT(IN) :: HGT_pft(mp)       ! canopy height
REAL, INTENT(IN) :: HgtAboveSnow(mp)  ! computed canopy height above snow

!local_vars:
REAL :: Hgt_eff(mp)
REAL :: FracOfCanopyAboveSnow(mp)

!Fraction Of Canopy Above Snow
FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)

! LAI decreases due to snow:
reducedLAIdue2snow = lai_pft * FracOfCanopyAboveSnow

RETURN
END SUBROUTINE LAI_eff

END MODULE hruff_eff_LAI_mod_cbl
