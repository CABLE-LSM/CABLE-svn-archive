MODULE cbl_LAI_eff_mod

IMPLICIT NONE

PUBLIC LAI_eff

CONTAINS

!Computes Effective LAI of exposed canopy given effect of snow present
!variable formerly known as canopy%vlaiw
SUBROUTINE LAI_eff( mp, lai_pft, Hgt_PFT, HgtAboveSnow,                        &
                    reducedLAIdue2snow )

  !re-decl input args
INTEGER  :: mp
REAL :: lai_pft(mp)
REAL :: Hgt_PFT(mp)
REAL :: HgtAboveSnow(mp)
REAL :: reducedLAIdue2snow(mp)

!local_vars:
REAL :: Hgt_eff(mp)
REAL :: FracOfCanopyAboveSnow(mp)

!Fraction Of Canopy Above Snow
FracOfCanopyAboveSnow = HgtAboveSnow/ MAX( 0.01, Hgt_PFT)

! LAI decreases due to snow:
reducedLAIdue2snow = lai_pft * FracOfCanopyAboveSnow

END SUBROUTINE LAI_eff

END MODULE cbl_LAI_eff_mod
