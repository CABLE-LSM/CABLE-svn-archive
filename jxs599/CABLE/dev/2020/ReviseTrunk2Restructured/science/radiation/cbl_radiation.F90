!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Computes radiation absorbed by canopy and soil surface
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

MODULE cbl_radiation_module

  IMPLICIT NONE

  PUBLIC radiation
  PRIVATE

CONTAINS

!newSUBROUTINE radiation( ssnow, veg, air, met, rad, canopy, sunlit_veg_mask,&
!new  !constants
!new  clai_thresh, Csboltz, Cemsoil, Cemleaf, Ccapp &
!new)
!new
!new    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
!new         veg_parameter_type, soil_snow_type,         &
!new         air_type, mp, mf, r_2
!new
!newIMPLICIT NONE
!newlogical :: sunlit_veg_mask(mp)
!new!constants
!newreal :: CLAI_thresh
!newreal :: CSboltz
!newreal :: Cemsoil
!newreal :: Cemleaf
!newreal :: Ccapp
!new
!new    TYPE (canopy_type),   INTENT(IN) :: canopy
!new    TYPE (air_type),      INTENT(IN) :: air
!new    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow
!new    TYPE (met_type),      INTENT(INOUT) :: met
!new    TYPE (radiation_type),INTENT(INOUT) :: rad
!new
!new    TYPE (veg_parameter_type), INTENT(IN) :: veg
!new
!new    REAL, DIMENSION(mp) ::                                                      &
!new         cf1, &      ! (1.0 - rad%transb * cexpkdm) / (extkb + extkdm(:,b))
!new         cf3, &      ! (1.0 - rad%transb * cexpkbm) / (extkb + extkbm(:,b))
!new         cf2n, &     ! exp(-extkn * vlai) (nitrogen)
!new         emair, &    ! air emissivity
!new         flpwb, &    ! black-body long-wave radiation
!new         flwv, &     ! vegetation long-wave radiation (isothermal)
!new         dummy, dummy2
!new
!new    INTEGER :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
!new
!new    INTEGER, SAVE :: call_number =0
!new
!new    call_number = call_number + 1
!new
!new    ! Relative leaf nitrogen concentration within canopy:
!new    cf2n = EXP(-veg%extkn * canopy%vlaiw)
!new
!new    rad%transd = 1.0
!new
!new    WHERE (canopy%vlaiw > cLAI_thresh )    ! where vegetation exists....
!new
!new       ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
!new       ! leaf SW transmittance and REFLectance);
!new       ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
!new       rad%transd = EXP(-rad%extkd * canopy%vlaiw)
!new
!new    END WHERE
!new
!new    ! Define fraction of SW beam tranmitted through canopy:
!new    !C!jhan: check rel. b/n extkb, extkbm,transb,cexpkbm def. cable_albedo, qsabbs
!new    !! vh_js !!
!new    dummy2 = MIN(rad%extkb * canopy%vlaiw,30.) ! vh version to avoid floating underflow !
!new    dummy = EXP(-dummy2)
!new    ! dummy2 = -rad%extkb * canopy%vlaiw
!new    ! dummy = EXP(dummy2)
!new    rad%transb = REAL(dummy)
!new
!new    ! Define longwave from vegetation:
!new    flpwb = CSboltz * (met%tvrad) ** 4
!new    flwv = Cemleaf * flpwb
!new
!new    rad%flws = CSboltz*Cemsoil* ssnow%tss **4
!new
!new    ! Define air emissivity:
!new    emair = met%fld / flpwb
!new
!new    rad%gradis = 0.0 ! initialise radiative conductance
!new    rad%qcan = 0.0   ! initialise radiation absorbed by canopy
!new
!new    WHERE (canopy%vlaiw > CLAI_thresh )
!new
!new       ! Define radiative conductance (Leuning et al, 1995), eq. D7:
!new       rad%gradis(:,1) = ( 4.0 * Cemleaf / (Ccapp * air%rho) ) * flpwb        &
!new            / (met%tvrad) * rad%extkd                              &
!new            * ( ( 1.0 - rad%transb * rad%transd ) /                &
!new            ( rad%extkb + rad%extkd )                              &
!new            + ( rad%transd - rad%transb ) /                        &
!new            ( rad%extkb - rad%extkd ) )
!new
!new       rad%gradis(:,2) = ( 8.0 * Cemleaf / ( Ccapp * air%rho ) ) *            &
!new            flpwb / met%tvrad * rad%extkd *                        &
!new            ( 1.0 - rad%transd ) / rad%extkd - rad%gradis(:,1)
!new
!new       ! Longwave radiation absorbed by sunlit canopy fraction:
!new       rad%qcan(:,1,3) = (rad%flws - flwv ) * rad%extkd *                       &
!new            ( rad%transd - rad%transb ) / ( rad%extkb - rad%extkd )&
!new            + ( emair- Cemleaf ) * rad%extkd * flpwb *            &
!new            ( 1.0 - rad%transd * rad%transb )                      &
!new            / ( rad%extkb + rad%extkd )
!new
!new       ! Longwave radiation absorbed by shaded canopy fraction:
!new       rad%qcan(:,2,3) = ( 1.0 - rad%transd ) *                                 &
!new            ( rad%flws + met%fld - 2.0 * flwv ) - rad%qcan(:,1,3)
!new
!new    END WHERE
!new
!new    ! Convert radiative conductance from m/s to mol/m2/s:
!new    rad%gradis=SPREAD(air%cmolar, 2, mf)*rad%gradis
!new    rad%gradis = MAX(1.0e-3_r_2,rad%gradis)
!new
!new    ! Update extinction coefficients and fractional transmittance for
!new    ! leaf transmittance and REFLection (ie. NOT black leaves):
!new    ! Define qcan for short wave (par, nir) for sunlit leaf:
!new    ! UM recieves met%fsd(:,b) forcing. assumed for offline that USED met%fsd(:,b) = 1/2* INPUT met%fsd
!new    DO b = 1, 2 ! 1 = visible, 2 = nir radiaition
!new
!new       WHERE (sunlit_veg_mask) ! i.e. vegetation and sunlight are present
!new
!new          cf1 = ( 1.0 - rad%transb * rad%cexpkdm(:,b) ) /                       &
!new               ( rad%extkb + rad%extkdm(:,b) )
!new          cf3 = (1.0 - rad%transb * rad%cexpkbm(:,b)) /                         &
!new               ( rad%extkb + rad%extkbm(:,b) )
!new
!new          ! scale to real sunlit flux
!new          rad%qcan(:,1,b) = met%fsd(:,b) * (                                    &
!new               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
!new               * rad%extkdm(:,b) * cf1                             &
!new               + rad%fbeam(:,b) * ( 1.0-rad%reffbm(:,b) ) *        &
!new               rad%extkbm(:,b) * cf3                               &
!new               + rad%fbeam(:,b) * ( 1.0 - veg%taul(:,b)            &
!new               - veg%refl(:,b) ) * rad%extkb                       &
!new               * ( ( 1-rad%transb ) / rad%extkb - ( 1 -            &
!new               rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )
!new
!new          ! Define qcan for short wave (par, nir) for shaded leaf:
!new          rad%qcan(:,2,b) = met%fsd(:,b) * (                                    &
!new               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
!new               * rad%extkdm(:,b) *                                 &
!new               ( ( 1.0 - rad%cexpkdm(:,b) ) / rad%extkdm(:,b)      &
!new               - cf1 ) + rad%fbeam(:,b) * ( 1. - rad%reffbm(:,b) ) &
!new               * rad%extkbm(:,b) * ( ( 1.0 - rad%cexpkbm(:,b) ) /  &
!new               rad%extkbm(:,b) - cf3 ) - rad%fbeam(:,b) *          &
!new               ( 1.0 - veg%taul(:,b) -veg%refl(:,b)) * rad%extkb   &
!new               * ( ( 1 - rad%transb ) / rad%extkb -                &
!new               ( 1 - rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )
!new
!new       END WHERE
!new
!new    END DO
!new
!new    rad%qssabs = 0.
!new
!new    WHERE (sunlit_veg_mask) ! i.e. vegetation and sunlight are present
!new
!new       ! Calculate shortwave radiation absorbed by soil:
!new       ! (av. of transmitted NIR and PAR through canopy)*SWdown
!new       rad%qssabs = met%fsd(:,1) * (                                            &
!new            rad%fbeam(:,1) * ( 1. - rad%reffbm(:,1) ) *                 &
!new            EXP( -MIN(rad%extkbm(:,1) * canopy%vlaiw,20.) ) +           &
!new            ( 1. - rad%fbeam(:,1) ) * ( 1. - rad%reffdf(:,1) ) *        &
!new            EXP( -MIN(rad%extkdm(:,1) * canopy%vlaiw,20.) ) )           &
!new            + met%fsd(:,2) * ( rad%fbeam(:,2) * ( 1. - rad%reffbm(:,2) )&
!new            * rad%cexpkbm(:,2) + ( 1. - rad%fbeam(:,2) ) *              &
!new            ( 1. - rad%reffdf(:,2) ) * rad%cexpkdm(:,2) )
!new
!new       ! Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C:
!new       rad%scalex(:,1) = ( 1.0 - rad%transb * cf2n ) / ( rad%extkb + veg%extkn )
!new
!new       ! LAI of big leaf, sunlit, shaded, respectively:
!new       rad%fvlai(:,1) = ( 1.0 - rad%transb ) / rad%extkb
!new       rad%fvlai(:,2) = canopy%vlaiw - rad%fvlai(:,1)
!new
!new    ELSEWHERE ! i.e. either vegetation or sunlight are NOT present
!new
!new       ! Shortwave absorbed by soil/snow surface:
!new       rad%qssabs = ( 1.0 - ssnow%albsoilsn(:,1) ) * met%fsd(:,1) +             &
!new            ( 1.0 - ssnow%albsoilsn(:,2) ) * met%fsd(:,2)
!new
!new       rad%scalex(:,1) = 0.0
!new       rad%fvlai(:,1) = 0.0
!new       rad%fvlai(:,2) = canopy%vlaiw
!new
!new    END WHERE
!new
!new    rad%scalex(:,2) = (1.0 - cf2n) / veg%extkn - rad%scalex(:,1)
!new
!new    ! Total energy absorbed by canopy:
!new    rad%rniso = SUM(rad%qcan, 3)
!new
!new  END SUBROUTINE radiation

  SUBROUTINE radiation( ssnow, veg, air, met, rad, canopy )

    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type, soil_snow_type,         &
         air_type, mp, mf, r_2

USE cable_other_constants_mod, ONLY : CLAI_thresh => LAI_thresh
USE cable_other_constants_mod, ONLY : CRad_thresh => Rad_thresh
USE cable_phys_constants_mod, ONLY :  CEMsoil => EMsoil 
USE cable_phys_constants_mod, ONLY :  CEMleaf => EMleaf 
USE cable_phys_constants_mod, ONLY :  CSboltz => Sboltz
USE cable_phys_constants_mod, ONLY :  CCapp   => Capp   
    TYPE (canopy_type),   INTENT(IN) :: canopy
    TYPE (air_type),      INTENT(IN) :: air
    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow
    TYPE (met_type),      INTENT(INOUT) :: met
    TYPE (radiation_type),INTENT(INOUT) :: rad

    TYPE (veg_parameter_type), INTENT(IN) :: veg

    REAL, DIMENSION(mp) ::                                                      &
         cf1, &      ! (1.0 - rad%transb * cexpkdm) / (extkb + extkdm(:,b))
         cf3, &      ! (1.0 - rad%transb * cexpkbm) / (extkb + extkbm(:,b))
         cf2n, &     ! exp(-extkn * vlai) (nitrogen)
         emair, &    ! air emissivity
         flpwb, &    ! black-body long-wave radiation
         flwv, &     ! vegetation long-wave radiation (isothermal)
         dummy, dummy2

    LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation

    INTEGER :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave

    INTEGER, SAVE :: call_number =0

    call_number = call_number + 1

    ! Define vegetation mask:
    mask = canopy%vlaiw > CLAI_THRESH .AND.                                    &
         ( met%fsd(:,1)+met%fsd(:,2) ) > CRAD_THRESH

    ! Relative leaf nitrogen concentration within canopy:
    cf2n = EXP(-veg%extkn * canopy%vlaiw)

    rad%transd = 1.0

    WHERE (canopy%vlaiw > CLAI_THRESH )    ! where vegetation exists....

       ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance);
       ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
       rad%transd = EXP(-rad%extkd * canopy%vlaiw)

    END WHERE

    ! Define fraction of SW beam tranmitted through canopy:
    !C!jhan: check rel. b/n extkb, extkbm,transb,cexpkbm def. cable_albedo, qsabbs
    !! vh_js !!
    dummy2 = MIN(rad%extkb * canopy%vlaiw,30.) ! vh version to avoid floating underflow !
    dummy = EXP(-dummy2)
    ! dummy2 = -rad%extkb * canopy%vlaiw
    ! dummy = EXP(dummy2)
    rad%transb = REAL(dummy)

    ! Define longwave from vegetation:
    flpwb = Csboltz * (met%tvrad) ** 4
    flwv = CEMLEAF * flpwb

    rad%flws = Csboltz*CEMSOIL* ssnow%tss **4

    ! Define air emissivity:
    emair = met%fld / flpwb

    rad%gradis = 0.0 ! initialise radiative conductance
    rad%qcan = 0.0   ! initialise radiation absorbed by canopy

    WHERE (canopy%vlaiw > CLAI_THRESH )

       ! Define radiative conductance (Leuning et al, 1995), eq. D7:
       rad%gradis(:,1) = ( 4.0 * CEMLEAF / (CCAPP * air%rho) ) * flpwb        &
            / (met%tvrad) * rad%extkd                              &
            * ( ( 1.0 - rad%transb * rad%transd ) /                &
            ( rad%extkb + rad%extkd )                              &
            + ( rad%transd - rad%transb ) /                        &
            ( rad%extkb - rad%extkd ) )

       rad%gradis(:,2) = ( 8.0 * CEMLEAF / ( CCAPP * air%rho ) ) *            &
            flpwb / met%tvrad * rad%extkd *                        &
            ( 1.0 - rad%transd ) / rad%extkd - rad%gradis(:,1)

       ! Longwave radiation absorbed by sunlit canopy fraction:
       rad%qcan(:,1,3) = (rad%flws - flwv ) * rad%extkd *                       &
            ( rad%transd - rad%transb ) / ( rad%extkb - rad%extkd )&
            + ( emair- CEMLEAF ) * rad%extkd * flpwb *            &
            ( 1.0 - rad%transd * rad%transb )                      &
            / ( rad%extkb + rad%extkd )

       ! Longwave radiation absorbed by shaded canopy fraction:
       rad%qcan(:,2,3) = ( 1.0 - rad%transd ) *                                 &
            ( rad%flws + met%fld - 2.0 * flwv ) - rad%qcan(:,1,3)

    END WHERE

    ! Convert radiative conductance from m/s to mol/m2/s:
    rad%gradis=SPREAD(air%cmolar, 2, mf)*rad%gradis
    rad%gradis = MAX(1.0e-3_r_2,rad%gradis)

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and REFLection (ie. NOT black leaves):
    ! Define qcan for short wave (par, nir) for sunlit leaf:
    ! UM recieves met%fsd(:,b) forcing. assumed for offline that USED met%fsd(:,b) = 1/2* INPUT met%fsd
    DO b = 1, 2 ! 1 = visible, 2 = nir radiaition

       WHERE (mask) ! i.e. vegetation and sunlight are present

          cf1 = ( 1.0 - rad%transb * rad%cexpkdm(:,b) ) /                       &
               ( rad%extkb + rad%extkdm(:,b) )
          cf3 = (1.0 - rad%transb * rad%cexpkbm(:,b)) /                         &
               ( rad%extkb + rad%extkbm(:,b) )

          ! scale to real sunlit flux
          rad%qcan(:,1,b) = met%fsd(:,b) * (                                    &
               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
               * rad%extkdm(:,b) * cf1                             &
               + rad%fbeam(:,b) * ( 1.0-rad%reffbm(:,b) ) *        &
               rad%extkbm(:,b) * cf3                               &
               + rad%fbeam(:,b) * ( 1.0 - veg%taul(:,b)            &
               - veg%refl(:,b) ) * rad%extkb                       &
               * ( ( 1-rad%transb ) / rad%extkb - ( 1 -            &
               rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )

          ! Define qcan for short wave (par, nir) for shaded leaf:
          rad%qcan(:,2,b) = met%fsd(:,b) * (                                    &
               ( 1.0 - rad%fbeam(:,b) ) * ( 1.0 - rad%reffdf(:,b) )&
               * rad%extkdm(:,b) *                                 &
               ( ( 1.0 - rad%cexpkdm(:,b) ) / rad%extkdm(:,b)      &
               - cf1 ) + rad%fbeam(:,b) * ( 1. - rad%reffbm(:,b) ) &
               * rad%extkbm(:,b) * ( ( 1.0 - rad%cexpkbm(:,b) ) /  &
               rad%extkbm(:,b) - cf3 ) - rad%fbeam(:,b) *          &
               ( 1.0 - veg%taul(:,b) -veg%refl(:,b)) * rad%extkb   &
               * ( ( 1 - rad%transb ) / rad%extkb -                &
               ( 1 - rad%transb**2 ) / ( rad%extkb + rad%extkb ) ) )

       END WHERE

    END DO

    rad%qssabs = 0.

    WHERE (mask) ! i.e. vegetation and sunlight are present

       ! Calculate shortwave radiation absorbed by soil:
       ! (av. of transmitted NIR and PAR through canopy)*SWdown
       rad%qssabs = met%fsd(:,1) * (                                            &
            rad%fbeam(:,1) * ( 1. - rad%reffbm(:,1) ) *                 &
            EXP( -MIN(rad%extkbm(:,1) * canopy%vlaiw,20.) ) +           &
            ( 1. - rad%fbeam(:,1) ) * ( 1. - rad%reffdf(:,1) ) *        &
            EXP( -MIN(rad%extkdm(:,1) * canopy%vlaiw,20.) ) )           &
            + met%fsd(:,2) * ( rad%fbeam(:,2) * ( 1. - rad%reffbm(:,2) )&
            * rad%cexpkbm(:,2) + ( 1. - rad%fbeam(:,2) ) *              &
            ( 1. - rad%reffdf(:,2) ) * rad%cexpkdm(:,2) )

       ! Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C:
       rad%scalex(:,1) = ( 1.0 - rad%transb * cf2n ) / ( rad%extkb + veg%extkn )

       ! LAI of big leaf, sunlit, shaded, respectively:
       rad%fvlai(:,1) = ( 1.0 - rad%transb ) / rad%extkb
       rad%fvlai(:,2) = canopy%vlaiw - rad%fvlai(:,1)

    ELSEWHERE ! i.e. either vegetation or sunlight are NOT present

       ! Shortwave absorbed by soil/snow surface:
       rad%qssabs = ( 1.0 - ssnow%albsoilsn(:,1) ) * met%fsd(:,1) +             &
            ( 1.0 - ssnow%albsoilsn(:,2) ) * met%fsd(:,2)

       rad%scalex(:,1) = 0.0
       rad%fvlai(:,1) = 0.0
       rad%fvlai(:,2) = canopy%vlaiw

    END WHERE

    rad%scalex(:,2) = (1.0 - cf2n) / veg%extkn - rad%scalex(:,1)

    ! Total energy absorbed by canopy:
    rad%rniso = SUM(rad%qcan, 3)

  END SUBROUTINE radiation


END MODULE cbl_radiation_module
