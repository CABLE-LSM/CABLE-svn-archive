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
! Purpose: Calculates surface albedo, including from snow covered surface
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b (but was previously in cable_radiation)
!
!
! ==============================================================================

MODULE cable_albedo_module

  USE cable_data_module, ONLY : ialbedo_type, point2constants

  IMPLICIT NONE

  PUBLIC surface_albedo
  PRIVATE

  TYPE(ialbedo_type) :: C


CONTAINS


  SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)

    USE cable_common_module
    USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &
         canopy_type, met_type, radiation_type,      &
         soil_snow_type, mp, r_2, nrb
use cbl_rhoch_module, only : calc_rhoch
use cbl_snow_albedo_module, only : surface_albedosn

    TYPE (canopy_type),INTENT(IN)       :: canopy
    TYPE (met_type),INTENT(INOUT)       :: met
    TYPE (radiation_type),INTENT(INOUT) :: rad
    TYPE (soil_snow_type),INTENT(INOUT) :: ssnow

    TYPE (veg_parameter_type),INTENT(INOUT)  :: veg
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil

    REAL(r_2), DIMENSION(mp)  ::                                                &
         dummy2, & !
         dummy

    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: c1, rhoch

    LOGICAL, DIMENSION(mp)  :: mask ! select points for calculation

    INTEGER :: b    !rad. band 1=visible, 2=near-infrared, 3=long-wave

    ! END header

    CALL point2constants(C)

    IF (.NOT. ALLOCATED(c1)) &
         ALLOCATE( c1(mp,nrb), rhoch(mp,nrb) )


    CALL surface_albedosn(ssnow, veg, met, soil)

    rad%cexpkbm = 0.0
    rad%extkbm  = 0.0
    rad%rhocbm  = 0.0

    ! Initialise effective conopy beam reflectance:
    rad%reffbm = ssnow%albsoilsn
    rad%reffdf = ssnow%albsoilsn
    rad%albedo = ssnow%albsoilsn

    ! Define vegetation mask:
    mask = canopy%vlaiw > C%LAI_THRESH .AND.                                    &
         ( met%fsd(:,1) + met%fsd(:,2) ) > C%RAD_THRESH

call calc_rhoch( c1,rhoch, mp, nrb, veg%taul, veg%refl )

    ! Update extinction coefficients and fractional transmittance for
    ! leaf transmittance and reflection (ie. NOT black leaves):
    !---1 = visible, 2 = nir radiaition
    DO b = 1, 2

       rad%extkdm(:,b) = rad%extkd * c1(:,b)

       !--Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)

       !---Calculate effective diffuse reflectance (fraction):
       WHERE( canopy%vlaiw > 1e-2 )                                             &
            rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssnow%albsoilsn(:,b)             &
            - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2

       !---where vegetated and sunlit
       WHERE (mask)

          rad%extkbm(:,b) = rad%extkb * c1(:,b)

          ! Canopy reflection (6.21) beam:
          rad%rhocbm(:,b) = 2. * rad%extkb / ( rad%extkb + rad%extkd )          &
               * rhoch(:,b)

          ! Canopy beam transmittance (fraction):
          dummy2 = MIN(rad%extkbm(:,b)*canopy%vlaiw, 20.)
          dummy  = EXP(-dummy2)
          rad%cexpkbm(:,b) = REAL(dummy)

          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rad%rhocbm(:,b) + (ssnow%albsoilsn(:,b)             &
               - rad%rhocbm(:,b))*rad%cexpkbm(:,b)**2

       END WHERE

       ! Define albedo:
       WHERE( canopy%vlaiw> C%LAI_THRESH )                                      &
            rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*rad%reffdf(:,b) +           &
            rad%fbeam(:,b) * rad%reffbm(:,b)

    END DO


  END SUBROUTINE surface_albedo

END MODULE cable_albedo_module
