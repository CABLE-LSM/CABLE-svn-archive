MODULE cbl_albedo_mod

  IMPLICIT NONE

  PUBLIC surface_albedo
  PRIVATE

CONTAINS

SUBROUTINE surface_albedo(ssnow, veg, met, rad, soil, canopy)

!subrs called
USE cbl_rhoch_module, ONLY : calc_rhoch
USE cbl_snow_albedo_module, ONLY : surface_albedosn

USE cable_other_constants_mod, ONLY :  CLAI_thresh => LAI_thresh
USE cable_other_constants_mod, ONLY :  CRad_thresh => Rad_thresh
    USE cable_common_module
    USE cable_def_types_mod, ONLY : veg_parameter_type, soil_parameter_type,    &
         canopy_type, met_type, radiation_type,      &
         soil_snow_type, mp, r_2, nrb

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
    mask = canopy%vlaiw > CLAI_THRESH .AND.                                    &
         ( met%fsd(:,1) + met%fsd(:,2) ) > CRAD_THRESH

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
       WHERE( canopy%vlaiw> CLAI_THRESH )                                      &
            rad%albedo(:,b) = ( 1. - rad%fbeam(:,b) )*rad%reffdf(:,b) +           &
            rad%fbeam(:,b) * rad%reffbm(:,b)

    END DO


  END SUBROUTINE surface_albedo


END MODULE cbl_albedo_mod
