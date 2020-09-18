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

MODULE cbl_init_radiation_module

  IMPLICIT NONE

  PUBLIC init_radiation
  PRIVATE

!FUDGED local pars -masks tuned at other times - review conssitency!!
real :: Ccoszen_tols_huge  ! 1e-4 * threshold cosine of sun's zenith angle, below which considered SUNLIT
real :: Ccoszen_tols_tiny  ! 1e-4 * threshold cosine of sun's zenith angle, below which considered SUNLIT

CONTAINS

  SUBROUTINE init_radiation( met, rad, veg, canopy, veg_mask, sunlit_veg_mask)

    ! Alternate version of init_radiation that uses only the
    ! zenith angle instead of the fluxes. This means it can be called
    ! before the cable albedo calculation.
    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type, nrb, mp
    USE cable_common_module

USE cbl_spitter_module, ONLY : Spitter
USE cbl_rhoch_module, ONLY : calc_rhoch
USE cable_math_constants_mod, ONLY : CPI => PI
USE cable_math_constants_mod, ONLY : CPI180 => PI180
USE cable_other_constants_mod, ONLY : CLAI_thresh => LAI_thresh
USE cable_other_constants_mod, ONLY : CGauss_W => Gauss_W

    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (met_type),       INTENT(INOUT) :: met

    TYPE (canopy_type),    INTENT(IN)    :: canopy

    TYPE (veg_parameter_type), INTENT(INOUT) :: veg

logical :: veg_mask(mp) 
logical :: sunlit_veg_mask(mp) 

    REAL, DIMENSION(nrb) ::                                                     &
         cos3       ! cos(15 45 75 degrees)
    REAL, DIMENSION(mp,nrb) ::                                                  &
         xvlai2,  & ! 2D vlai
         xk         ! extinct. coef.for beam rad. and black leaves

    REAL, DIMENSION(mp) ::                                                      &
         xphi1,   & ! leaf angle parmameter 1
         xphi2      ! leaf angle parmameter 2

    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE ::                                  &
                                ! subr to calc these curr. appears twice. fix this
         c1,      & !
         rhoch



    INTEGER :: ictr


    !CALL point2constants( C )

    IF(.NOT. ALLOCATED(c1) ) ALLOCATE( c1(mp,nrb), rhoch(mp,nrb) )

    cos3 = COS(CPI180 * (/ 15.0, 45.0, 75.0 /))

    ! See Sellers 1985, eq.13 (leaf angle parameters):
    WHERE ( veg_mask )
       xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
       xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
    END WHERE

    ! 2 dimensional LAI
    xvlai2 = SPREAD(canopy%vlaiw, 2, 3)

    ! Extinction coefficient for beam radiation and black leaves;
    ! eq. B6, Wang and Leuning, 1998
    WHERE (xvlai2 > CLAI_THRESH) ! vegetated
       xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
    ELSEWHERE ! i.e. bare soil
       xk = 0.0
    END WHERE

    WHERE ( veg_mask )

       ! Extinction coefficient for diffuse radiation for black leaves:
       rad%extkd = -LOG( SUM(                                                   &
            SPREAD( CGAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
            / canopy%vlaiw

    ELSEWHERE ! i.e. bare soil
       rad%extkd = 0.7
    END WHERE


CALL calc_rhoch( c1,rhoch, mp, nrb, veg%taul, veg%refl )

    IF( .NOT. cable_runtime%um) THEN

       ! Define beam fraction, fbeam:
       rad%fbeam(:,1) = spitter(mp, cpi, met%doy, met%coszen, met%fsd(:,1))
       rad%fbeam(:,2) = spitter(mp, cpi, met%doy, met%coszen, met%fsd(:,2))

       ! coszen is set during met data read in.

       WHERE (met%coszen <1.0e-2)
          rad%fbeam(:,1) = 0.0
          rad%fbeam(:,2) = 0.0
       END WHERE

    ENDIF

    ! In gridcells where vegetation exists....
    WHERE ( sunlit_veg_mask )

       ! SW beam extinction coefficient ("black" leaves, extinction neglects
       ! leaf SW transmittance and REFLectance):
       rad%extkb = xphi1 / met%coszen + xphi2

    ELSEWHERE ! i.e. bare soil
       rad%extkb = 0.5
    END WHERE

    WHERE ( ABS(rad%extkb - rad%extkd)  < 0.001 )
       rad%extkb = rad%extkd + 0.001
    END WHERE

    WHERE( met%coszen < 1.e-6 )
       ! higher value precludes sunlit leaves at night. affects
       ! nighttime evaporation - Ticket #90
       rad%extkb=1.0e5
    END WHERE

  END SUBROUTINE init_radiation

END MODULE cbl_init_radiation_module
