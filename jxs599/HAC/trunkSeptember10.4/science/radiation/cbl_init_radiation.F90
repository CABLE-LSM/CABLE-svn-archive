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

  USE cable_data_module, ONLY : irad_type, point2constants

  IMPLICIT NONE

  PUBLIC init_radiation
  PRIVATE

  TYPE ( irad_type ) :: C


CONTAINS

  SUBROUTINE init_radiation( met, rad, veg, canopy )

    ! Alternate version of init_radiation that uses only the
    ! zenith angle instead of the fluxes. This means it can be called
    ! before the cable albedo calculation.
    USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
         veg_parameter_type, nrb, mp
    USE cable_common_module

    TYPE (radiation_type), INTENT(INOUT) :: rad
    TYPE (met_type),       INTENT(INOUT) :: met

    TYPE (canopy_type),    INTENT(IN)    :: canopy

    TYPE (veg_parameter_type), INTENT(INOUT) :: veg

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


    CALL point2constants( C )

    IF(.NOT. ALLOCATED(c1) ) ALLOCATE( c1(mp,nrb), rhoch(mp,nrb) )

    cos3 = COS(C%PI180 * (/ 15.0, 45.0, 75.0 /))

    ! See Sellers 1985, eq.13 (leaf angle parameters):
    WHERE (canopy%vlaiw > C%LAI_THRESH)
       xphi1 = 0.5 - veg%xfang * (0.633 + 0.33 * veg%xfang)
       xphi2 = 0.877 * (1.0 - 2.0 * xphi1)
    END WHERE

    ! 2 dimensional LAI
    xvlai2 = SPREAD(canopy%vlaiw, 2, 3)

    ! Extinction coefficient for beam radiation and black leaves;
    ! eq. B6, Wang and Leuning, 1998
    WHERE (xvlai2 > C%LAI_THRESH) ! vegetated
       xk = SPREAD(xphi1, 2, 3) / SPREAD(cos3, 1, mp) + SPREAD(xphi2, 2, 3)
    ELSEWHERE ! i.e. bare soil
       xk = 0.0
    END WHERE

    WHERE (canopy%vlaiw > C%LAI_THRESH ) ! vegetated

       ! Extinction coefficient for diffuse radiation for black leaves:
       rad%extkd = -LOG( SUM(                                                   &
            SPREAD( C%GAUSS_W, 1, mp ) * EXP( -xk * xvlai2 ), 2) )       &
            / canopy%vlaiw

    ELSEWHERE ! i.e. bare soil
       rad%extkd = 0.7
    END WHERE


    CALL calc_rhoch( veg, c1, rhoch )

    ! Canopy REFLection of diffuse radiation for black leaves:
    DO ictr=1,nrb

       rad%rhocdf(:,ictr) = rhoch(:,ictr) *  2. *                                &
            ( C%GAUSS_W(1) * xk(:,1) / ( xk(:,1) + rad%extkd(:) )&
            + C%GAUSS_W(2) * xk(:,2) / ( xk(:,2) + rad%extkd(:) )&
            + C%GAUSS_W(3) * xk(:,3) / ( xk(:,3) + rad%extkd(:) ) )

    ENDDO

    IF( .NOT. cable_runtime%um) THEN

       ! Define beam fraction, fbeam:
       rad%fbeam(:,1) = spitter(met%doy, met%coszen, met%fsd(:,1))
       rad%fbeam(:,2) = spitter(met%doy, met%coszen, met%fsd(:,2))

       ! coszen is set during met data read in.

       WHERE (met%coszen <1.0e-2)
          rad%fbeam(:,1) = 0.0
          rad%fbeam(:,2) = 0.0
       END WHERE

    ENDIF

    ! In gridcells where vegetation exists....

    !!vh !! include RAD_THRESH in condition
    WHERE (canopy%vlaiw > C%LAI_THRESH .AND. met%coszen > 1.e-6 )

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
