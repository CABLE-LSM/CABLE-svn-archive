!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Computes radiation absorbed by canopy and soil surface
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

MODULE cable_radiation_module

   USE cable_data_module, ONLY : irad_type, point2constants
 
   IMPLICIT NONE

   PUBLIC init_radiation, radiation
   PRIVATE

  TYPE ( irad_type ) :: C 


CONTAINS

SUBROUTINE init_radiation(met,rad,veg, canopy, c1, rhoch, xk) 

USE cbl_rhoch_module,   ONLY : calc_rhoch
USE cbl_spitter_module, ONLY : Spitter
USE cable_math_constants_mod,  ONLY: cpi => pi
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
   
    REAL :: c1(mp,nrb), rhoch(mp,nrb)

   LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation

   INTEGER :: ictr
  
   
   CALL point2constants( C ) 
   
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

   mask = canopy%vlaiw > C%LAI_THRESH  .AND.                                   &
          ( met%fsd(:,1) + met%fsd(:,2) ) > C%RAD_THRESH

CALL calc_rhoch( c1,rhoch, mp, nrb, veg%taul, veg%refl )

   ! Canopy REFLection of diffuse radiation for black leaves:
   DO ictr=1,nrb
     
     rad%rhocdf(:,ictr) = rhoch(:,ictr) *                                      &
                          ( C%GAUSS_W(1) * xk(:,1) / ( xk(:,1) + rad%extkd(:) )&
                          + C%GAUSS_W(2) * xk(:,2) / ( xk(:,2) + rad%extkd(:) )&
                          + C%GAUSS_W(3) * xk(:,3) / ( xk(:,3) + rad%extkd(:) ) )

   ENDDO
   
   IF( .NOT. cable_runtime%um) THEN
   
      ! Define beam fraction, fbeam:
      rad%fbeam(:,1) = spitter(mp, cpi, INT(met%doy), met%coszen, met%fsd(:,1))
      rad%fbeam(:,2) = spitter(mp, cpi, INT(met%doy), met%coszen, met%fsd(:,2))
      ! coszen is set during met data read in.
   
      WHERE (met%coszen <1.0e-2)
         rad%fbeam(:,1) = 0.0
         rad%fbeam(:,2) = 0.0
      END WHERE
   
   ENDIF
   
   ! In gridcells where vegetation exists....
   WHERE (canopy%vlaiw > C%LAI_THRESH)    
      
      ! SW beam extinction coefficient ("black" leaves, extinction neglects
      ! leaf SW transmittance and REFLectance):
      rad%extkb = xphi1 / met%coszen + xphi2
   
   ELSEWHERE ! i.e. bare soil
      rad%extkb = 0.5
   END WHERE
   
   WHERE ( abs(rad%extkb - rad%extkd)  < 0.001 )
      rad%extkb = rad%extkd + 0.001
   END WHERE
   
   WHERE(rad%fbeam(:,1) < C%RAD_THRESH )
      rad%extkb=30.0         ! keep cexpkbm within real*4 range (BP jul2010)
   END WHERE
!borrowed from cbl_iniit_radiation IN CABLE3 - bc there it is moved out oof Albedo and into iniit_radiation
!so here we are missing this definition   
rad%extkbm(:,:) = 0.0
rad%extkdm(:,:) = 0.0

call EffectiveExtinctCoeffs( rad%extkbm, rad%extkdm,               &
                             mp, nrb, mask,                        &
                             rad%extkb, rad%extkd, c1 )
   
END SUBROUTINE init_radiation

! ------------------------------------------------------------------------------

SUBROUTINE radiation( ssnow, veg, air, met, rad, canopy )
   
   USE cable_def_types_mod, ONLY : radiation_type, met_type, canopy_type,      &
                                   veg_parameter_type, soil_snow_type,         &
                                   air_type, mp, mf, r_2
                                       
   USE cable_common_module, only : cable_runtime, cable_user

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
      xx1,tssp    ! 
      
   REAL(r_2), DIMENSION(mp) ::                                                 &
      dummy, dummy2
   
   LOGICAL, DIMENSION(mp)    :: mask   ! select points for calculation
   
   INTEGER :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
   
   INTEGER, SAVE :: call_number =0
   
   REAL s1,s2,s3,step

   call_number = call_number + 1

   ! Define vegetation mask:
   mask = canopy%vlaiw > C%LAI_THRESH .AND.                                    &
          ( met%fsd(:,1)+met%fsd(:,2) ) > C%RAD_THRESH 

   ! Relative leaf nitrogen concentration within canopy:
   cf2n = EXP(-veg%extkn * canopy%vlaiw)
   
   rad%transd = 1.0
   
   WHERE (canopy%vlaiw > C%LAI_THRESH )    ! where vegetation exists....
            
      ! Diffuse SW transmission fraction ("black" leaves, extinction neglects
      ! leaf SW transmittance and REFLectance);
      ! from Monsi & Saeki 1953, quoted in eq. 18 of Sellers 1985:
      rad%transd = EXP(-rad%extkd * canopy%vlaiw)

   END WHERE

   ! Define fraction of SW beam tranmitted through canopy:
   dummy2 = -rad%extkb * canopy%vlaiw
   dummy = EXP(dummy2)
   rad%transb = REAL(dummy)

   ! Define longwave from vegetation:
   flpwb = C%sboltz * (met%tk) ** 4
   flwv = C%EMLEAF * flpwb

   rad%flws = C%sboltz*C%EMSOIL* ssnow%tss **4
   
   ! Define air emissivity:
   emair = met%fld / flpwb
    
   rad%gradis = 0.0 ! initialise radiative conductance
   rad%qcan = 0.0   ! initialise radiation absorbed by canopy
   
   WHERE (canopy%vlaiw > C%LAI_THRESH )

      ! Define radiative conductance (Leuning et al, 1995), eq. D7:
      rad%gradis(:,1) = ( 4.0 * C%EMLEAF / (C%CAPP * air%rho) ) * flpwb        &
                        / (met%tk) * rad%extkd                                 &
                        * ( ( 1.0 - rad%transb * rad%transd ) /                &
                        ( rad%extkb + rad%extkd )                              &
                        + ( rad%transd - rad%transb ) /                        &
                        ( rad%extkb - rad%extkd ) )
      
      rad%gradis(:,2) = ( 8.0 * C%EMLEAF / ( C%CAPP * air%rho ) ) *            &
                        flpwb / met%tk * rad%extkd *                           &
                        ( 1.0 - rad%transd ) / rad%extkd - rad%gradis(:,1)
      
      ! Longwave radiation absorbed by sunlit canopy fraction:
      rad%qcan(:,1,3) = (rad%flws - flwv ) * rad%extkd *                       &
                        ( rad%transd - rad%transb ) / ( rad%extkb - rad%extkd )&
                        + ( emair- C%EMLEAF ) * rad%extkd * flpwb *            &
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
                   EXP( -rad%extkbm(:,1) * canopy%vlaiw ) +                    &
                   ( 1. - rad%fbeam(:,1) ) * ( 1. - rad%reffdf(:,1) ) *        &
                   EXP( -rad%extkdm(:,1) * canopy%vlaiw ) )                    &
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

! ------------------------------------------------------------------------------
Subroutine EffectiveExtinctCoeffs( EffExtCoeff_beam, EffExtCoeff_dif, mp, nrb, &
                                   sunlit_veg_mask,                        &
                                   ExtCoeff_beam, ExtCoeff_dif, c1 )
implicit none
integer :: mp                   !total number of "tiles"  
integer :: nrb                  !number of radiation bands [per legacy=3, but really=2 VIS,NIR. 3rd dim was for LW]

REAL :: EffExtCoeff_beam(mp,nrb)!Effective Extinction co-efficient for Direct Beam component of SW radiation
REAL :: EffExtCoeff_dif(mp,nrb) !Effective Extinction co-efficient for Diffuse component of SW radiation

REAL :: c1(mp,nrb)
logical :: sunlit_veg_mask(mp)  !combined mask - BOTH sunlit and vegetated
REAL :: ExtCoeff_beam(mp)       !"raw" Extinction co-efficient for Direct Beam component of SW radiation
REAL :: ExtCoeff_dif(mp)        !"raw"Extinction co-efficient for Diffuse component of SW radiation

EffExtCoeff_dif = 0.
call EffectiveExtinctCoeff( EffExtCoeff_dif, mp, ExtCoeff_dif, c1 )

EffExtCoeff_beam = 0.
call EffectiveExtinctCoeff( EffExtCoeff_beam, mp, ExtCoeff_beam, c1, &
                            sunlit_veg_mask )
End Subroutine EffectiveExtinctCoeffs

! modified k diffuse(6.20)(for leaf scattering)
subroutine EffectiveExtinctCoeff(Eff_ExtCoeff, mp, ExtCoeff, c1, mask )
implicit none
integer :: mp 
real :: Eff_ExtCoeff(mp,2) 
real :: ExtCoeff(mp) 
real :: c1(mp,2) 
logical, optional :: mask(mp) 
integer :: i, b

DO i = 1,mp
  DO b = 1, 2
    !IF mask is present we are doing the beam component then: 
    if( present(mask)) then 
      !then ONLY IF it is sunlit and vegetated -else default 
      if( mask(i) ) Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)
    else         
      Eff_ExtCoeff(i,b) = ExtCoeff(i) * c1(i,b)          
    endif

  enddo
enddo

End subroutine EffectiveExtinctCoeff




END MODULE cable_radiation_module
