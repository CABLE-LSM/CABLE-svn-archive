!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
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
! Purpose: Calculate roughness lengths as a function of soil and canopy 
!          parameters
!
! Contact: Eva.Kowalczyk@csiro.au
!
! History: No significant changes since v1.4b except change to cope with 
!          split timestep in ACCESS (zref_uv, zref_tq)
!
!
! ==============================================================================

MODULE cable_roughness_module
   
   USE cable_data_module, ONLY : irough_type, point2constants
   
   IMPLICIT NONE
   
   TYPE ( irough_type ) :: C 
   PRIVATE
   PUBLIC ruff_resist

CONTAINS

!-------------------------------------------------------------------------------

SUBROUTINE ruff_resist(veg, rough, ssnow, canopy)

   ! m.r. raupach, 24-oct-92
   ! see: Raupach, 1992, BLM 60 375-395
   !      MRR notes "Simplified wind model for canopy", 23-oct-92
   !      MRR draft paper "Simplified expressions...", dec-92
   ! modified to include resistance calculations by Ray leuning 19 Jun 1998  

   USE cable_common_module, ONLY : cable_runtime, cable_user,                  &
                                   ktau => ktau_gl,                            &
                                   kend => kend_gl, knode_gl
   USE cable_def_types_mod, ONLY : veg_parameter_type, roughness_type,         &
                                   soil_snow_type, canopy_type, mp  
   USE cable_diag_module
   
   TYPE(roughness_type), INTENT(INOUT) :: rough
!pass these vars from types in cbm
   TYPE (canopy_type),   INTENT(INOUT) :: canopy
   TYPE(soil_snow_type), INTENT(IN)    :: ssnow
   TYPE (veg_parameter_type),  INTENT(INOUT) :: veg

   REAL, DIMENSION(mp) ::                                                      &
      xx,      & ! =C%CCD*LAI; working variable 
      dh         ! d/h where d is zero-plane displacement
                                    ! tiles belonging to the same grid
   REAL, DIMENSION(mp) ::                                                      &
         za_uv,   & ! level of lowest atmospheric model layer
         za_tq      ! level of lowest atmospheric model layer
   
   REAL, DIMENSION(mp) ::                                                      &
      Eff_SnowDensity, & 
      Eff_SnowDepth, & 
      Eff_LAI, &
      Eff_height
   
   REAL, PARAMETER ::                                                          &
      z0soil = 1.e-6    ! Roughness length of bare soil (m)
   
!jhan: make these scalars
   REAL, DIMENSION(mp) ::                                                      &
      term2, term3, term5, term6 ! for aerodyn resist. calc.

   REAL, DIMENSION(mp) :: xx_term

   LOGICAL, DIMENSION(mp) :: BareSoil_mask

   INTEGER, SAVE :: iDiag1=0      

!jhan: do same thing here as in cable_air,i.e.bring subr into local module 
   CALL point2constants( C ) 

!...............................................................................

   ! Reference height zref is height above the displacement height
   ! za_uv,tq are forcing
!jhan: isnt this very artificial, imposing ahard minimum      
   rough%zref_uv = MAX( 3.5, rough%za_uv )
   rough%zref_tq = MAX( 3.5, rough%za_tq )

   call cable_diag( iDiag1, "za_uv", mp, kend, ktau,                   &
                    knode_gl, "za_uv", rough%za_uv )

!...............................................................................

   ! Roughness length of bare snow (m): ?BUT we are calcuating this everywhere
!jhan: this stops roughness length of bare snow from going negative
!jhan: isnt it better to just range this from 0:calc'ed   
   ! in this case we consider snow density doesnt go above 20 i.e.(0:20) 
   Eff_SnowDepth = MIN( ssnow%snowd, 20. )

!jhan: why not just put z0soil =0
!jhan: i.e.  rough%z0soilsn = z0soil - 5.e-8 * ssnow%snowd
!jhan: but why 5.e-8
!fair enough is less friction than soil but by how much (5.e-8?)

!snow depth can be zero -> rough%z0soilsn = z0soil = 1.e-6 
!either way if eff_snowd=20, %z0soiln= 2.e1 * 5e-8 = 1.e-7 

   ! Adjusted Roughness length of bare snow (m), i.e.(9:10)e-7
   rough%z0soilsn = z0soil - 5.e-8 * Eff_SnowDepth

!jhan: this will never happen!! delete
   ! Restricted Roughness length of bare snow (m):
   rough%z0soilsn = MAX( rough%z0soilsn, 1.e-8 )

!...............................................................................

   ! Set canopy height above snow level:
   
   ! ssnow%ssdn = snow density
   ! in this case we consider snow density doesnt go below 100  
   Eff_SnowDensity = MAX( ssnow%ssdnn, 100. )  

   ! Set canopy height above snow level:
   ! if there is no snow = veg%hc
   ! veg%hc = roughness height of canopy (veg-snow) ?comment? 
!jhan: BUt surely this can be 0. or bare snow (z0soilsn) if trees are 
! buried under snow 
!jhan: ?1.2?
!BUT this isn't dimensionally sound
   rough%hruff = veg%hc - 1.2 * ssnow%snowd / Eff_SnowDensity  
   
!jhan: USE report_max subr
   ! rough%hruff is limited here to not be less than 0.01 
   rough%hruff = MAX( 0.01, rough%hruff ) 
   
!...............................................................................

!jhan: can this effective height be used to get hruff above  i.e. 
   !rough%hruff = Eff_height - 1.2 * ssnow%snowd / Eff_SnowDensity  
   Eff_height = MAX( 0.01, veg%hc )
   
   ! LAI decreases due to snow and vegetation fraction:
!jhan: chheck this
   ! veg%vlai is the input LAI
!jhan: LAI modified by roughness height relative to real height ? why?
!ratio of exposed to buried affects LAI same ratio
   canopy%vlaiw = veg%vlai * rough%hruff / Eff_height 

   ! By default rghLAI = snow adjusted LAI 
   canopy%rghlai = canopy%vlaiw

!...............................................................................
   
   ! We need to compute a rghLAI to include in vlaiw that  
   ! Where ~NO snow AND is not broadleaf evergreen? forest
   ! rghLAI does not exceed 3 Why?
   WHERE( ssnow%snowd .LT. 0.001 .AND. veg%iveg .NE. 1 )                       &
      canopy%rghlai = MIN( 3., canopy%rghlai )

!...............................................................................

!jhan: put LOGICAL mask.criteria might be elsewhere?
!jhan: How can second case (.OR.)ever happen 
!jhan: z0soilsn=0 if %ssnowd =0, BUT then restricted to MIN of .1e-8 
!jhan: hruff = veg height if no snow, AND restricted to MIN of .01 
!jhan: canopy height above snow LT roughness of bare snow (which should never happen as roughness length ~10% of height) 
  
   ! Effective LAI to consider in calc of friction velocity 
   Eff_LAI = canopy%vlaiw * 0.5

   ! Mask of Bare Soil surfaces
   BareSoil_mask = canopy%vlaiw .LT. 0.01 .OR.                                 &
                   rough%hruff .LT. rough%z0soilsn 

   ! set exposed height to zero when bare soil anyway 
   WHERE( BareSoil_mask )                                                      &
      rough%hruff = 0.0
  
   ! set Effective LAI to consider for VEGETATED SURFACEs
   WHERE( .NOT. BareSoil_mask )                                                &
      Eff_LAI = canopy%rghlai * 0.5
     
!...............................................................................

   ! Friction velocity/windspeed at canopy height
   ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
   ! (C%USUHM set in physical_constants module):
   ! Effective LAI to consider here
   rough%usuh = SQRT( C%CSD + C%CRD * Eff_LAI )
   rough%usuh = MIN( rough%usuh, C%USUHM )

!...............................................................................

!jhan:explain     
   xx_term = MAX( Eff_LAI, 0.0005 ) 
   xx = SQRT( C%CCD * xx_term )
    
   ! Displacement height/canopy height:
   ! eq.8 Raupach 1994, BLM, vol 71, p211-216
   dh = 1.0 - ( 1.0 - EXP( -xx ) ) / xx

   ! Extinction coefficient for wind profile in canopy:
   ! eq. 3.14, SCAM manual (CSIRO tech report 132)
   rough%coexp = rough%usuh / ( C%VONK * C%CCW_C * ( 1.0 - dh ) )

!-------------------------------------------------------------------------------

   ! These initializations are over-written on vegetated surfaces

   ! zero-plane displacement
   rough%disp = 0.0
   
   ! SCALAR Roughness sublayer depth (ground=origin)
   rough%zruffs = 0.0
   
   ! set to roughness length here EQV to bare snow
   rough%z0m = rough%z0soilsn

   ! eq. 3.54, SCAM manual (CSIRO tech report 132)
   rough%rt0us = 0.0  

   ! resistance from disp to hruf
   rough%rt1usa = 0.0 
   
   ! resist fr hruf to zruffs (zref if zref<zruffs)
   rough%rt1usb = 0.0

!...............................................................................

   ! set roughness AND resistance(s) for VEGETATED SURFACEs
   WHERE( .NOT. BareSoil_mask )
      
      ! Calculate zero-plane displacement:
      rough%disp = dh * rough%hruff
       
      ! Calcualte roughness length:
      rough%z0m = ( (1.0 - dh) * EXP( LOG( C%CCW_C ) - 1. + 1. / C%CCW_C       &
                  - C%VONK / rough%usuh ) ) * rough%hruff
       
      term2  = EXP( 2 * C%CSW * canopy%rghlai *                          &
                     ( 1 - rough%disp / rough%hruff ) )

      term3  = C%A33**2 * C%CTL * 2 * C%CSW * canopy%rghlai
      
      term5  = MAX( ( 2. / 3. ) * rough%hruff / rough%disp, 1.0 )
      
      term6 =  EXP( 3. * rough%coexp * ( rough%disp / rough%hruff -1. ) )
      
      ! eq. 3.54, SCAM manual (CSIRO tech report 132)
      rough%rt0us  = term5 * ( C%ZDLIN * LOG(                            &
                     C%ZDLIN * rough%disp / rough%z0soilsn ) +                 &
                     ( 1 - C%ZDLIN ) )                                         &
                     * ( EXP( 2 * C%CSW * canopy%rghlai )  -  term2 )    &
                     / term3  
      
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
      rough%zruffs = rough%disp + rough%hruff * C%A33**2 * C%CTL / C%VONK /    &
                     term5
      
      ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
      rough%rt1usa = term5 * ( term2 - 1.0 ) / term3
      
      rough%rt1usb = term5 * ( MIN( rough%zref_tq + rough%disp,          &
                     rough%zruffs ) - rough%hruff ) /                          &
                     ( C%A33**2 * C%CTL * rough%hruff )

      rough%rt1usb = MAX( rough%rt1usb, 0.0 ) ! in case zrufs < rough%hruff
    
    END WHERE




END SUBROUTINE ruff_resist

!-------------------------------------------------------------------------------

FUNCTION exponentialGrowth( usuh ) RESULT (zexp)

   REAL, INTENT(IN) :: usuh
   REAL :: factor, zexp 

         factor = LOG( C%CCW_C ) - 1. + ( 1. / C%CCW_C )             &
                            - ( C%VONK / usuh )

         zexp= EXP( factor )

END FUNCTION exponentialGrowth 

!-------------------------------------------------------------------------------

END MODULE cable_roughness_module
