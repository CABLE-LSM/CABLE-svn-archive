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
! Purpose: Calculates surface exchange fluxes through the solution of surface 
!          energy balance and its interaction with plant physiology. Specific 
!        representation of the transport of scalars within a canopy is included.
!
! Called from: cbm
!
! Contact: Yingping.Wang@csiro.au and Eva.Kowalczyk@csiro.au
!
! History: Revision of canopy temperature calculation (relative to v1.4b) 
!          Reorganisation of code (dryLeaf, wetLeaf, photosynthesis subroutines
!          taken out of define_canopy)
!
!
! ==============================================================================

MODULE cable_canopy_module
   
   USE cable_data_module, ONLY : icanopy_type, point2constants 
   
   IMPLICIT NONE
   
   PUBLIC define_canopy
   PRIVATE
   
   TYPE( icanopy_type ) :: C
  
     
CONTAINS
 

SUBROUTINE define_canopy(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)
   USE cable_def_types_mod
   !USE cable_radiation_module
   !USE cable_air_module
   !USE cable_common_module   
   !USE cable_roughness_module
   USE CABLE_canopy_main_module

   TYPE (balances_type), INTENT(INOUT)  :: bal
   TYPE (radiation_type), INTENT(INOUT) :: rad
   TYPE (roughness_type), INTENT(INOUT) :: rough
   TYPE (air_type), INTENT(INOUT)       :: air
   TYPE (met_type), INTENT(INOUT)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy

   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT)    :: veg
   REAL, INTENT(IN)               :: dels ! integration time setp (s)

   ! END header
   
   CALL cable_canopy_main(bal,rad,rough,air,met,dels,ssnow,soil,veg, canopy)

END SUBROUTINE define_canopy

! -----------------------------------------------------------------------------


    
END MODULE cable_canopy_module
