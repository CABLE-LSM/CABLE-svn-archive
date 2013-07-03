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
! Purpose: Fills CABLE type 'air' with appropriate values calculating 
!          temperature dependent physical constants
!
! Called from: cbm, define_canopy
!
! Contact: Yingping.Wang@csiro.au
!
! History: No significant change from v1.4b
!
!
! ==============================================================================

!air% type variables calculated here and used elsewhere in model
!IN _CANOPY
!air%rho
!air%cmolar
!air%visc
!air%rlam  
!air%dsatdk
!air%psyc
!air%epsi
!
!IN _RAD
!air%rho 
!air%cmolar


MODULE cable_air_module

   ! local pointers to global constants defined in
   USE cable_data_module
   !USE cable_data_module, ONLY : iair_type, point2constants
   
   IMPLICIT NONE

   PUBLIC define_air
   PRIVATE

   ! local pointers setup to to global constants defined in cable_data
   ! to avoid carrying ALL constants througout modules unecessarily
   TYPE iair_type
      REAL, POINTER ::                                                         &
         ! physical constants
         TFRZ, RMAIR, RGAS,                                                    &
         TETENA, TETENB, TETENC,                                               &
         CAPP, RMH2O, HL
   END TYPE iair_type

   TYPE( iair_type ) :: C


CONTAINS

SUBROUTINE define_air(met,air)

   USE cable_def_types_mod,   ONLY : air_type, met_type, mp    
   USE cable_common_module,   ONLY : cable_runtime, cable_user, ktau_gl,       & 
                                     open_code_log, report_min 

   TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
   TYPE (met_type), INTENT(IN)    :: met ! meteorological variables
   
   ! local vatiables 
    REAL, DIMENSION(mp) ::                                                     &
      es,               & ! sat vapour pressure (mb)   
      volm,             & ! molar volume (m3 mol-1)
      qsat,             & ! saturation specific humidity
      above_freezing,   & ! 4 convenience comp. met%tvair - C%TFRZ 
      Ratio_water2Air,  & ! 4 convenience comp. Ratio 
      Teten_ratio,      & ! 4 convenience comp. Ratio 
      Teten_exp           ! 4 convenience comp. Exponential
 
! END header
   
   ! local ptrs to constants defined in cable_data_module
   CALL air_type_ptr 
 
   ! 4 Concenience calculate: 
   ! a few things depend on relative tvair above freezing 
   above_freezing = met%tvair - C%TFRZ 
   Ratio_water2Air = C%RMH2O / C%RMAIR
   Teten_ratio = C%TETENB *  above_freezing / ( C%TETENC +  above_freezing ) 
   Teten_exp = EXP( Teten_ratio ) 

   ! Calculate saturation vapour pressure
   es =  C%TETENA * Teten_exp
   
   ! Calculate conversion factor from from m/s to mol/m2/s
   air%cmolar = met%pmb * 100.0 / (C%RGAS * (met%tvair))
             
   ! saturation specific humidity
   qsat = Ratio_water2Air * es / met%pmb

!jhan: move to driver  
CALL open_code_log( "code_log.txt" ) 


!jhan: why 1.3
   ! Calculate dry air density:
   air%rho = MIN(1.3,C%RMAIR * air%cmolar)
   CALL report_min( "air%rho", "1.3", "RMAIR * cmolar", "" )
  
!jhan: why 100.0
   ! molar volume (m^3/mol)
   volm = C%RGAS * (met%tvair) / (100.0 * met%pmb)

   ! latent heat for water (j/kg)
   air%rlam= C%HL
   
   ! saturation specific humidity
   air%qsat = (C%RMH2O / C%RMAIR) * es / met%pmb
   
   ! d(qsat)/dT ((kg/kg)/K)
   air%epsi = (air%rlam / C%CAPP) * Ratio_water2Air * es * C%TETENB *          &
              C%TETENC / ( C%TETENC + above_freezing ) ** 2 / met%pmb
              
!jhan: why 1.35 etc
   ! air kinematic viscosity (m^2/s)
   air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * above_freezing )

CALL report_min( "air%visc", "1.0", "1.35+.0092+  ", "please explain ..." )
   
   ! psychrometric constant
    air%psyc = met%pmb * 100.0 * C%CAPP / air%rlam / Ratio_water2Air 

   
   ! d(es)/dT (mb/K)
   air%dsatdk = 100.0*(C%TETENA*C%TETENB*C%TETENC)/( above_freezing +      &
                C%TETENC)**2 * Teten_exp

  
END SUBROUTINE define_air

!===============================================================================

SUBROUTINE air_type_ptr
      
   ! local pointers to global constants defined in
   USE cable_data_module, ONLY : PHYS
   
   C%TFRZ   => PHYS%TFRZ
   C%RMAIR  => PHYS%RMAIR
   C%RGAS   => PHYS%RGAS            
   C%TETENA => PHYS%TETENA 
   C%TETENB => PHYS%TETENB 
   C%TETENC => PHYS%TETENC 
   C%CAPP   => PHYS%CAPP
   C%RMH2O  => PHYS%RMH2O
   C%HL     => PHYS%HL 

END SUBROUTINE air_type_ptr

!===============================================================================

END MODULE cable_air_module
















