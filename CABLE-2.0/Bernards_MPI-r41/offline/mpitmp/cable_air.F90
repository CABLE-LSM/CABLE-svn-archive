
!===COPYRIGHT==================================================================
! The source codes are part of the australian 
! Community Atmosphere Biosphere Land Exchange (CABLE) model. 
! Please register online at xxx and sign the agreement before use 
! contact: whox@xxxx.yyy about registration user agreement
!==============================================================================


!==============================================================================
! Name: cable_air_module 
! Purpose: calculate air properties for 
! CALLed from: executed PROGRAM 
! MODULEs used:  cable_def_types_mod
!                cable_common_module
! 
! 
! CALLs: None      
!
!
! Major contribution: land surface modeling team, CSIRO, Aspendale
!            
!==============================================================================


!==============================================================================
! changes since version release on 
! changes made by who on date
!
!==============================================================================


MODULE cable_air_module

   ! local pointers to global constants defined in
   USE cable_data_module, ONLY : iair_type, point2constants
   
   IMPLICIT NONE

   PUBLIC define_air
   PRIVATE

   TYPE( iair_type ) :: C


CONTAINS


SUBROUTINE define_air(met,air)

   USE cable_def_types_mod,          ONLY : air_type, met_type, mp    
   USE cable_common_module,   ONLY : cable_runtime, cable_user, ktau_gl 

   TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
   TYPE (met_type), INTENT(IN)    :: met ! meteorological variables
   
   ! local vatiables 
   REAL, DIMENSION(mp)     :: es ! sat vapour pressure (mb)   
   INTEGER :: i 
      
   ! END header
   
   ! local ptrs to constants defined in cable_data_module
   call point2constants( C )
   
   ! Calculate saturation vapour pressure
   es = C%TETENA * EXP( C%TETENB * ( met%tvair - C%TFRZ )                     &
        / ( C%TETENC + ( met%tvair - C%TFRZ ) ) )
   
   ! Calculate conversion factor from from m/s to mol/m2/s
   air%cmolar = met%pmb * 100.0 / (C%RGAS * (met%tvair))
   
   ! Calculate dry air density:
   air%rho = MIN(1.3,C%RMAIR * air%cmolar)
   
   ! molar volume (m^3/mol)
   air%volm = C%RGAS * (met%tvair) / (100.0 * met%pmb)
   
   ! latent heat for water (j/kg)
   air%rlam= C%HL
   
   ! saturation specific humidity
   air%qsat = (C%RMH2O / C%RMAIR) * es / met%pmb
   
   ! d(qsat)/dT ((kg/kg)/K)
   air%epsi = (air%rlam / C%CAPP) * (C%RMH2O / C%RMAIR) * es * C%TETENB *     &
              C%TETENC / ( C%TETENC + (met%tvair - C%TFRZ) ) ** 2 / met%pmb
   
   ! air kinematic viscosity (m^2/s)
   air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - C%TFRZ) )
   
   ! psychrometric constant
   air%psyc = met%pmb * 100.0 * C%CAPP * C%RMAIR / air%rlam / C%RMH2O
   
   ! d(es)/dT (mb/K)
   air%dsatdk = 100.0*(C%TETENA*C%TETENB*C%TETENC)/((met%tvair-C%TFRZ) +      &
                C%TETENC)**2 * EXP( C%TETENB * ( met%tvair-C%TFRZ ) /         &
                ( (met%tvair-C%TFRZ) + C%TETENC) )
  
END SUBROUTINE define_air



END MODULE cable_air_module
