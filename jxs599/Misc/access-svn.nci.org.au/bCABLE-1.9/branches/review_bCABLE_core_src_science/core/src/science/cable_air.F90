
!===COPYRIGHT==================================================================
! The source codes are part of the australian 
! Community Atmosphere Biosphere Land Exchange (CABLE) model. 
! Please register online at xxx and sign the agreement before use 
! contact: whox@xxxx.yyy about registration user agreement
!==============================================================================


!==============================================================================
! Name: air_module 
! Purpose: calculate air properties for 
! CALLed from: executed PROGRAM 
! MODULEs used:  physical_constants
!                define_types
!                define_dimensions
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


MODULE air_module
   IMPLICIT NONE
   PRIVATE
   PUBLIC define_air


CONTAINS


SUBROUTINE define_air(met,air)
   USE physical_constants
   USE define_types,          ONLY : air_type, met_type    
   USE define_dimensions,     ONLY : mp, r_1
   USE cable_common_module,   ONLY : cable_runtime, cable_user, ktau_gl 
   USE cable_diag_module,     ONLY : cable_stat

   TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
   TYPE (met_type), INTENT(IN)    :: met ! meteorological variables
   
   ! local vatiables 
   REAL(r_1), DIMENSION(mp)     :: es ! sat vapour pressure (mb)   
   INTEGER :: i 
      
   ! END header
  
   
   IF( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
      CALL cable_stat('define_air')

    ! Calculate saturation vapour pressure
    es = tetena * EXP(tetenb * (met%tvair-tfrz)/(tetenc + (met%tvair-tfrz)))
    
    ! Calculate conversion factor from from m/s to mol/m2/s
    air%cmolar = met%pmb * 100.0 / (rgas * (met%tvair))
    
    ! Calculate dry air density:
    air%rho = MIN(1.3,rmair * air%cmolar)
    
    ! molar volume (m^3/mol)
    air%volm = rgas * (met%tvair) / (100.0 * met%pmb)
    
    ! latent heat for water (j/kg)
    air%rlam= hl
    
    ! saturation specific humidity
    air%qsat = (rmh2o / rmair) * es / met%pmb
    
    ! d(qsat)/dT ((kg/kg)/K)
    air%epsi = (air%rlam / capp) * (rmh2o / rmair) * es * tetenb * tetenc / &
               (tetenc + (met%tvair - tfrz)) ** 2 / met%pmb
    
    ! air kinematic viscosity (m^2/s)
    air%visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met%tvair - tfrz))
    
    ! psychrometric constant
    air%psyc = met%pmb * 100.0 * capp * rmair / air%rlam / rmh2o
    
    ! d(es)/dT (mb/K)
    air%dsatdk = 100.0*(tetena*tetenb*tetenc)/((met%tvair-tfrz)+tetenc)**2 &
                 * EXP(tetenb*(met%tvair-tfrz)/((met%tvair-tfrz) + tetenc))
  
END SUBROUTINE define_air


END MODULE air_module
