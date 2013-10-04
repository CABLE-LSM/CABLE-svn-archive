!===COPYRIGHT==================================================================
! The source codes are part of the australian 
! Community Atmosphere Biosphere Land Exchange (CABLE) model. 
! Please register online at xxx and sign the agreement before use 
! contact: whox@xxxx.yyy about registration user agreement
! ==============================================================================


!==================================================================
! Name: XXXX
! Purpose:
! modules used:
! Major contribution: land surface modeling team, CSIRO, Aspendale
! input  file: list of file names
! output file: list of file names
!==================================================================
! changes since version release on 
! changes made by who on date
!
!==================================================================


module air_module
   implicit none
   public define_air
   private

contains


!==================================================================
! Name: XXXX
! Type: input/output or data declaration or initilization or simulation
! Purpose:
! related subroutines:
! Major contribution: land surface modeling team, CSIRO, Aspendale
! input  data: list of variables including dimensions
! output data: list of variables including dimensions
! input  file: list of file names
! output file: list of file names
!==================================================================
!==================================================================
! changes since version release on 
! changes made by who on date
!
!==================================================================


   subroutine define_air( a, met_tvair, met_pmb, cp )
      use cable_common_module, only : cable_user
      use cable_diag_module, only : cable_stat
      use cable_data_module, only : air_auto_type, physical_constants
      use define_dimensions, only : mp, r_1
      implicit none
      !___input args___
      ! type air_auto_type declared in cable_data. 
      ! calculated here, used elsewhere in CABLE, then flushed 
      type (air_auto_type), intent(out) :: a
      
      ! passed variables declared in cable_data, in type met_in_type. 
      ! input forcing, used here and elsewhere in CABLE
      real, dimension(mp), intent(in) :: met_tvair, & ! within canopy air temperature (oK)
                                      :: met_pmb      ! surface air pressure (mbar)

      ! type physical_constants declared and defined in cable_data. 
      type (physical_constants), intent(in) :: cp

      !___local variables___
      real, dimension(mp)     :: &
         es,   &    !sat vapour pressure (mb)   
         !used in Mk3L driver. make local  
         volm, &   
         qsat

         
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
           call cable_stat('define_air')

         ! Calculate saturation vapour pressure
         es = cp%tetena * EXP(cp%tetenb * (met_tvair - cp%tfrz)/(cp%tetenc + (met_tvair - cp%tfrz)))
     
         ! Calculate conversion factor from from m/s to mol/m2/s
         a%air_cmolar = met_pmb * 100.0 / (cp%rgas * (met_tvair))
         
         ! Calculate dry air density:
         a%air_rho = MIN(1.3,cp%rmair * a%air_cmolar)
         
         ! molar volume (m^3/mol)
         volm = cp%rgas * (met_tvair) / (100.0 * met_pmb)
   
         ! latent heat for water (j/kg)
         !a%air%rlam = (2501.0 - 2.38 * (met%tvair- tfrz)) * 1000.0
         a%air_rlam = cp%hl
         ! saturation specific humidity
         qsat = (cp%rmh2o / cp%rmair) * es / met_pmb
   
         ! d(qsat)/dT ((kg/kg)/K)
         a%air_epsi = (a%air_rlam / cp%capp) * (cp%rmh2o / cp%rmair) * es * cp%tetenb * cp%tetenc / &
              &  (cp%tetenc + (met_tvair - cp%tfrz)) ** 2 / met_pmb
         
         ! air kinematic viscosity (m^2/s)
         a%air_visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (met_tvair - cp%tfrz))
         
         ! psychrometric constant
         a%air_psyc = met_pmb * 100.0 * cp%capp * cp%rmair / a%air_rlam / cp%rmh2o
         
         ! d(es)/dT (mb/K)
         a%air_dsatdk = 100.0*(cp%tetena*cp%tetenb*cp%tetenc)/((met_tvair - cp%tfrz) + cp%tetenc)**2 &
              * EXP(cp%tetenb*(met_tvair-cp%tfrz)/((met_tvair-cp%tfrz) + cp%tetenc))
!
      return
   end subroutine define_air


end module air_module


