
MODULE air_module
   use define_dimensions, only : mp, r_1
   IMPLICIT NONE
   PUBLIC define_air
   PRIVATE

   type air_in
      real(r_1), dimension(:), pointer     :: &
         met_tvair, &
         met_pmb
   end type air_in

   type air_out
      real(r_1), dimension(:), pointer     :: &
         air_cmolar, &
         air_rho, &
         air_rlam, &
         air_dsatdk, &
         air_epsi, & 
         air_visc, & 
         air_psyc 
   end type air_out


CONTAINS


   SUBROUTINE define_air(met,air)
      use define_types, only : air_type, met_type    
      implicit none
      TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
      TYPE (met_type), INTENT(IN)  :: met ! meteorological variables
      
      call define_air_wrapper(met,air)

      return
   END SUBROUTINE define_air




   SUBROUTINE define_air_wrapper(met,air)
      use define_types, only : air_type, met_type    
      use define_dimensions, only : mp, r_1
      implicit none
      TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
      TYPE (met_type), INTENT(IN)  :: met ! meteorological variables
      type (air_in) :: i
      type (air_out) :: o

      i%met_tvair    => met%tvair 
      i%met_pmb       => met%pmb 

      o%air_epsi     = 0.; o%air_visc     = 0.
      o%air_psyc     = 0.; o%air_dsatdk   = 0.   
      
      !filename =    
      !open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
      !   iostat=gopenstatus, form="unformatted", position='append' )
      !if(gopenstatus==gok) then
      !      write (713942) var1
      !else
      !   write (*,*) filename//'.bin',' NOT open for write. Error'
      !endif
      !close(713942)

      call new_define_air(i,o)

      !debrief air type and return
      return
   END SUBROUTINE define_air_wrapper
      



   subroutine new_define_air(i,o)
      use cable_common_module, only : cable_user
      use cable_diag_module, only : cable_stat
      use cable_data_module, only : const
      use define_dimensions, only : mp, r_1
      implicit none
      type (air_in), intent(in) :: i
      type (air_out), intent(out) :: o
      integer :: j 
      real(r_1), dimension(mp)     :: &
         es,   &    !sat vapour pressure (mb)   
         !used in Mk3L driver. make local  
         volm, &   
         qsat

      real(r_1), pointer     :: &
         CAPP, &  
         HL, &      
         RGAS, &  
         RMAIR, & 
         RMH2O, & 
         TFRZ, &  
         TETENA, & 
         TETENB, & 
         TETENC 




         CAPP    => const%phys%capp  
         HL      => const%phys%hl       
         RGAS    => const%phys%rgas  
         RMAIR   => const%phys%rmair 
         RMH2O   => const%phys%rmh2o 
         TFRZ    => const%phys%tfrz  
         TETENA  => const%phys%tetena  
         TETENB  => const%phys%tetenb  
         TETENC  => const%phys%tetenc  

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
           call cable_stat('define_air')
    
         ! Calculate saturation vapour pressure
         es = TETENA * EXP(TETENB * (i%met_tvair-TFRZ)/(TETENC + (i%met_tvair-TFRZ)))
     
              !do i=1,mp
              !print *, 'jhan:air', ktau_gl, i, met%tvair(i) 
              !enddo
         
         ! Calculate conversion factor from from m/s to mol/m2/s
         o%air_cmolar = i%met_pmb * 100.0 / (RGAS * (i%met_tvair))
         
         ! Calculate dry air density:
         o%air_rho = MIN(1.3,RMAIR * o%air_cmolar)
         
         ! molar volume (m^3/mol)
         volm = RGAS * (i%met_tvair) / (100.0 * i%met_pmb)
   
         ! latent heat for water (j/kg)
         !o%air%rlam = (2501.0 - 2.38 * (met%tvair- tfrz)) * 1000.0
         o%air_rlam= HL
         ! saturation specific humidity
         qsat = (RMH2O / RMAIR) * es / i%met_pmb
   
         ! d(qsat)/dT ((kg/kg)/K)
         o%air_epsi = (o%air_rlam / CAPP) * (RMH2O / RMAIR) * es * TETENB * TETENC / &
              &  (TETENC + (i%met_tvair - TFRZ)) ** 2 / i%met_pmb
         
         ! air kinematic viscosity (m^2/s)
         o%air_visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (i%met_tvair - TFRZ))
         
         ! psychrometric constant
         o%air_psyc = i%met_pmb * 100.0 * CAPP * RMAIR / o%air_rlam / RMH2O
         
         ! d(es)/dT (mb/K)
         o%air_dsatdk = 100.0*(TETENA*TETENB*TETENC)/((i%met_tvair-TFRZ)+TETENC)**2 &
              * EXP(TETENB*(i%met_tvair-TFRZ)/((i%met_tvair-TFRZ) + TETENC))
   end subroutine new_define_air

END MODULE air_module

 
 !   local decs
 !  
 !  RHS
 !  LHS - check where else these are used


 !     !used in canopy. Mk3L driver 
 !     air%epsi 
 !     air%visc 
 !     air%psyc 
 !     air%dsatdk 

 !     !used in canopy & radiation. Mk3L driver 
 !     air%cmolar  

 !     !used in canopy & radiation. Mk3L. UM drivers 
 !     air%rho     

 !     !used in canopy & radiation. Mk3L. UM drivers. cable_checks 
 !     air%rlam   
  
