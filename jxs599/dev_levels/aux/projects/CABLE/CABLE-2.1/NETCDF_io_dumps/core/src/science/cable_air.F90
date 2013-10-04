!jhan:include knode capability - and how much can i keep? test this on vayu with UM aswell
module air_module
   implicit none
   PUBLIC define_air
   PRIVATE

   interface define_air
      module procedure  define_air_data
   end interface define_air
   
CONTAINS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   SUBROUTINE define_air_data(met,air)
      use cable_common_module, only : cable_user
      use define_types, only : air_type, met_type    
      use cable_data_module, only : const, air_in, air_out, air_const
      use cable_dumps_module, only :   air_in_dump, &
                                       air_out_dump, &
                                       air_in_read

      implicit none
      TYPE (air_type), INTENT(INOUT) :: air ! air_type variables
      TYPE (met_type), INTENT(IN)  :: met ! meteorological variables
      type (air_in),dimension(:), allocatable, save :: i
      type (air_out), save :: o
      type (air_const), save :: c
      integer, save :: n_call =1
      integer :: j
      character(len=20) :: ncfile, tag 
      !jhan: work these out properly
      real, parameter, dimension(5) :: gauss_wt = (/.6,.8,1.,1.2,1.4/)
         ncfile = "air_in.nc"

      if(n_call == 1) then      
         allocate(i(cable_user%AIR_IN_J))     
         do j=1, cable_user%AIR_IN_J     
            i(j)%met_tvair    => met%tvair 
            i(j)%met_pmb      => met%pmb
            i(j)%met_tvair    = i(j)%met_tvair * gauss_wt(j) 
            i(j)%met_pmb      = i(j)%met_pmb  * gauss_wt(j) 
         enddo
         o%air_cmolar   => air%cmolar
         o%air_rho      => air%rho          
         o%air_rlam     => air%rlam          
         o%air_epsi     => air%epsi          
         o%air_visc     => air%visc         
         o%air_psyc     => air%psyc         
         o%air_dsatdk   => air%dsatdk       
    
         c%capp    => const%phys%capp  
         c%hl      => const%phys%hl       
         c%rgas    => const%phys%rgas  
         c%rmair   => const%phys%rmair 
         c%rmh2o   => const%phys%rmh2o 
         c%tfrz    => const%phys%tfrz  
         c%tetena  => const%phys%tetena  
         c%tetenb  => const%phys%tetenb  
         c%tetenc  => const%phys%tetenc  
      endif
     
      do j=1, cable_user%AIR_IN_J     
!         write(tag,11) j 
!   11    format(i2.2)  
!         ncfile = trim("air_in"//trim(tag))//".nc"

         if( cable_user%AIR_IN_DUMP ) &    
            call air_in_dump(i(j), n_call, cable_user%AIR_IN_J, ncfile,j)
         
         if( cable_user%AIR_IN_READ ) &    
            call air_in_read(i(j), n_call)
   
         call define_air_model(i(j),o,c)
   
         if( cable_user%AIR_OUT_DUMP ) &    
            call air_out_dump(o, n_call)
      enddo
     
      n_call = n_call + 1

      return
   END SUBROUTINE define_air_data
      

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine define_air_model(i,o,c)
      use cable_common_module, only : cable_user
      use cable_diag_module, only : cable_stat
      use cable_data_module, only : const, air_in, air_out, air_const
      use define_dimensions, only : mp, r_1
      implicit none
      !type (air_in), intent(in), dimension(2) :: i
      type (air_in), intent(in) :: i
      type (air_out), intent(out) :: o
      type (air_const),intent(in) :: c
      integer :: j 
      real, dimension(mp)     :: &
         es,   &    !sat vapour pressure (mb)   
         !used in Mk3L driver. make local  
         volm, &   
         qsat

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
           call cable_stat('define_air')
    
         ! Calculate saturation vapour pressure
         es = c%tetena * EXP(c%tetenb * (i%met_tvair-c%tfrz)/(c%tetenc + (i%met_tvair-c%tfrz)))
     
         ! Calculate conversion factor from from m/s to mol/m2/s
         o%air_cmolar = i%met_pmb * 100.0 / (c%rgas * (i%met_tvair))
         
         ! Calculate dry air density:
         o%air_rho = min(1.3,c%rmair * o%air_cmolar)
         
         ! molar volume (m^3/mol)
         volm = c%rgas * (i%met_tvair) / (100.0 * i%met_pmb)
   
         ! latent heat for water (j/kg)
         !o%air%rlam = (2501.0 - 2.38 * (met%tvair- tfrz)) * 1000.0
         o%air_rlam= c%hl
         ! saturation specific humidity
         qsat = (c%rmh2o / c%rmair) * es / i%met_pmb
   
         ! d(qsat)/dT ((kg/kg)/K)
         o%air_epsi = (o%air_rlam / c%capp) * (c%rmh2o / c%rmair) * es * c%tetenb * c%tetenc / &
              &  (c%tetenc + (i%met_tvair - c%tfrz)) ** 2 / i%met_pmb
         
         ! air kinematic viscosity (m^2/s)
         o%air_visc = 1e-5 * MAX(1.0, 1.35 + 0.0092 * (i%met_tvair - c%tfrz))
         
         ! psychrometric constant
         o%air_psyc = i%met_pmb * 100.0 * c%capp * c%rmair / o%air_rlam / c%rmh2o
         
         ! d(es)/dT (mb/K)
         o%air_dsatdk = 100.0*(c%tetena*c%tetenb*c%tetenc)/((i%met_tvair-c%tfrz)+c%tetenc)**2 &
              * EXP(c%tetenb*(i%met_tvair-c%tfrz)/((i%met_tvair-c%tfrz) + c%tetenc))
      return
   end subroutine define_air_model

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end module air_module
 
