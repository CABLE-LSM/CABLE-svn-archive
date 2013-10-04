module cable_wrappers_module
   implicit none
   PUBLIC define_air
   PRIVATE

   interface define_air
      module procedure  define_air_cable, define_air_comp
   end interface define_air
   
CONTAINS

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ wrapper subrs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!jhan: we can put this into  eventually
   subroutine define_air_comp(mp,kend)
      use cable_data_module, only : air_in, air_out
      use cable_dumps_module, only :   dump_air_in_ncdf, &
                                       read_air_in_ncdf
      use air_module, only : define_air
      implicit none
      integer, intent(in) :: mp
      integer, intent(in) :: kend 
      type (air_in),dimension(2),save :: i
      type (air_out), save :: o
      integer, save :: n_call =1
      integer :: j

      if(n_call == 1) then
         allocate( i(1)%met_tvair(mp) )
         allocate( i(1)%met_pmb(mp) )    
         
         allocate( o%air_cmolar(mp) )  
         allocate( o%air_rho(mp) )       
         allocate( o%air_rlam(mp) )       
         allocate( o%air_epsi(mp) )       
         allocate( o%air_visc(mp) )      
         allocate( o%air_psyc(mp) )      
         allocate( o%air_dsatdk(mp) )    
      endif    
      
      do j=1, kend
         !call read_air_in(i(1))

         call new_define_air(i(1),o)

         !call dump_air_out(o)
      enddo

      n_call = n_call + 1

      return
   end subroutine define_air_comp
      
end module cable_wrappers_module
