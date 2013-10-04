
!========================================================================!
!=== PURPOSE: to enable accessibility of fundamental vars in global   ===!
!=== model thru out program (mainly CABLE).                           ===!
!=== USE: use module in subroutines (SRs) at top level of global model===!
!--- 2 define cable_timestep_ data with vars peculiar to global model ===!
!=== optionally call alias_ SR to give name smore familiar to cable   ===!
!=== people. then again use this module in  any subsequent SR in which===!
!=== you want to access this data.                                    ===!  
!========================================================================!

module cable_common_module
   implicit none 

   !---allows reference to "gl"obal timestep in run (from atm_step)
   !---total number of timesteps, and processing node 
   integer, save :: cable_gltimestep_i, cable_gltimestep_tot, &
         cable_gltimestep_width, cable_glnode_i
   integer, save :: ktau_gl, kend_gl, knode_gl, kwidth_gl
   
   !---CABLE runtime switches def in this type
   type kbl_internal_switches
      logical :: um, um_explicit, um_implicit, um_radiation, um_hydrology
      logical :: offline, mk3l
   end type kbl_internal_switches 

   type (kbl_internal_switches), save :: cable_runtime

   !---CABLE runtime switches def in this type
   type kbl_user_switches
      character(len=200) :: VEG_PARS_FILE 
      character(len=20) :: DIAG_SOIL_RESP 
      character(len=200) :: LEAF_RESPIRATION
      character(len=200) :: FWSOIL_SWITCH 
      character(len=3) :: RUN_DIAG_LEVEL 
   end type kbl_user_switches

   type (kbl_user_switches), save :: cable_user

interface jhprint
   module procedure  jhprint_r, jhprint_2r, jhprint_2rr, jhprint_i, jhprint_2i
end interface jhprint

!   !---parameters, tolerances, etc. could be set in _directives.h
!   real, parameter :: RAD_TOLS = 1.0e-2

!jhan:temporary measure. improve hiding
!   real, dimension(:,:), pointer,save :: c1, rhoch
      
   contains
   
   subroutine alias_cable_timestep()
      implicit none 
      ktau_gl = cable_gltimestep_i
      kwidth_gl = cable_gltimestep_width
      kend_gl = cable_gltimestep_tot
      knode_gl = cable_glnode_i
      return
   end subroutine alias_cable_timestep

subroutine jhprint_i(com, var)
   chaRACTER(len=*) :: com
   integer:: var
   print *, '' 
   print *, 'jhan: ',com, var
   print *, '' 
end subroutine jhprint_i

subroutine jhprint_2i(com, var)
   chaRACTER(len=*) :: com
   integer:: var(:)
   print *, '' 
   print *, 'jhan: ',com, var
   print *, '' 
end subroutine jhprint_2i

subroutine jhprint_r(com, var)
   chaRACTER(len=*) :: com
   real :: var
   print *, '' 
   print *, 'jhan: ',com, var
   print *, '' 
end subroutine jhprint_r

subroutine jhprint_2r(com, var)
   chaRACTER(len=*) :: com
   real :: var(:)
   print *, '' 
   print *, 'jhan: ',com, var
   print *, '' 
end subroutine jhprint_2r

subroutine jhprint_2rr(com, var)
   chaRACTER(len=*) :: com
   real :: var(:,:)
   print *, '' 
   print *, 'jhan: ',com, var
   print *, '' 
end subroutine jhprint_2rr







!   subroutine cable_switch_status(diag_message, diag_var )
!      implicit none
!      character(len=*), intent(in) :: diag_message, diag_var 
!      character(len=68) :: message
!         message=trim(trim('CABLE_log:')//trim(diag_message))            &
!                  //' '//trim(diag_var)
!         if(ktau_gl==1 .and. knode_gl==0 ) then
!            print *, message
!         endif   
!      return
!   end subroutine cable_switch_status
!   subroutine calc_rhoch(veg) 
!      use define_types
!      use other_constants
!      implicit none
!      type (veg_parameter_type), intent(inout) :: veg
!
!      allocate( c1(mp,nrb), rhoch(mp,nrb) )
!
!!jhan:UM uses rad%extkn instead of veg%extkn, which should be read from par. file anyway
!!jhan:cahnge Mk3l to read veg%taul like UM 
!#if defined(CABLE_17TILES) || defined(CABLE_9TILES) 
!#else
!         veg%taul(:,1) = taul(1)
!         veg%taul(:,2) = taul(2)
!         veg%refl(:,1) = refl(1) 
!         veg%refl(:,2) = refl(2) 
!#endif                  
!         c1(:,1) = SQRT(1. - veg%taul(:,1) - veg%refl(:,1))
!         c1(:,2) = SQRT(1. - veg%taul(:,2) - veg%refl(:,2))
!         c1(:,3) = 1.
!          
!         ! Canopy reflection black horiz leaves (eq. 6.19 in Goudriaan and van Laar, 1994):
!         rhoch = (1.0 - c1) / (1.0 + c1)
!      return
!   end subroutine calc_rhoch 

end module cable_common_module

