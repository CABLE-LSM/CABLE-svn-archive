!!this code illustrates the dynamics of using arrays in another function which have been declared in a module outside of the function  - pointers need to be used

module generic_pars
   implicit none
   integer,parameter :: gs=4  !default 4 byte kind !!compiler specific value
   integer,parameter :: gd=8  !default 8 byte kind !!compiler specific value
   integer(gs),parameter :: gchlen=99 !default max character length
   integer(gs), parameter :: ok=0
   integer(gs) :: openstatus,alloctest=1
   integer(gs), parameter :: deriv_method=0  !0 => derivative constant at end; 1 => constant at beginning
   integer(gs), parameter :: HOLES=-1  !0 => derivative constant at end; 1 => constant at beginning
   integer(gs), parameter :: POLES=1  !0 => derivative constant at end; 1 => constant at beginning
   real(gs), parameter :: maxvalue=10000.0
   real(gs), parameter :: gminvalue=-10000.0
end module generic_pars

module var_mod 
   use generic_pars,  only : gs, alloctest
   implicit none

   type var_type
      real, dimension(10) :: arv
   end type var_type

   type var_typeb 
      real, dimension(10) :: arv
   end type var_typeb

   type (var_type) :: var
   type (var_typeb) :: var
!      contains
!      subroutine falloc 
!         implicit none
!            allocate( ar(mp), stat=alloctest ) !this works!!
!         return
!      end subroutine falloc
end module var_mod


