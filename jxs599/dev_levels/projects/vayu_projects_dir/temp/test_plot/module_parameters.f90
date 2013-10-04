
module generic_pars
   implicit none
   integer,parameter :: gs=4  !default 4 byte kind !!compiler specific value
   integer,parameter :: gd=8  !default 8 byte kind !!compiler specific value
   integer(gs),parameter :: gchlen=99 !default max character length
   integer(gs), parameter :: gok=0
   integer(gs) :: gopenstatus,galloctest=1
   integer(gs), parameter :: gderiv_method=0 !0=>deriv=constant at end; 1=>constant at beginning
   integer(gs), parameter :: gHOLES=-1  
   integer(gs), parameter :: gPOLES=1  
   real(gs), parameter :: gpi= 3.14
   real(gs), parameter :: gmaxvalue=10000.0
   real(gs), parameter :: gminvalue=-10000.0
end module generic_pars



