
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

module data_pars 
   use generic_pars
   implicit none
!   real(gs), parameter :: correction_factor = (0.38/0.7)
   integer(gs),parameter :: gNfiles=2  !number of mfps=num. of files !!dummy
   integer(gs),parameter :: gNFobs=3  !number of mfps=num. of files !!dummy
   integer(gs),parameter :: gNx_i=8  !number of mfps=num
   integer(gs),parameter :: gNy_i=100  !number of Gammas at each mfp
   integer(gs),parameter :: gNx_f=100  !number of mfps to extrapolate to 
   integer(gs),parameter :: gNy_f=100  !number of Gammas at each mfp after extrapolation
   character (len=gchlen), dimension(gNfiles), parameter :: datafiles = (/"Input/grid.dat.0.001","Input/grid.dat.0.005"/)
   real(gs), parameter :: gGJ21_propto=2.82  
   real(gs), parameter, dimension(gNFobs) :: gFobs = (/ 0.080, 0.040, 0.005 /)  
   real(gs), parameter, dimension(gNFobs) :: gFobserr = (/ 0.007, 0.040, 0.005 /)  
end module data_pars



