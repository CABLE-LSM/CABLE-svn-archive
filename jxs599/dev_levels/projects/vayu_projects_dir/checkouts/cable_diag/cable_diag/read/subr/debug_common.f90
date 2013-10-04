
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_common
   implicit none
   !--- defined in debug_directives.f90
   integer, parameter :: i_d = DEFi_d
   integer, parameter :: r_1 = DEFr_1
   !--- base of filename read by subr. read_args from file input.dat
   !--- which is created by script cable_diag.pl after reading files created
   !--- host program (i.e. UM, Mk3L, offline CABLE)
   character(len=30) :: filename 
   !--- t_window = averaging window for time series
   !--- x_window = variance of Gaussian filter which data is convoluted with
   !--- in units of # cells 
   integer :: t_window, x_window
   !--- do you want to write/plot ?
   logical :: write_flag, plot_flag
   !--- vars read by subr. read_args -also  from file input.dat
   !--- Nvars = # vars contained in binary file output by host program 
   !--- dimx = typically # landpoints over which the var is specified at each timestep 
   !--- dimy = # timesteps
   integer(i_d) :: Nvars, dimx, dimy 
   !--- array of the varibale names, appears in plot/written text file
   character(len=30), dimension(:), allocatable :: ar_varname
   !--- array to read in the 1D vector of Nvars * dimx  
   real(r_1), dimension(:), allocatable :: ar_Nvars
   !--- array which seperates Nvars in slices of the 3D (cubic) array
   !--- each row contains dimx values, dimy columns NB. this may be incorrect wrt to std 
   !--- fortran practice, but you get the idea 
   real(r_1), dimension(:,:,:), allocatable :: ar_data
end module debug_common

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

