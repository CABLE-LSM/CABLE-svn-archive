
module cat_common
   implicit none
   
   integer :: logu ! unit number for log file
   !--- base of filename read by subr. read_args from file input.dat
   !--- which is created by script cable_diag.pl after reading files created
   !--- host program (i.e. UM, Mk3L, offline CABLE)
   character(len=2000) :: filename, n_nodes 
    
   !--- t_window = averaging window for time series
   !--- x_window = variance of Gaussian filter which data is convoluted with
   !--- in units of # cells 
   integer :: t_window, x_window
   
   !--- do you want to write/plot ?
   logical :: write_flag, plot_flag
   
   !--- vars read by subr. read_args -also  from file input.dat
   !--- Nvars = # vars contained in binary file output by host program 
   !--- dimx = typically # landpoints over which the var is specified at 
   !--- each timestep 
   !--- dimy = # timesteps
   integer :: Nvars, dimx, dimy, i_nodes 
   integer :: dimyQ1, dimyQ2, dimyQ3,dimyQ4

   !--- array of the varibale names, appears in plot/written text file
   character(len=30), dimension(:), allocatable :: ar_varname
   
   !--- array to read in the 1D vector of Nvars * dimx  
   real, dimension(:), allocatable :: ar_Nvars
   
   !--- array which seperates Nvars in slices of the 3D (cubic) array
   !--- each row contains dimx values, dimy columns NB. this may be incorrect 
   !--- wrt to std fortran practice, but you get the idea 
   real, dimension(:,:,:), allocatable :: ar_data
   integer, dimension(:), allocatable :: dimx_i
   integer, save :: dimx_tot=0 

end module cat_common

