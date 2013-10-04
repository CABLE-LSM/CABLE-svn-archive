!#test diff
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

module debug_common
   implicit none
   integer, parameter :: i_d = DEFi_d
   integer, parameter :: r_1 = DEFr_1
   character(len=30) :: filename1, filename2, filename3, filename4 
   integer :: t_window, x_window
   logical :: write_flag, plot_flag
   integer(i_d) :: Nvars, dimx, dimy 
   character(len=30), dimension(:), allocatable :: ar_varname
   real(r_1), dimension(:), allocatable :: ar_Nvars
   real(r_1), dimension(:,:,:), allocatable :: ar_data
   real, dimension(:,:), allocatable :: cmp_data
end module debug_common

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

