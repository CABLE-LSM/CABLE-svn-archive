
#include "debug_directives.f90"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program debug 
   use debug_common,   only : filename1, filename2, filename3, filename4, Nvars
   use debug_read_mod, only : read_args, read_txt_file, read_dat_file
   use debug_comp_mod, only : x_data, al_data, comp_data, de_data
   use debug_plot_mod, only : plot_data
   implicit none
   integer :: i
      call read_args
      i=0 
      call read_txt_file( trim(filename1) )
      call x_data()
      call al_data()
      call read_dat_file( trim(filename1) )
      call comp_data( i )
      i = i + Nvars
      call de_data()

      call read_txt_file( trim(filename2) )
      call al_data()
      call read_dat_file( trim(filename2) )
      call comp_data( i )
      i = i + Nvars
      call de_data()

      call read_txt_file( trim(filename3) )
      call al_data()
      call read_dat_file( trim(filename3) )
      call comp_data( i )
      i = i + Nvars
      call de_data()

      call read_txt_file( trim(filename4) )
      call al_data()
      call read_dat_file( trim(filename4) )
      call comp_data( i )
      i = i + Nvars
      call de_data()
 
      call plot_data()
   stop
end program debug 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !=======================================================================!
      !=======================================================================!
                                                                                

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!





