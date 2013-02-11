
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ potentially adaptable to anything however current usage for diagnostic+++!
!+++ of global climate models which operate over many landpoints, producing+++!
!+++ a vector(len=landpoints) of values at each timestep. general motiv.   +++!    
!+++ is to be able to output desired vars from GCM without expensive print +++!
!+++ to std out, which then dumps messily into single file. the tools (of  +++!
!+++ which this program is a part) allows the user to dump any var. from   +++!
!+++ the GCM(host) in inexpensive binary format, inito a desired location, +++!
!+++ by including cable_diag.f90 in their build, including this module in  +++!
!+++ subr. from which they want to output vars (i.e. 'use cable_diag_mod'))+++!
!+++ , and then include statement in subr. 'call cable_diag(...args...)    +++!
!+++ [see cable_diag.f90 for details]. binary data can then be cleanly     +++!
!+++ retrieved by calling this program from the command line, specifying   +++!
!+++ desired behaviour through given args. see cable_diag.pl for details.  +++!   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

#include "debug_directives.f90" !fpp directives

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program debug 
   use debug_common, only : filename1, filename2, write_flag, plot_flag,      &
         Nvars,dimy, olddata,newdata     
   use debug_read_mod   !only called funcs are public in this mod 
!   use debug_write_mod  !only called funcs are public in this mod
   implicit none
   integer :: t,j       !primarily counters, double as args also
      !======================================================================!
      !=== read perl script interp. (input.dat) of command line args      ===!
      !--- which determine behaviour of program. which file to process,   ===!
      !=== plot/write text file, how to smooth the data                   ===! 
      !======================================================================!
      call read_args
      
      !======================================================================!
      !=== read info about the spec. binary data which was created by the ===!
      !--- host so we know how many vars are contained within, how many   ===!
      !=== points there are at each timestep, how many timesteps.         ===! 
!      if(write_flag .or. plot_flag) then 
      !======================================================================!
         call read_txt_file( trim(filename1) )
         call read_dat_file( trim(filename1),olddata )
         call read_dat_file( trim(filename2),newdata )
!      endif
#define read_diff_only
#ifdef read_diff_only
         call comp_diff( olddata, newdata )
      stop
#endif
!      !======================================================================!
!      !=== if command line arg plot_flag==.true., for each var. in binary ===! 
!      !=== "filename".bin, w x points smoothed by x_window Gaussian       ===!
!      !=== debug_dat(x,y) x=# of the var. in data,y=total # timesteps     ===!
!      !======================================================================!
!#ifdef testxs
!      !======================================================================!
!      !=== for testing purposes plot first timestep to screen             ===!  
!      !======================================================================!
!      j=1;t=1
!      if(plot_flag) call debug_dat(j, t ) 
!#else
!      !======================================================================!
!      !=== else for timesteps averaged by t_window.                       ===!  
!      !======================================================================!
!      do j=1,Nvars
!         if(plot_flag) call debug_dat(j, dimy ) 
!      enddo
!#endif
!      !======================================================================!
!      !=== write text file "unsmoothed" if write_flag==.true.             ===!  
!      !======================================================================!
!      if(write_flag) call write_txt_file( trim(filename, old_data) )
!   stop
end program debug 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!




!=============================================================================!
!=== subr. to read command line args interpreted by perl script.comm. line ===!
!=== see top description of program for further explanation                ===! 
!=============================================================================!
   
subroutine read_args
   use debug_common
   implicit none
   integer(i_d), parameter :: gok=0
   integer(i_d) :: gopenstatus
      open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
         if(gopenstatus==gok) then
            read(1,*), filename1
            read(1,*), filename2
            !read(1,*), write_flag
            !read(1,*), plot_flag
            !read(1,*), t_window
            !read(1,*), x_window
         else
            write (*,*), 'input.dat',' NOT found to read'
         endif
      close(1)
   return
end subroutine read_args
   
!=============================================================================!
!=== subr. manages desired 'smoothing' over landpoints and timesteps, and  ===!
!=== calls more generic subrs. to compute this data from read binary       ===! 
!=============================================================================!

!subroutine debug_dat( var_j, n_tsteps )
!   use debug_common
!!   use debug_comp_mod
!!   use debug_plot_mod
!   implicit none
!   integer,intent(in) :: var_j      ! # of var. in data
!   integer,intent(in) :: n_tsteps   ! total # timesteps to consider
!   integer, save :: n_call          ! how many times as we been here
!   integer ::i, t                   ! counters
!      n_call =  1
!      !---  average over t_window
!      if( t_window > 1 ) then 
!         do i=1,n_tsteps, t_window
!            ary_twindow = real(0)   ! initialize each window in series    
!            do t=i,(i-1)+t_window
!               !---  compute data per time_step in t_window
!               call comp_data( var_j, t )
!               ary_twindow = ary_twindow + ar_y    ! sum within window
!            enddo
!            ar_y = ary_twindow / real(t_window)    ! average in window
!            !--- smooth in x direction
!            if(x_window>1) then
!               call smooth_av( dimx, ar_x, ar_y, x_window ) 
!            endif
!            !--- plot each effective timestep
!            call plot_data( var_j, n_call )
!            n_call = n_call + 1
!         enddo
!      else
!         do i=1,n_tsteps
!            !---  compute data per time_step 
!            call comp_data( var_j, i )
!            !--- smooth in x direction
!            if(x_window>1) then
!               call smooth_av( dimx, ar_x, ar_y, x_window ) 
!            endif
!            !--- plot each effective timestep
!            call plot_data( var_j, n_call )
!            n_call = n_call + 1
!         enddo
!      endif
!   return   
!end subroutine debug_dat
      
subroutine comp_diff( )
   use debug_common
   real(r_1) :: sum_data
      diff_data = olddata-newdata
      sum_diff = sum(diff_data)
      print *, ''
      print *, 'summed difference between old and new data :sum_diff'
      print *, sum_diff
      print *, ''
   return
end subroutine comp_diff
!=======================================================================!
!=======================================================================!






