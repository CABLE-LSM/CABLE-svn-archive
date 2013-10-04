
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


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program debug 
   use debug_common, only : filename, write_flag, plot_flag, Nvars,dimy     
   use debug_read_mod   !only called funcs are public in this mod 
   use debug_write_mod  !only called funcs are public in this mod
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
      if(write_flag .or. plot_flag) then 
      !======================================================================!
         call read_txt_file( trim(filename) )
         call read_dat_file( trim(filename) )
      endif

      !======================================================================!
      !=== write text file "unsmoothed" if write_flag==.true.             ===!  
      !======================================================================!
      if(write_flag) call write_txt_file( trim(filename) )
   stop
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
   integer, parameter :: gok=0
   integer :: gopenstatus
      open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
         if(gopenstatus==gok) then
            read(1,*), filename
            read(1,*), write_flag
            read(1,*), plot_flag
            read(1,*), t_window
            read(1,*), x_window
         else
            write (*,*), 'input.dat',' NOT found to read'
         endif
      close(1)
   return
end subroutine read_args
   
 
!=======================================================================!
!=======================================================================!






