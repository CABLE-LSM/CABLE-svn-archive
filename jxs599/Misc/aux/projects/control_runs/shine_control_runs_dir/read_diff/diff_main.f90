
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program debug 
   use debug_common, only : filename1, filename2, write_flag, plot_flag,      &
         Nvars,dimy, olddata,newdata     
   use debug_read_mod   !only called funcs are public in this mod 
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

         call comp_diff( )

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
   
   open(unit=1,file='input.dat', status="unknown",action="read", &
            iostat=gopenstatus )
      if(gopenstatus==gok) then
         read(1,*), filename1
         read(1,*), filename2
      else
         write (*,*), 'input.dat',' NOT found to read'
      endif
   close(1)
   return
end subroutine read_args
   
     
subroutine comp_diff( )
   use debug_common
   real :: sum_data
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






