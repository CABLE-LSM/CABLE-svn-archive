
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program debug 
   use debug_common, only : filename, n_nodes     
   use cat_nodes_mod   
   implicit none
      !======================================================================!
      !=== read perl script interp. (input.dat) of command line args      ===!
      !--- which determine behaviour of program. which file to process,   ===!
      !=== plot/write text file, how to smooth the data                   ===! 
      !======================================================================!
      call read_args

      !======================================================================!
      !======================================================================!
      call cat_nodes() 

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
   integer(i_d), parameter :: gok=0
   integer(i_d) :: gopenstatus
      open(unit=1,file='input.dat', status="unknown",action="read", iostat=gopenstatus )
         if(gopenstatus==gok) then
            read(1,*), filename
            read(1,*), n_nodes 
         else
            write (*,*), 'input.dat',' NOT found to read'
         endif
      close(1)
   return
end subroutine read_args

