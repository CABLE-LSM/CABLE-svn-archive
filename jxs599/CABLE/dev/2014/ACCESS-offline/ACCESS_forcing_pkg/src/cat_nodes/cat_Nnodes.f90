
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program  catter
   use cat_common, only : filename, n_nodes, i_nodes, logu     
   use cat_nodes_mod   
   implicit none
   character(len=*), parameter :: logfile = '/home/599/jxs599/catter.log'
      ! get CLI args, basename of file & how many nodes
      IF( IARGC() > 0 ) THEN
         CALL GETARG(1, filename)
         CALL GETARG(2, n_nodes)
      ENDIF
      logu=444
      open( unit=logu,file=logfile, status="unknown",action="write", &
            position='append' )
         write(logu,*) ""
         write(logu,*) ""
         write(logu,*) "Fortran executable: catter - concatenating multiples files"
         write(logu,*) "called with args: " 
         write(logu,*) "filename: ", trim(filename)
         write(logu,*) "N nodes: ", trim(n_nodes)
         write(logu,*) "executing.... " 

      ! reformat #nodes for fortran 
      READ(n_nodes,12) i_nodes
   12 FORMAT(I3.3)
      
      ! call driver of "model" components
      call cat_nodes() 
      
         write(logu,*) "Fortran executable: catter"
         write(logu,*) "Finished." 
      
      close(unit=logu)
   stop
end program catter 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!



