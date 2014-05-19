
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

program  catter
   use cat_common, only : filename, n_nodes, i_nodes     
   use cat_nodes_mod   
   implicit none

      ! get CLI args, basename of file & how many nodes
      IF( IARGC() > 0 ) THEN
         CALL GETARG(1, filename)
         CALL GETARG(2, n_nodes)
      ENDIF
      
      ! reformat #nodes for fortran 
      READ(n_nodes,12) i_nodes
   12 FORMAT(I3.3)
      
      ! call driver of "model" components
      call cat_nodes() 

   stop
end program catter 
      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!



