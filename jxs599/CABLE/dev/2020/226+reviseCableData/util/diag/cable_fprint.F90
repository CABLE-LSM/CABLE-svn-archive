!H!MODULE cable_fprint_module
!H!
!H!  interface cable_fprintf
!H!    module procedure cable_fprintf0, cable_fprintf1, cable_fprintf2
!H!  end interface cable_fprintf
!H!
!H!CONTAINS
!H!
!H!! Writes to std output stereamm(=6) 
!H!!===============================================================================
!H!
!H!SUBROUTINE cable_fprintf0( basename, L_fprint )
!H!  USE cable_common_module, only : knode_gl, ktau_gl
!H!  implicit none  
!H!  
!H!  character(len=*) :: basename     ! subr name 
!H!  logical :: L_fprint
!H!
!H!  if( .NOT. L_fprint ) return
!H!   
!H!  call cable_print0( 6, basename, ktau_gl, knode_gl )
!H!                             
!H!END SUBROUTINE cable_fprintf0
!H!
!H!SUBROUTINE cable_print0(idiag, basename, ktau, mype )
!H!  implicit none  
!H!  integer :: idiag
!H!  character(len=*) :: basename     ! subr name 
!H!  integer :: ktau, mype 
!H!
!H!19 format(  "CABLE_LSM: ", A20, " Finished."   )
!H!  write (idiag, 19)  basename
!H!
!H!END SUBROUTINE cable_print0
!H!
!H!!===============================================================================
!H!
!H!! Writes to file
!H!!===============================================================================
!H!
!H!SUBROUTINE cable_fprintf1( iDiag, basename, knode, ktau, L_fprint )
!H!  USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
!H!  USE cable_common_module, only : 
!H!  implicit none  
!H!
!H!  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
!H!  integer :: ktau, knode
!H!  ! writes file per processor (basename+node)
!H!  character(len=*) :: basename     ! subr name 
!H!  logical :: L_fprint
!H!  ! LOCAL vars
!H!  integer, SAVE :: pDiag=3713      ! give unique SEED per module procedure 
!H!
!H!  if( .NOT. L_fprint ) return
!H!
!H!  ! Returns unique unit=iDiag and modified basename
!H!  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode )
!H!   
!H!  call cable_print1( idiag, basename, ktau, knode )
!H!                             
!H!END SUBROUTINE cable_fprintf1
!H!
!H!SUBROUTINE cable_print1(idiag, basename, ktau, mype )
!H!  implicit none  
!H!  integer :: idiag
!H!  character(len=*) :: basename     ! subr name 
!H!  integer :: ktau, mype 
!H!
!H!19 format(  "CABLE_LSM: ", A20, " finished @", I8.1, " on ",I3.1  )
!H!  write (idiag, 19)  basename, ktau, mype 
!H!
!H!END SUBROUTINE cable_print1
!H!
!H!!===============================================================================
!H!
!H!! Writes to file
!H!!===============================================================================
!H!SUBROUTINE cable_fprintf2( iDiag, basename, msg, L_fprint )
!H!  USE cable_fFile_module, ONLY : open_file_per_node, fprintf_dir
!H!  USE cable_common_module, only : knode_gl, ktau_gl
!H!  implicit none  
!H!  
!H!  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
!H!  ! writes file per processor (basename+node)
!H!  character(len=*) :: basename     ! subr name 
!H!  character(len=*) :: msg ! message 
!H!  logical :: L_fprint
!H!  ! LOCAL vars
!H!  integer, SAVE :: pDiag=5713      ! give unique SEED per module procedure 
!H!
!H!  if( .NOT. L_fprint ) return
!H!
!H!  ! Returns unique unit=iDiag and modified basename
!H!  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
!H!   
!H!  call cable_print2( idiag, basename, msg, ktau_gl, knode_gl )
!H!                             
!H!END SUBROUTINE cable_fprintf2
!H!
!H!SUBROUTINE cable_print2(idiag, basename, msg, ktau, mype )
!H!  implicit none  
!H!  integer :: idiag
!H!  character(len=*) :: basename     ! subr name 
!H!  character(len=*) :: msg ! message 
!H!  integer :: ktau, mype 
!H!
!H!  write (idiag, *) ""  
!H!  write (idiag, *) msg  
!H!  write (idiag, *) &
!H!    "lsm_id set CABLE. call succeeded however return and call JULES" 
!H!  write (idiag, *) ""  
!H!
!H!END SUBROUTINE cable_print2
!H!
!H!
!H!
!H!END MODULE cable_fprint_module
!H!
!H!
