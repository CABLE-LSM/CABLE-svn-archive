!H!MODULE cable_Pyfprint_module
!H!  USE cable_fFile_module
!H!
!H!   interface cable_Pyfprintf
!H!      module procedure cable_Pyfprintf1, cable_Pyfprintf2
!H!   end interface cable_Pyfprintf
!H!
!H!CONTAINS
!H!
!H!! 1-D REAL
!H!!==========================================================================!
!H!
!H!SUBROUTINE cable_Pyfprintf1( iDiag, basename, var1, dimx, L_fprint )
!H!  USE cable_fFile_module
!H!  USE cable_common_module, only : knode_gl, ktau_gl
!H!  implicit none  
!H!  ! IN vars
!H!  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
!H!  integer :: dimx   ! 1-D length
!H!  real, dimension(dimx) :: var1    ! var CALLed    
!H!  ! writes file per processor (basename+node)
!H!  character(len=*) :: basename     ! filename based on var
!H!  logical :: L_fprint
!H!  ! LOCAL vars
!H!  integer, SAVE :: pDiag=1713      ! give unique SEED per module procedure 
!H!
!H!  if( .NOT. L_fprint ) return
!H!
!H!  ! Returns unique unit=iDiag and modified basename
!H!  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
!H!   
!H!  call cable_Pyprint1( idiag, dimx, var1, ktau_gl )
!H!                             
!H!END SUBROUTINE cable_Pyfprintf1
!H!
!H!
!H!SUBROUTINE cable_Pyprint1(idiag,dimx,field, ktau)
!H!  implicit none  
!H!  integer :: idiag
!H!  integer :: dimx
!H!  real, dimension(dimx) :: field
!H!  integer :: ktau 
!H!  
!H!  integer :: j
!H!
!H!  do j=1, dimx  
!H!      write (iDiag,*) field(j) 
!H!  enddo
!H!
!H!END SUBROUTINE cable_Pyprint1
!H!
!H!! 2-D REAL
!H!!==========================================================================!
!H!
!H!SUBROUTINE cable_Pyfprintf2( iDiag, basename, var1, dimx, dimy, L_fprint )
!H!  USE cable_fFile_module
!H!  USE cable_common_module, only : knode_gl, ktau_gl
!H!  implicit none  
!H!  ! IN vars
!H!  integer :: iDiag  ! f^n creates unique unidID to be returned to calling point
!H!  integer :: dimx, dimy   ! 2-D length
!H!  real, dimension(dimx,dimy) :: var1    ! var CALLed    
!H!  ! writes file per processor (basename+node)
!H!  character(len=*) :: basename     ! filename based on var
!H!  logical :: L_fprint
!H!  ! LOCAL vars
!H!  integer, SAVE :: pDiag=2713      ! give unique SEED per module procedure 
!H!
!H!  if( .NOT. L_fprint ) return
!H!
!H!  ! Returns unique unit=iDiag and modified basename
!H!  call open_file_per_node( iDiag, pDiag, fprintf_dir, basename, knode_gl )
!H!   
!H!  call cable_Pyprint2( iDiag, dimx, dimy, var1, ktau_gl )
!H!                             
!H!END SUBROUTINE cable_Pyfprintf2
!H!
!H!
!H!SUBROUTINE cable_Pyprint2(idiag,dimx,dimy, field, ktau)
!H!  implicit none  
!H!  integer :: idiag
!H!  integer :: dimx, dimy
!H!  real, dimension(dimx, dimy) :: field
!H!  integer :: ktau 
!H!  
!H!  integer :: i,j
!H!
!H!  do i=1, dimx  
!H!  do j=1, dimy  
!H!      write (iDiag,*) field(i,j) 
!H!  enddo
!H!  enddo
!H!
!H!END SUBROUTINE cable_Pyprint2
!H!
!H!END MODULE cable_Pyfprint_module
!H!
!H!
