!H!MODULE cable_fFile_module
!H!  use cable_common_module
!H!  
!H!  IMPLICIT NONE
!H!  INTEGER, PARAMETER :: gok=0
!H!  INTEGER :: galloctest=1
!H!
!H!  !CABLE_LSM: intro'd quick writing capbility. remove from here. keep for ref
!H!  character(len=200) :: fprintf_dir_root = "/"
!H!  
!H!  character(len=15) :: unique_subdir = "./"
!H!
!H!  logical :: L_cable_fprint   = .FALSE.,    &  
!H!             L_cable_Pyfprint = .FALSE. 
!H!              
!H!  character(len=300) :: fprintf_dir
!H!
!H!CONTAINS
!H!
!H!SUBROUTINE open_file_per_node( iDiag,pDiag, dir, basename, node, fbasename )
!H!  ! use cable_common_module
!H!   integer :: iDiag, pDiag, node
!H!   integer :: gopenstatus = 1
!H!   character(len=*)  :: dir 
!H!   character(len=*)  :: basename
!H!   character(len=*), optional  :: fbasename
!H!   character(len=300) :: infilename
!H!   character(len=30) :: chnode
!H!
!H!   write(chnode,10) node
!H!10 format(i3.3)   
!H!   infilename=trim( trim(dir)//trim(basename)//trim(chnode) )
!H!    
!H!   IF(iDiag==0) tHEN
!H!      pDiag = pDiag+2  
!H!      iDiag=pDiag
!H!         
!H!      call open_iDiag( iDiag, infilename, gopenstatus) 
!H!      
!H!   ENDIF
!H!   
!H!   !jhan:check if file is open      
!H!   if(gopenstatus==gok) then
!H!      fbasename = infilename
!H!      return
!H!   else
!H!      write (*,*) infilename,' NOT open for write. Error(open_file_per_node)'
!H!      write (*,*) '???Check that the path exists???'
!H!      STOP
!H!   endif
!H!
!H!END SUBROUTINE open_file_per_node
!H!
!H!subroutine qprint( iDiag, infilename) 
!H!
!H!!   use cable_common_module
!H!   integer :: iDiag 
!H!   character(len=*) :: infilename
!H!   character(len=300) :: ffilename
!H!   
!H!   ffilename=trim( trim(infilename)// '.txt' )
!H!
!H!   open( unit=iDiag, file=ffilename, status="unknown", &
!H!     action="write", form="formatted", &
!H!     position='append' )
!H!
!H!End subroutine qprint 
!H!
!H!
!H!subroutine open_iDiag( iDiag, infilename, gopenstatus) 
!H!
!H!!   use cable_common_module
!H!   integer :: iDiag 
!H!   integer :: gopenstatus
!H!   character(len=*) :: infilename
!H!   character(len=300) :: ffilename
!H!   
!H!   ffilename=trim( trim(infilename)// '.txt' )
!H!
!H!   open( unit=iDiag, file=ffilename, status="unknown", &
!H!     action="write", iostat=gopenstatus, form="formatted", &
!H!     position='append' )
!H!
!H!End subroutine open_iDiag
!H!
!H!END MODULE cable_fFile_module
!H!
!H!
