
MODULE cable_fprint_module
   IMPLICIT NONE
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
 
   interface cable_fprintf
      module procedure cable_fprintf1!, cable_fprintf2
   end interface cable_fprintf
   
CONTAINS

! writes text files 
SUBROUTINE cable_fprintf1( iDiag, basename, dimx, dimy, timestep, node, &
                        vname1, var1 )
   integer, intent(inOUT) :: iDiag 
   integer, SAVE :: pDiag=713 
   integer, intent(in) :: dimx, dimy, timestep,node
   integer, save :: gopenstatus = 1
   integer, intent(in) :: var1
   !integer, intent(in), dimension(:) :: var1
   !real, intent(in), dimension(:) :: var1
   integer :: i=0
   character(len=*), intent(in) :: basename, vname1
   character(len=6) :: filename
   character(len=30) :: chnode
 
      write(chnode,10) node
   10 format(i2.2)   
      filename=trim(trim(basename)//trim(chnode))
    
      IF(iDiag==0) tHEN
         pDiag = pDiag+2  
         iDiag=pDiag
         
         open(unit=iDiag+1,file=trim(filename//'.txt'),status="unknown", &
           action="write", iostat=gopenstatus, form="formatted", &
           position='append' )
      
      ENDIF
         
   if(gopenstatus==gok) then
         write (iDiag+1,*) trim(vname1), var1
   else
      write (*,*) filename//'.txt',' NOT open for write. Error'
   endif

   if (timestep == dimy) & 
      close(iDiag+1)
     
                             
END SUBROUTINE cable_fprintf1


END MODULE cable_fprint_module
 
!==========================================================================!


