!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: handles additional, dynamically decided diagnostic output from model.
!          permanently used for bitwise identical testing. more applications 
!          will follow.   
!
! Contact: Jhan.Srbinovsky@csiro.au
!
! History: Currently stripped down version of cable_diag here. will be 
!          re-implemented in time.
!
! ==============================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ x is typically the number of landpoints(tiles). binary file is        +++!
!+++ then appended every timestep with the new foo(x_i)                    +++!
!+++                                                                       +++! 
!+++ CALL syntax:                                                          +++!  
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++! 
!+++ output binaries can be interpreted from the command line              +++!
!+++ using a suite of tools. Currently, only zero_diff.ksh is supported.   +++!  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


MODULE cable_diag_module
   IMPLICIT NONE
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
  
   !--- subrs overloaded to respond to call cable_diag 
   INTERFACE cable_diag
      MODULE PROCEDURE cable_diag1
   END INTERFACE cable_diag
  
CONTAINS

!==========================================================================!
! cable_diag1/2/3 call subrs to write filename.dat which contains description
! of data and format etc., and filename.bin containing the data   
!==========================================================================!

SUBROUTINE cable_diag1( iDiag, basename, dimx, dimy, timestep, node, &
                        vname1, var1 )
   integer, intent(inOUT) :: iDiag 
   integer, SAVE :: pDiag=713 
   integer, intent(in) :: dimx, dimy, timestep,node
   real, intent(in), dimension(:) :: var1
   integer :: Nvars=1 !this WAS input
   integer :: i=0
   character(len=*), intent(in) :: basename, vname1
   character(len=30) :: filename, chnode
 
      IF(iDiag==0) tHEN
         pDiag = pDiag+2  
         iDiag=pDiag
      ENDIF
         
      write(chnode,10) node
   10 format(i3.3)   
      filename=trim(trim(basename)//trim(chnode))
      
      if (timestep == 1) & 
         call cable_diag_desc1( iDiag, trim(filename), dimx, dimy, vname1 )
      
      call cable_diag_data1( iDiag, trim(filename), dimx, timestep, dimy, &
                             var1 )
END SUBROUTINE cable_diag1

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc1( iDiag, filename, dimx, dimy, vname1 )

   integer, intent(in) :: iDiag,dimx,dimy 
   integer, PARAMETER :: Nvars=1
   character(len=*), intent(in) :: filename, vname1
   integer, save :: gopenstatus = 1

     open(unit=iDiag,file=filename//'.dat', status="replace", &
          action="write", iostat=gopenstatus )
     
      if(gopenstatus==gok) then
            write (iDiag,*) 'Number of var(s): '
            write (iDiag,*) Nvars
            write (iDiag,*) 'Name of var(s): '
            write (iDiag,7139) vname1 
 7139       format(a)            
            write (iDiag,*) 'dimension of var(s) in x: '
            write (iDiag,*) dimx 
            write (iDiag,*) 'dimension of var(s) in y: '
            write (iDiag,*) dimy 
      else
         write (*,*), filename//'.dat',' Error: unable to write'
      endif
      
   close(iDiag)
  
END SUBROUTINE cable_diag_desc1


SUBROUTINE cable_diag_data1( iDiag, filename, dimx, timestep, kend, var1  )

   integer, intent(in) :: iDiag, dimx, timestep, kend
   integer, PARAMETER :: Nvars=1
   real, intent(in), dimension(:) :: var1
   character(len=*), intent(in) :: filename
   integer, save :: gopenstatus = 1

   if (timestep == 1)  then 
      open(unit=iDiag+1,file=filename//'.bin',status="unknown", &
           action="write", iostat=gopenstatus, form="unformatted", &
           position='append' )
   endif   
 
   if(gopenstatus==gok) then
         write (iDiag+1) var1
   else
      write (*,*) filename//'.bin',' NOT open for write. Error'
   endif

   if (timestep == kend) & 
      close(iDiag+1)

END SUBROUTINE cable_diag_data1

!=============================================================================!
!=============================================================================!

!RL: add cable_diag2 for 3d array (tiles on proc, sm_level, timesteps)

SUBROUTINE cable_diag2( iDiag, basename, dimx, kend, dimz, timestep, &
                        node, vname2, var2, start_run )
   integer, intent(inOUT) :: iDiag 
   integer, SAVE :: pDiag=713 
   integer, intent(in) :: dimx, kend, dimz, timestep, node
   real, intent(in), dimension(:,:) :: var2
   logical, intent(in) :: start_run
   integer :: Nvars=1 !this WAS input
   integer :: i=0
   integer :: dimy
   character(len=*), intent(in) :: basename, vname2
   character(len=30) :: filename, chnode
 
   dimy = kend-timestep+1

      IF(iDiag==0) THEN
         pDiag = pDiag+2  
         iDiag=pDiag
      ENDIF
         
      write(chnode,10) node
   10 format(i3.3)   
      filename=trim(trim(basename)//trim(chnode))
      
      if (start_run) & 
      call cable_diag_desc2( iDiag, trim(filename), dimx, &
      	   		      dimy, dimz, vname2, timestep )
      
      call cable_diag_data2( iDiag, trim(filename),dimx, dimz,  &
                              timestep, kend, var2, start_run )
END SUBROUTINE cable_diag2

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc2( iDiag, filename, dimx, &
	   		     dimy, dimz, vname2, timestep )

   integer, intent(in) :: iDiag, dimx, dimy, dimz, timestep
   integer, PARAMETER :: Nvars=1
   character(len=*), intent(in) :: filename, vname2
   integer, save :: gopenstatus = 1

     open(unit=iDiag,file=filename//'.dat', status="replace", &
          action="write", iostat=gopenstatus )
     
      if(gopenstatus==gok) then
            write (iDiag,*) 'Number of var(s): '
            write (iDiag,*) Nvars
            write (iDiag,*) 'Name of var(s): '
            write (iDiag,7139) vname2 
 7139       format(a)   
            write (iDiag,*) 'Timestep is: '
            write (iDiag,*)  timestep        
            write (iDiag,*) 'dimension of var(s) in x: '
            write (iDiag,*) dimx 
            write (iDiag,*) 'dimension of var(s) in y: '
            write (iDiag,*) dimy
            write (iDiag,*) 'dimension of var(s) in z: '
            write (iDiag,*) dimz 
      else
         write (*,*), filename//'.dat',' Error: unable to write'
      endif
      
   close(iDiag)
  
END SUBROUTINE cable_diag_desc2


SUBROUTINE cable_diag_data2( iDiag, filename, dimx, dimz,	&
	   		      timestep, kend, var2, start_run  )

   integer, intent(in) :: iDiag, dimx, dimz, timestep, kend
   integer, PARAMETER :: Nvars=1
   real, intent(in), dimension(:,:) :: var2
   character(len=*), intent(in) :: filename
   logical, intent(in) :: start_run
   integer, save :: gopenstatus = 1
   integer :: i,k
 
   if (start_run)  then 
      open(unit=iDiag+1,file=filename//'.bin',status="unknown", &
           action="write", iostat=gopenstatus, form="unformatted", &
           position='append' )
   endif   
 
   if(gopenstatus==gok) then
      do i=1,dimx
      	 do k=1,dimz
            write (iDiag+1) var2(i,k)
	 enddo
	 enddo
   else
      write (*,*) filename//'.bin',' NOT open for write. Error'
   endif

   if (timestep == kend) & 
      close(iDiag+1)

END SUBROUTINE cable_diag_data2

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*), intent(in) :: routname
      if(knode_gl==1) & 
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat


END MODULE cable_diag_module




!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

MODULE cable_diag_read_mod !RL: changed to read 3d variable

   IMPLICIT NONE

CONTAINS
  
   SUBROUTINE cable_diagRead( iDiag, basename, dimx, dimy, dimz, timestep, node, &
                        vname, fdata )

      integer :: iDiag, dimx, dimy, dimz, timestep, node
      ! dimx = typically #landpoints over which the var is specified per timestep 
      ! dimy = # timesteps
      ! dimz = # soil levels
      !INTEGER ::dimx, dimy, dimz, iDiag 

      ! passed filename (per field, per processor, at present)
      character(len=*), intent(in) :: basename, vname
      
      character(len=50) :: filename, chnode
 
      ! field (3D at present - one time, one spatial, one vertical ) 
      REAL, DIMENSION(:,:), POINTER :: fdata

 
      write(chnode,10) node
   10 format(i3.3)   
      filename=trim(trim(basename)//trim(chnode))
   
      ! read the binary data and store in 2nd arg 
      CALL read_dat_file(iDiag, TRIM(filename), fdata, dimx, dimz, timestep, dimy)

   END SUBROUTINE cable_diagRead   

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

SUBROUTINE read_dat_file( iDiag, filename, fdata, dimx, dimz, timestep, kend)

   INTEGER :: iDiag
   INTEGER, INTENT(IN) :: dimx, dimz, timestep, kend 
   REAL, DIMENSION(:,:), POINTER :: fdata
   CHARACTER(LEN=*), INTENT(IN) :: filename
   
   INTEGER, PARAMETER :: gok=0
   INTEGER, SAVE :: gopenstatus
   LOGICAL, SAVE :: first_call = .TRUE. 
   INTEGER :: i,k

      IF (first_call) THEN
         OPEN(UNIT=iDiag, FILE=trim(filename)//'.bin', STATUS="unknown", ACTION="read", &
               IOSTAT=gopenstatus, FORM="unformatted" )
         first_call= .FALSE.
      ENDIF   

         IF(gopenstatus==gok) THEN
	    DO i=1,dimx
               DO k=1,dimz
               	  READ(iDiag), fdata(i,k)
	       ENDDO 
	    ENDDO
     
         ELSE
            WRITE (*,*), trim(filename)//'.bin',' NOT found for read'
            STOP
     
         ENDIF

!    CLOSE(iDiag)

END SUBROUTINE read_dat_file 

!==========================================================================!
      

END MODULE cable_diag_read_mod 




