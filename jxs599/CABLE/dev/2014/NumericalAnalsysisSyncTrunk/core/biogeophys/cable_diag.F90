!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.cawcr.gov.au/projects/access/cable
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
   integer :: iDiag 
   integer, SAVE :: pDiag=713 
   integer :: dimx, dimy, timestep,node
   real, dimension(:) :: var1
   integer :: Nvars=1 !this WAS input
   integer :: i=0
   character(len=*) :: basename, vname1
   character(len=30) :: filename, chnode
 
      IF(iDiag==0) tHEN
         pDiag = pDiag+2  
         iDiag=pDiag
      ENDIF
         
      write(chnode,10) node
   10 format(i2.2)   
      filename=trim(trim(basename)//trim(chnode))
      
      if (timestep == 1) & 
         call cable_diag_desc1( iDiag, trim(filename), dimx, dimy, vname1 )
      
      call cable_diag_data1( iDiag, trim(filename), dimx, timestep, dimy, &
                             var1 )
END SUBROUTINE cable_diag1

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc1( iDiag, filename, dimx, dimy, vname1 )

   integer :: iDiag,dimx,dimy 
   integer, PARAMETER :: Nvars=1
   character(len=*) :: filename, vname1
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

   integer :: iDiag, dimx, timestep, kend
   integer, PARAMETER :: Nvars=1
   real, dimension(:) :: var1
   character(len=*) :: filename
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

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*) :: routname
      if(knode_gl==1) & 
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat


END MODULE cable_diag_module



