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
   IMPLICIT NONe
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
  
   !--- subrs overloaded to respond to call cable_diag 
   INTERFACE cable_diag
      MODULE PROCEDURE cable_diag_zero
   END INTERFACE cable_diag
  
CONTAINS

!==========================================================================!
! cable_diag1/2/3 call subrs to write filename.dat which contains description
! of data and format etc., and filename.bin containing the data   
!==========================================================================!

SUBROUTINE cable_diag_zero( Nvars, basename, dimx, dimy, timestep, node,       &
                        vname1, var1, level )
   INTEGER, INTENT(IN) :: Nvars,dimx, dimy, timestep,node
   REAL, INTENT(IN), DIMENSION(:) :: var1
   INTEGER :: i=0, fstatus
   INTEGER, DIMENSION(13) :: fvalues 
   CHARACTER(LEN=*), INTENT(IN) :: basename, vname1, level
   CHARACTER(LEN=30) :: filename, chnode
 
      WRITE(chnode,10) node
   10 FORMAT(i2.2)   
      filename=TRIM(TRIM(basename)//TRIM(chnode))
   
      fstatus = STAT( filename, fvalues ) 
      write(6,*) 'file status', fstatus 
      stop 
      ! if first timestep in run and file already exists stop 
      IF( timestep == 1 ) THEN
         WRITE(6,*) 'CABLE_log: file already exists', filename
         STOP
      ENDIF   
      ! if NOT first timestep in run and file size is too large write WARNING 
      IF (timestep == 1)                                                       & 
      IF( timestep == 1 .AND. knode_gl==1 )                                                       & 
      
      ! if NOT first timestep in run and file size is too large write WARNING 
      IF (timestep == 1)                                                       & 
         CALL cable_diag_desc1( Nvars, TRIM(filename), dimx, dimy, vname1 )
      
      CALL cable_diag_data1( Nvars, TRIM(filename), dimx, timestep, dimy,      &
                             var1 )
END SUBROUTINE cable_diag1

!=============================================================================!
!=============================================================================!

SUBROUTINE cable_diag_desc1( Nvars, filename, dimx, dimy, vname1 )

   INTEGER, INTENT(IN) :: Nvars,dimx,dimy 
   CHARACTER(LEN=*), INTENT(IN) :: filename, vname1
   INTEGER, SAVE :: gopenstatus = 1

     OPEN(UNIT=713941,FILE=filename//'.dat', STATUS="replace",                 &
          ACTION="write", IOSTAT=gopenstatus )
     
      IF(gopenstatus==gok) THEN
            WRITE (713941,*) 'Number of var(s): '
            WRITE (713941,*) Nvars
            WRITE (713941,*) 'Name of var(s): '
            WRITE (713941,7139) vname1 
 7139       FORMAT(a)            
            WRITE (713941,*) 'dimension of var(s) in x: '
            WRITE (713941,*) dimx 
            WRITE (713941,*) 'dimension of var(s) in y: '
            WRITE (713941,*) dimy 
      ELSE
         WRITE (*,*), filename//'.dat',' Error: unable to write'
      ENDIF
      
   CLOSE(713941)
  
END SUBROUTINE cable_diag_desc1


SUBROUTINE cable_diag_data1( Nvars, filename, dimx, timestep, kend, var1  )

   INTEGER, INTENT(IN) :: Nvars, dimx, timestep, kend
   REAL, INTENT(IN), DIMENSION(:) :: var1
   CHARACTER(LEN=*), INTENT(IN) :: filename
   INTEGER, SAVE :: gopenstatus = 1

   IF (timestep == 1)  THEN 
      OPEN(UNIT=713942,FILE=filename//'.bin',STATUS="unknown",                 &
           ACTION="write", IOSTAT=gopenstatus, FORM="unformatted",             &
           POSITION='append' )
   ENDIF   
 
   IF(gopenstatus==gok) THEN
         WRITE (713942) var1
   ELSE
      WRITE (*,*) filename//'.bin',' NOT open for write. Error'
   ENDIF

   IF (timestep == kend)                                                       & 
      CLOSE(713942)

END SUBROUTINE cable_diag_data1

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   USE cable_common_module, ONLY : ktau_gl, knode_gl

   CHARACTER(LEN=*), INTENT(IN) :: routname
      IF(knode_gl==1)                                                          & 
         WRITE(6,*) 'CABLE_log @  ', routname, ktau_gl

END SUBROUTINE cable_stat


END MODULE cable_diag_module



