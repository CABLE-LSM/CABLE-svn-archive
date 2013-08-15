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

!==========================================================================!
!--- cable generic print status
!==========================================================================!

SUBROUTINE cable_stat( routname)
   use cable_common_module, only : ktau_gl, knode_gl

   character(len=*), intent(in) :: routname
      if(knode_gl==1) & 
         write(6,*) 'CABLE@  ', routname, ktau_gl

END SUBROUTINE cable_stat

!==========================================================================!
!--- cable status NaN
!==========================================================================!


SUBROUTINE cable_NaN(fname,field,logn)

real, dimension(:,:) :: field
character(len=*), dimension(:) :: fname
integer, intent(in) :: logn
!real, dimension(:) :: mfield

integer :: i,j
logical :: isnan 

logical,save :: NoNan
integer :: n,m
character(len=*), PARAMETER :: & 
      logfmt1 = "( 'Element(',I6.1,') of field ',A,' is NaN' )"  , &  
      logfmt21 = "( 'Field ',A,' is clear of NaNs' )"    


   n = size(fname)
   m = size(field,1 )
   
   do i=1, n  
      NoNan = .TRUE.
      do j=1, m  
         if (isnan(  field(i,j) )) then 
            write(logn,*) "CABLE_log: "
            write(logn,logfmt1) j,fname(i)  
            NoNan = .FALSE.
         end if 
      enddo
   
      if(NoNaN) then 
         write(logn,*) "CABLE_log: "
         write(logn,logfmt21) fname(i) 
      end if 
   
   enddo

END SUBROUTINE cable_NaN


LOGICAL FUNCTION isnan(var) 
   real var 

   if (var .ne. var) then 
      isnan = .true. 
   else 
      isnan = .false. 
   end if 

   return 

END FUNCTION isnan 

!logical function isinf(a) 
!real a 
!
!!if ((a*0).ne.0) then 
!!isinf = .true. 
!!else 
!!isinf = .false. 
!!end if 
!!return 
!!end 


SUBROUTINE cable_farray( mp, CheckNames, CheckFields, &
                  n1,f1, n2,f2, n3,f3, n4,f4, n5,f5, n6,f6, n7,f7, & 
                  n8,f8, n9,f9, n10,f10, n11,f11, n12,f12, n13,f13, &
                  n14,f14, n15,f15, n16,f16, n17,f17, n18,f18, & 
                  n19,f19, n20,f20, n21,f21, n22,f22, n23,f23, & 
                  n24,f24, n25,f25, n26,f26, n27,f27, n28,f28, & 
                  n29,f29, n30,f30, n31,f31, n32,f32, n33,f33, & 
                  n34,f34, n35,f35, n36,f36, n37,f37, n38,f38, & 
                  n39,f39, n40,f40, n41,f41, n42,f42, n43,f43, & 
                  n44,f44, n45,f45, n46,f46, n47,f47, n48,f48, & 
                  n49,f49, n50,f50 & 
               )

   USE cable_common_module, ONLY : fnames => cable_farray_names, &
                                   ffields => cable_farray_fields, &
                                   fn_max => cable_farray_nmax

   integer :: mp
    
   character(len=*), optional :: & 
                        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, & 
                        n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
                        n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
                        n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
                        n41, n42, n43, n44, n45, n46, n47, n48, n49, n50  

   real, dimension(:), optional :: & 
                        f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, & 
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,& 
                        f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,& 
                        f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,& 
                        f41, f42, f43, f44, f45, f46, f47, f48, f49, f50
   
   character(len=30), dimension(:), allocatable :: CheckNames 
   real, dimension(:,:), allocatable :: CheckFields

   integer :: i, k

   allocate( ffields( fn_max, mp) ) 
   i = 0

   if( present (n1) .AND. present(f1) ) then
      CALL fill_farray( n1, f1, i )   
   else
      print *, "CABLE_log: cable_farray missing dummy args"
      return
   endif
   
   if( present (n2) .AND. present(f2) ) CALL fill_farray( n2, f2, i )   
   if( present (n3) .AND. present(f3) ) CALL fill_farray( n3, f3, i )   
   if( present (n4) .AND. present(f4) ) CALL fill_farray( n4, f4, i )   
   if( present (n5) .AND. present(f5) ) CALL fill_farray( n5, f5, i )   
   if( present (n6) .AND. present(f6) ) CALL fill_farray( n6, f6, i )   
   if( present (n7) .AND. present(f7) ) CALL fill_farray( n7, f7, i )   
   if( present (n8) .AND. present(f8) ) CALL fill_farray( n8, f8, i )   
   if( present (n9) .AND. present(f9) ) CALL fill_farray( n9, f9, i )   
   if( present (n10) .AND. present(f10) ) CALL fill_farray( n10, f10, i )   
   if( present (n11) .AND. present(f11) ) CALL fill_farray( n11, f11, i )   
   if( present (n12) .AND. present(f12) ) CALL fill_farray( n12, f12, i )   
   if( present (n13) .AND. present(f13) ) CALL fill_farray( n13, f13, i )   
   if( present (n14) .AND. present(f14) ) CALL fill_farray( n14, f14, i )   
   if( present (n15) .AND. present(f15) ) CALL fill_farray( n15, f15, i )   
   if( present (n16) .AND. present(f16) ) CALL fill_farray( n16, f16, i )   
   if( present (n17) .AND. present(f17) ) CALL fill_farray( n17, f17, i )   
   if( present (n18) .AND. present(f18) ) CALL fill_farray( n18, f18, i )   
   if( present (n19) .AND. present(f19) ) CALL fill_farray( n19, f19, i )   
   if( present (n20) .AND. present(f20) ) CALL fill_farray( n20, f20, i )   
   if( present (n21) .AND. present(f21) ) CALL fill_farray( n21, f21, i )   
   if( present (n22) .AND. present(f22) ) CALL fill_farray( n22, f22, i )   
   if( present (n23) .AND. present(f23) ) CALL fill_farray( n23, f23, i )   
   if( present (n24) .AND. present(f24) ) CALL fill_farray( n24, f24, i )   
   if( present (n25) .AND. present(f25) ) CALL fill_farray( n25, f25, i )   
   if( present (n26) .AND. present(f26) ) CALL fill_farray( n26, f26, i )   
   if( present (n27) .AND. present(f27) ) CALL fill_farray( n27, f27, i )   
   if( present (n29) .AND. present(f29) ) CALL fill_farray( n29, f29, i )   
   if( present (n30) .AND. present(f30) ) CALL fill_farray( n30, f30, i )   
   if( present (n31) .AND. present(f31) ) CALL fill_farray( n31, f31, i )   
   if( present (n32) .AND. present(f32) ) CALL fill_farray( n32, f32, i )   
   if( present (n33) .AND. present(f33) ) CALL fill_farray( n33, f33, i )   
   if( present (n34) .AND. present(f34) ) CALL fill_farray( n34, f34, i )   
   if( present (n35) .AND. present(f35) ) CALL fill_farray( n35, f35, i )   
   if( present (n36) .AND. present(f36) ) CALL fill_farray( n36, f36, i )   
   if( present (n37) .AND. present(f37) ) CALL fill_farray( n37, f37, i )   
   if( present (n38) .AND. present(f38) ) CALL fill_farray( n38, f38, i )   
   if( present (n39) .AND. present(f39) ) CALL fill_farray( n39, f39, i )   
   if( present (n40) .AND. present(f40) ) CALL fill_farray( n40, f40, i )   
   if( present (n41) .AND. present(f41) ) CALL fill_farray( n41, f41, i )   
   if( present (n42) .AND. present(f42) ) CALL fill_farray( n42, f42, i )   
   if( present (n50) .AND. present(f50) ) CALL fill_farray( n50, f50, i )   
   
   allocate( CheckNames(i) )
   allocate( CheckFields(i,mp) )

   do k=1,i
      CheckNames(k) = fnames(k) 
      CheckFields(k,:) = ffields(k,:)
   enddo        

   deallocate( ffields )
    
END SUBROUTINE cable_farray


subroutine fill_farray( n, f, i )   
   
   USE cable_common_module, ONLY : fnames => cable_farray_names, &
                                   ffields => cable_farray_fields, &
                                   fn_max => cable_farray_nmax
   
   character(len=*) :: n
   real, dimension(:) :: f
   integer :: i
      
      i=i+1
      fnames(i) = n
      ffields(i,:) = f
   
end subroutine fill_farray 


END MODULE cable_diag_module



