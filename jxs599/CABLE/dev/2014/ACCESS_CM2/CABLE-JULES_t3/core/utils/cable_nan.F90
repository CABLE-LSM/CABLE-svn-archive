MODULE cable_nan_module
   implicit none
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1

   interface cable_NaN
      module procedure cable_NaN1, cable_NaN2
   end interface cable_NaN

CONTAINS

!--- cable status NaN
!==========================================================================!

SUBROUTINE cable_NaN1(fname,field,mype)

real, dimension(:,:) :: field
character(len=*), dimension(:) :: fname

integer, optional :: mype 
integer :: i,j
logical :: NoNaN, check
integer :: n,m

   n = size(fname)
   m = size(field,2)
   check = .FALSE.   

   do i=1, n  
      
      NoNaN = .TRUE.
      
      do j=1, m  

         call isnan(  field(i,j), check )

         if( check ) then 
            print *, ""
            print *, "CABLE_log: "
            if( present(mype) ) print *, "proc # ", mype
            print *, "   Element: ",j
            print *, "   of field ", fname(i)
            print *, "   is NaN"
            print *, "End CABLE_log: "
            print *, ""
            NoNaN = .FALSE.
         end if 

      enddo
      
      if(NoNaN) then 
         print *, ""
         print *, "CABLE_log: "
         if( present(mype) ) print *, "proc # ", mype
         print *, '   Field: ', fname(i)
         print *, "   is clear of NaNs"
         print *, "End CABLE_log: "
         print *, ""
      end if
       
   enddo

END SUBROUTINE cable_NaN1


SUBROUTINE cable_NaN2(fname,field,mype)

real, dimension(:,:,:) :: field
character(len=*), dimension(:) :: fname

integer, optional :: mype 
integer :: i,j,k
logical :: NoNaN, check
integer :: n,m,op

   n = size(fname)
   m = size(field,2)
   op = size(field,3)
   check = .FALSE.   

   do i=1, n  
      
      NoNaN = .TRUE.
      
      do j=1, m  
      
         do k=1, op  

            call isnan(  field(i,j,k), check )
   
            if( check ) then 
               print *, ""
               print *, "CABLE_log: "
               if( present(mype) ) print *, "proc # ", mype
               print *, "   Element: ",j,op
               print *, "   of field ", fname(i)
               print *, "   is NaN"
               print *, "End CABLE_log: "
               print *, ""
               NoNaN = .FALSE.
            end if 

         enddo
      
      enddo
      
      if(NoNaN) then 
         print *, ""
         print *, "CABLE_log: "
         if( present(mype) ) print *, "proc # ", mype
         print *, '   Field: ', fname(i)
         print *, "   is clear of NaNs"
         print *, "End CABLE_log: "
         print *, ""
      end if
       
   enddo

END SUBROUTINE cable_NaN2





subroutine isnan(var, check) 
   real :: var 
   logical :: check

   if (var .ne. var) then 
      check = .true. 
   else 
      check = .false. 
   end if 

end subroutine isnan


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


END MODULE cable_nan_module
