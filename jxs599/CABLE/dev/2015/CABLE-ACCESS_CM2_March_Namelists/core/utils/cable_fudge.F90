
MODULE cable_fudge_module

interface fudge_out
   module procedure fudge_out_r2D, fudge_out_r1D, fudge_out_r3D, fudge_out_i2D
End interface fudge_out

CONTAINS

SUBROUTINE fudge_out_i2D( i,j, var, varname, vzero, vval )
   ! interfaces on these
   integer :: i,j
   integer, dimension(:,:) :: var
   ! ft changes with interface
   character(len=*), parameter :: &
      ft = '(  "fudge: ", A10, "(", I2.1, ",", I2.1, X, ") = ", I1.1 )'
   
   character(len=*) :: varname
   logical :: vzero
   integer :: vval
   
   ! content changes with interface
   var = var(i,j) 
   if( (vzero) ) var = vval
   write (6, ft) varname,i, var(i,j)
End SUBROUTINE fudge_out_i2D 


SUBROUTINE fudge_out_r1D( i, var, varname, vzero, vval )
   ! interfaces on these
   integer :: i
   real, dimension(:) :: var
   ! ft changes with interface
   character(len=*), parameter :: &
      ft = '(  "fudge: ", A10, "(", I2.1, X, ") = ", F15.3 )'
   
   character(len=*) :: varname
   logical :: vzero
   real :: vval

   ! content changes with interface
   var = var(i) 
   if( (vzero) ) var = vval
   write (6, ft) varname,i, var(i)
End SUBROUTINE fudge_out_r1D 

SUBROUTINE fudge_out_r2D( i,j, var, varname, vzero, vval )
   ! interfaces on these
   integer :: i,j
   real, dimension(:,:) :: var
   ! ft changes with interface
   character(len=*), parameter :: &
      ft = '(  "fudge: ", A10, "(", I2.1, ",", I2.1, X, ") = ", F15.3 )'
   
   character(len=*) :: varname
   logical :: vzero
   real :: vval
   
   ! content changes with interface
   var = var(i,j) 
   if( (vzero) ) var = vval
   write (6, ft) varname,i,j, var(i,j)
End SUBROUTINE fudge_out_r2D 

SUBROUTINE fudge_out_r3D( i,j,k, var, varname, vzero, vval )
   ! interfaces on these
   integer :: i,j,k
   real, dimension(:,:,:) :: var
   ! ft changes with interface
   character(len=*), parameter :: &
      ft = '(  "fudge: ", A10, "(",  I2.1, ",",I2.1, ",", I2.1, X, ") = ", F15.3 )'
   
   character(len=*) :: varname
   logical :: vzero
   real :: vval
   
   ! content changes with interface
   var = var(i,j,k) 
   if( (vzero) ) var = vval
   write (6, ft) varname,i,j,k, var(i,j,k)
End SUBROUTINE fudge_out_r3D 



END MODULE cable_fudge_module
