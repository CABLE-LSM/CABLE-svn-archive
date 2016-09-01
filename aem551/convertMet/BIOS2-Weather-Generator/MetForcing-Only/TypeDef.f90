MODULE TypeDef
!------------------------------------------------------------------------------- 
! PRB: 18-06-2000
! This module explicitly defines the sizes of variable types
!-------------------------------------------------------------------------------
  implicit none
  save
! Define integer kind parameters to accommodate the range of numbers
! usually 
! associated with 4, 2, and 1 byte integers. 
  integer,parameter :: i4b = selected_int_kind(9)
  integer,parameter :: i2b = selected_int_kind(4)
  integer,parameter :: i1b = selected_int_kind(2)
! Define single and double precision real kind parameters: 
! * Kind(1.0)   defines sp as the machine's default size for single
! precision
! * Kind(1.0d0) defines dp as the machine's default size for double
! precision
  integer,parameter :: sp  = kind(1.0)
  integer,parameter :: dp  = kind(1.0d0)
! lgt is set to the default kind required for representing logical
! values. 
  integer,parameter :: lgt = kind(.true.)

END MODULE TypeDef
