MODULE TypeDef
  !-------------------------------------------------------------------------------
  ! This module explicitly defines the sizes of variable types
  !-------------------------------------------------------------------------------
  IMPLICIT NONE
  SAVE
  ! Define integer kind parameters to accommodate the range of numbers usually
  ! associated with 4, 2, and 1 byte integers.
  INTEGER,PARAMETER :: i4b = SELECTED_INT_KIND(9)
  INTEGER,PARAMETER :: i2b = SELECTED_INT_KIND(4)
  INTEGER,PARAMETER :: i1b = SELECTED_INT_KIND(2)
  ! Define single and double precision real kind parameters:
  ! * Kind(1.0)   defines sp as the machine's default size for single precision
  ! * Kind(1.0d0) defines dp as the machine's default size for double precision
  INTEGER,PARAMETER :: sp  = KIND(1.0)
  INTEGER,PARAMETER :: dp  = KIND(1.0d0)
  ! lgt is set to the default kind required for representing logical values.
  INTEGER,PARAMETER :: lgt = KIND(.TRUE.)

END MODULE TypeDef


