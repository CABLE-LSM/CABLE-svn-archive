#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

FUNCTION cube_create(shape) RESULT(cube)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Create a data cube
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: shape(:)  ! The shape of the cube

! Return type
  TYPE(data_cube) :: cube  ! The (empty) cube


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Allocate and populate the cube
!-----------------------------------------------------------------------------
  ALLOCATE(cube%shape(SIZE(shape)))
  ALLOCATE(cube%values(PRODUCT(shape)))

  cube%shape(:)  = shape(:)
  cube%values(:) = 0.0  ! Initialise values to 0

  RETURN

END FUNCTION cube_create
#endif
