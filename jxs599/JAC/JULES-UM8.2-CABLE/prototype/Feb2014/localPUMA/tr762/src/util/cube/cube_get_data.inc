#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE cube_get_data_scalar(cube, value)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a scalar
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: value  ! The scalar value to retrieve


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that the data in the cube is of the correct size to extract a scalar
!-----------------------------------------------------------------------------
  IF ( SIZE(cube%values) /= 1 )                                               &
    CALL log_fatal("cube_get_data_scalar",                                    &
                   "Cube must have size 1 to extract a scalar")

!-----------------------------------------------------------------------------
! Extract the scalar and return it
!-----------------------------------------------------------------------------
  value = cube%values(1)

  RETURN

END SUBROUTINE cube_get_data_scalar


SUBROUTINE cube_get_data_1d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 1d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 1 )                                                &
    CALL log_fatal("cube_get_data_1d",                                        &
                   "values has a different rank to the cube")

! Check that values is the correct size
  IF ( SIZE(values) /= cube%shape(1) )                                        &
    CALL log_fatal("cube_get_data_1d",                                        &
                   "values is not the same size as the cube")

!-----------------------------------------------------------------------------
! The cube's values are already stored as a 1d array
!-----------------------------------------------------------------------------
  values(:) = cube%values(:)

  RETURN

END SUBROUTINE cube_get_data_1d


SUBROUTINE cube_get_data_2d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 2d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 2 )                                                &
    CALL log_fatal("cube_get_data_2d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:2)) )                          &
    CALL log_fatal("cube_get_data_2d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:2))

  RETURN

END SUBROUTINE cube_get_data_2d


SUBROUTINE cube_get_data_3d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 3d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 3 )                                                &
    CALL log_fatal("cube_get_data_3d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:3)) )                          &
    CALL log_fatal("cube_get_data_3d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:3))

  RETURN

END SUBROUTINE cube_get_data_3d


SUBROUTINE cube_get_data_4d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 4d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:,:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 4 )                                                &
    CALL log_fatal("cube_get_data_4d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:4)) )                          &
    CALL log_fatal("cube_get_data_4d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:4))

  RETURN

END SUBROUTINE cube_get_data_4d


SUBROUTINE cube_get_data_5d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 5d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:,:,:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 5 )                                                &
    CALL log_fatal("cube_get_data_5d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:5)) )                          &
    CALL log_fatal("cube_get_data_5d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:5))

  RETURN

END SUBROUTINE cube_get_data_5d


SUBROUTINE cube_get_data_6d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 6d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:,:,:,:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 6 )                                                &
    CALL log_fatal("cube_get_data_6d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:6)) )                          &
    CALL log_fatal("cube_get_data_6d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:6))

  RETURN

END SUBROUTINE cube_get_data_6d


SUBROUTINE cube_get_data_7d(cube, values)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Get the data from a cube where the data is a 7d array
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(data_cube), INTENT(IN) :: cube  ! The cube containing the data

  REAL, INTENT(OUT) :: values(:,:,:,:,:,:,:)  ! The array to put values in


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Check that values is the right shape for the data in the cube
!-----------------------------------------------------------------------------
! First make sure that cube%shape has the correct number of dimensions
  IF ( SIZE(cube%shape) /= 7 )                                                &
    CALL log_fatal("cube_get_data_7d",                                        &
                   "values has a different rank to the cube")

! Check that all the dimensions of values have the correct size
  IF ( .NOT. ALL(SHAPE(values) == cube%shape(1:7)) )                          &
    CALL log_fatal("cube_get_data_7d",                                        &
                   "values is not the same shape as the cube")

!-----------------------------------------------------------------------------
! Reshape the cube's data to the expected shape
!-----------------------------------------------------------------------------
  values = RESHAPE(cube%values, cube%shape(1:7))

  RETURN

END SUBROUTINE cube_get_data_7d
#endif
