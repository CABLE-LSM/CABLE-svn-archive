#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_close(file)

  USE file_mod, ONLY : file_close

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given gridded file object
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file ! The file to close

! Work variables
  INTEGER :: i  ! Loop index


!-----------------------------------------------------------------------------


! Close the underlying file handle
  CALL file_close(file%fh)

! Deallocate the sizes of any variables
  DO i = 1,file%nvars
    IF ( ASSOCIATED(file%vars(i)%lev_sizes) ) THEN
      DEALLOCATE(file%vars(i)%lev_sizes)
      NULLIFY(file%vars(i)%lev_sizes)
    END IF
  END DO

  RETURN

END SUBROUTINE file_gridded_close
#endif
