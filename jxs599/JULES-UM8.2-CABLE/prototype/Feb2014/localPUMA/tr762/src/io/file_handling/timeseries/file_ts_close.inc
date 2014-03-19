#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ts_close(file)

  USE file_gridded_mod, ONLY : file_gridded_close

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Closes and frees any resources consumed by the given file object
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ts), INTENT(INOUT) :: file ! The file to close


! Work variables
  INTEGER :: i  ! Index variable


!-----------------------------------------------------------------------------

! Make sure the currently open file is closed
  CALL file_gridded_close(file%open_file)

! Free the memory associated with the files and file_times arrays if it was
! allocated (i.e. if not using time templating)
  IF ( .NOT. file%use_time_template ) THEN
    DEALLOCATE(file%files)
    NULLIFY(file%files)

    DEALLOCATE(file%file_times)
    NULLIFY(file%file_times)
  END IF

! Free the memory associated with the levels dimensions for variables
  DO i = 1,file%nvars
    IF ( ASSOCIATED(file%vars(i)%levels_dims) ) THEN
      DEALLOCATE(file%vars(i)%levels_dims)
      NULLIFY(file%vars(i)%levels_dims)
    END IF
  END DO

  RETURN

END SUBROUTINE file_ts_close
#endif
