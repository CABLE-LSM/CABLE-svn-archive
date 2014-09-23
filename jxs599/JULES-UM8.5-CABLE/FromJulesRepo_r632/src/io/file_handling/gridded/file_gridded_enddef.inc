#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_gridded_enddef(file)

  USE file_mod, ONLY : file_enddef

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Takes the file out of definition mode - no more dimensions or variables
!   may be defined after this
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file


!-----------------------------------------------------------------------------

! We just need to defer to the underlying file
  CALL file_enddef(file%fh)

  RETURN

END SUBROUTINE file_gridded_enddef
#endif
