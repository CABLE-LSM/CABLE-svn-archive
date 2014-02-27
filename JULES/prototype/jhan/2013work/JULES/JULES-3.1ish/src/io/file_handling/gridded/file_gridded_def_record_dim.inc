#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/gridded/file_gridded_def_record_dim.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

FUNCTION file_gridded_def_record_dim(file, dim_name) RESULT(dim_id)

  USE file_mod, ONLY : file_def_record_dim

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines a record dimension on the given file, returning the dimension id
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file
                                  ! The file to define the dimension on
  CHARACTER(len=*), INTENT(IN) :: dim_name
                                  ! The name of the record dimension

! Return type
  INTEGER :: dim_id               ! The dimension id


!-----------------------------------------------------------------------------

! Just defer to the underlying file
  dim_id = file_def_record_dim(file%fh, dim_name)

  RETURN

END FUNCTION file_gridded_def_record_dim
#endif
