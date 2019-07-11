#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ascii/file_ascii_introspect.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ascii_introspect(file)

  USE io_constants, ONLY : MODE_READ

  USE string_utils_mod, ONLY : str_starts_with, str_replace, str_split

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an open file in read mode, try to detect the dimensions and
!   variables in the file and define them on the file
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ascii), INTENT(INOUT) :: file  ! The file to detect items in

! Work variables
  CHARACTER(len=5000) :: line  ! The current line of the file

! Information about the dimensions found in the file
  INTEGER :: ndims
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_FILE)
  INTEGER :: dim_sizes(MAX_DIM_FILE)
  INTEGER :: dim_ids(MAX_DIM_FILE)  ! The ids of the dimensions

! Information about the record dimension
  LOGICAL :: has_record_dim
  CHARACTER(len=MAX_SDF_NAME_LEN) :: record_dim_name

! Information about the variables found in the file
  INTEGER :: nvars
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_names(MAX_VAR_FILE)
  INTEGER :: var_ndims(MAX_VAR_FILE)
  CHARACTER :: var_dim_names(MAX_VAR_FILE, MAX_DIM_VAR)

! Used in calls to str_split
  INTEGER :: nparts
  CHARACTER(len=200) :: parts(2)
  CHARACTER(len=200) :: dim_str

  INTEGER :: dummy  ! Throwaway value

  INTEGER :: var_dim_ids(MAX_DIM_VAR)  ! Used when defining variables
                                       ! The ids of the dimensions found


  INTEGER :: i,j,k  ! Index variables

  INTEGER :: error  ! Error indicator


!-----------------------------------------------------------------------------

! We can only do this in read mode
  IF ( file%mode /= MODE_READ )                                               &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Can only introspect files in read mode")

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  ndims          = 0
  has_record_dim = .FALSE.
  nvars          = 0

!-----------------------------------------------------------------------------
! Try to parse the header
! We can only do this for files that were created using this IO library and
! hence we know the format of the header...
!-----------------------------------------------------------------------------
! Start by rewinding the file to the start
  REWIND(file%unit, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error rewinding file " //                                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Find the line indicating start of dimensions
!-----------------------------------------------------------------------------
  DO
    READ(file%unit, "(A)", IOSTAT=error) line
! Check for end of file
    IF ( error < 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")
! Report a general read error
    IF ( error > 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "Error reading from file - likely bad file format " //   &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! We either keep reading until we get an error, or we exit when we find the
! dimensions line
    IF ( str_starts_with(line, "# Dimensions") ) EXIT
  END DO

!-----------------------------------------------------------------------------
! Read information about the dimensions
!-----------------------------------------------------------------------------
  DO
    READ(file%unit, "(A)", IOSTAT=error) line
! Check for end of file
    IF ( error < 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")
! Report a general read error
    IF ( error > 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "Error reading from file - likely bad file format " //   &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Chop the comment character off
    line = line(2:)

! If we reach an empty line, then we are done with dimensions
    IF ( LEN_TRIM(line) == 0 ) EXIT

! If we have a dimension that we can interpret, we should be able to split
! it into 2 parts about the '='
    CALL str_split(line, "=", nparts, parts)
    IF ( nparts /= 2 )                                                        &
! If we can't get two parts, then we have an unexpected format
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")

! If we get to here, we have successfully found a dimension
! The dimension name is the bit before the '=' (parts(1))
! The dimension length is the bit after the '=' (parts(2))

! Check if it is the record dimension by looking for UNLIMITED as length
    IF ( INDEX(parts(2), "UNLIMITED") > 0 ) THEN
      has_record_dim = .TRUE.
      record_dim_name = ADJUSTL(parts(1))
      CYCLE
    END IF

! If we get to here then we have a proper dimension
    ndims = ndims + 1
    dim_names(ndims) = ADJUSTL(parts(1))
! Try to read an integer out of the size part using list directed IO
    READ(parts(2), *, IOSTAT=error) dim_sizes(ndims)
    IF ( error /= 0 )                                                         &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")
  END DO


!-----------------------------------------------------------------------------
! Find the line indicating start of variables
!-----------------------------------------------------------------------------
  DO
    READ(file%unit, "(A)", IOSTAT=error) line
! Check for end of file
    IF ( error < 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")
! Report a general read error
    IF ( error > 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "Error reading from file - likely bad file format " //   &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! We either keep reading until we get an error, or we exit when we find the
! variables line
    IF ( str_starts_with(line, "# Variables") ) EXIT
  END DO

!-----------------------------------------------------------------------------
! Read information about the variables
!-----------------------------------------------------------------------------
  DO
    READ(file%unit, "(A)", IOSTAT=error) line
! Check for end of file
    IF ( error < 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")
! Report a general read error
    IF ( error > 0 )                                                          &
      CALL log_fatal("file_ascii_introspect",                                 &
                     "Error reading from file - likely bad file format " //   &
                     "(IOSTAT=" // TRIM(to_string(error)) // ")")

! If we've gone past the end of the header, we are done
    IF ( line(1:1) /= "#" ) EXIT

! Chop the comment character off
    line = line(2:)

! Ignore empty lines and lines with an '=' character in (attribute lines)
    IF ( LEN_TRIM(line) == 0 .OR. INDEX(line, "=") > 0 ) CYCLE

! If we have a variable that we can interpret, we should be able to split
! it into 2 parts about the '('
    CALL str_split(line, "(", nparts, parts)
    IF ( nparts /= 2 )                                                        &
! If we can't get two parts, then we have an unexpected format
      CALL log_fatal("file_ascii_introspect",                                 &
                     "File does not have expected format")

! If we get to here, we have successfully found a variable
    nvars = nvars + 1
!   parts(1) contains the variable name - all we need to do is left align it
    var_names(nvars) = ADJUSTL(parts(1))

! We want to get rid of the trailing ') ;' from parts(2) to get just the
! dimension names separated by commas
    dim_str = str_replace(parts(2), ") ;", " ")

! If there is no string left, then we have a variable with no non-record
! dimensions
    IF ( LEN_TRIM(dim_str) == 0 ) THEN
      var_ndims(nvars) = 0
      CYCLE
    END IF

! Now we can extract the dimensions by just splitting about commas
    CALL str_split(dim_str, ",", var_ndims(nvars), var_dim_names(nvars,:))
  END DO

!-----------------------------------------------------------------------------
! We have finished reading the file, so rewind it for consistency
!-----------------------------------------------------------------------------
  REWIND(file%unit, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_introspect",                                   &
                   "Error rewinding file " //                                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Now we know what dimensions and variables are in the file, we can define
! them
!-----------------------------------------------------------------------------
! First, dimensions
  DO i = 1,ndims
    dim_ids(i) = file_ascii_def_dim(file, dim_names(i), dim_sizes(i))
  END DO
! And record dimension
  IF ( has_record_dim )                                                       &
    dummy = file_ascii_def_record_dim(file, record_dim_name)

! Then define variables
  DO i = 1,nvars
! First, gather up the dimension ids
    DO j = 1,var_ndims(i)
! If the variable has a dimension that isn't in the file, error out
      IF ( .NOT. ANY(dim_names == var_dim_names(i,j)) )                       &
        CALL log_fatal("file_ascii_introspect",                               &
                       "Variable has a dimension not defined in file")

! Otherwise, get the corresponding dimension id
      DO k = 1,ndims
        IF ( var_dim_names(i,j) == dim_names(k) ) THEN
          var_dim_ids(j) = dim_ids(k)
          EXIT
        END IF
      END DO
    END DO

! Now we have the dimension ids for the variable, we can define it
! We throw away the id as we don't care what it is
    dummy = file_ascii_def_var(                                               &
      file, var_names(i), var_dim_ids(1:var_ndims(i)), has_record_dim         &
    )
  END DO

! Finally, we take the file out of define mode
  CALL file_ascii_enddef(file)

  RETURN

END SUBROUTINE file_ascii_introspect
#endif
