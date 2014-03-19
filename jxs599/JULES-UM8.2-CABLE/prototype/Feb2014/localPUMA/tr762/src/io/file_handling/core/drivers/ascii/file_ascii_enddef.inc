#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE file_ascii_enddef(file)

  USE io_constants, ONLY : MODE_READ, MODE_WRITE, ATTR_GLOBAL

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
  TYPE(file_ascii), INTENT(INOUT) :: file

! Work variables
  CHARACTER(len=100) :: var_dims_str  ! Used when writing headers to store
                                      ! variables dimension names as a
                                      ! string for writing

  INTEGER :: buffer_size  ! The size of the buffer
  INTEGER :: offset  ! Used in calculation of offsets
  INTEGER :: error  ! Error indicator

  INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------

  file%define_mode = .FALSE.

!-----------------------------------------------------------------------------
! Allocate the buffer - it must be big enough to hold one records worth of
! data for each variable
!-----------------------------------------------------------------------------
! To get the size of the buffer, just sum the sizes of the defined variables
  buffer_size = SUM(file%var_sizes(1:file%nvars))
  ALLOCATE(file%buffer(buffer_size), STAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("file_ascii_enddef",                                       &
                   "Error allocating memory - " //                            &
                   "(STAT=" // TRIM(to_string(error)) // ")")

! Now we know the size of the buffer, we can make a better guess at the
! record length we will need. For now, we'll assume no more than 30 characters
! per buffer element - this should be very conservative
  file%record_len = buffer_size * 30

!-----------------------------------------------------------------------------
! Calculate offsets for each variable
! It is assumed that they appear in the file in the order in which they
! were defined
!-----------------------------------------------------------------------------
  offset = 1
  DO i = 1,file%nvars
    file%var_offsets(i) = offset
    offset = offset + file%var_sizes(i)
  END DO


! Now we do different things depending on mode
  SELECT CASE ( file%mode )
    CASE ( MODE_READ )
! In read mode, we fill the buffer with the first data so we are ready to read
      CALL file_ascii_fill_buffer(file)

    CASE ( MODE_WRITE )
!-----------------------------------------------------------------------------
! In write mode, we write the headers describing the variables in the file
!
! The header is of the form:
!
! #
! # CREATED BY JULES LSM
! #
! # Global attributes:
! #     global_attr_1 = foobar ;
! #     global_attr_2 = baz ;
! #
! # Dimensions:
! #     dim1 = size1 ;
! #     dim2 = size2 ;
! #     dim3 = size3 ;
! #     dim4 = UNLIMITED ;
! #
! # Variables:
! #     variable1(dim1,dim2) ;
! #         variable1_attr_1 = foobaz ;
! #         variable1_attr_2 = blah ;
! #
! #     variable2(dim1) ;
! #
! #     variable3(dim1,dim2,dim3) ;
! #
!-----------------------------------------------------------------------------
! First write the created by header
      WRITE(file%unit, "(A)") "#"
      WRITE(file%unit, "(A)") "# CREATED BY JULES LSM"
      WRITE(file%unit, "(A)") "#"

!-----------------------------------------------------------------------------
! Write global attributes
!-----------------------------------------------------------------------------
      WRITE(file%unit, "(A)") "# Global attributes:"
      DO i = 1,file%nattrs
        IF ( file%attr_var_ids(i) == ATTR_GLOBAL )                            &
          WRITE(file%unit, "(A)") "#     " // TRIM(file%attr_values(i)) // " ;"
      END DO

      WRITE(file%unit, "(A)") "#"

!-----------------------------------------------------------------------------
! Write information about the dimensions in the file
!-----------------------------------------------------------------------------
      WRITE(file%unit, "(A)") "# Dimensions:"
      DO i = 1,file%ndims
        WRITE(file%unit, "(A)") "#     " // TRIM(file%dim_names(i)) //        &
                                " = " // TRIM(to_string(file%dim_sizes(i))) // " ;"

      END DO
! Write information about the record dimension
      IF ( file%has_record_dim )                                              &
        WRITE(file%unit, "(A)") "#     " // TRIM(file%record_dim_name) //     &
                                " = UNLIMITED ;"

      WRITE(file%unit, "(A)") "#"

!-----------------------------------------------------------------------------
! Write information about the variables in the file
!-----------------------------------------------------------------------------
      WRITE(file%unit, "(A)") "# Variables:"
      DO i = 1,file%nvars
! Build the dimension names string for the variable
        var_dims_str = ""
        DO j = 1,file%var_ndims(i)
          var_dims_str = TRIM(var_dims_str) // file%dim_names( file%var_dims(i,j) )

! Add a comma unless this is the last dimension that the variable has, including
! the addition of the record dim below
          IF ( j < file%var_ndims(i) .OR. file%has_record_dim )               &
            var_dims_str = TRIM(var_dims_str) // ","
        END DO
! Add the record dimension if it is in use
        IF ( file%has_record_dim ) THEN
          var_dims_str = TRIM(var_dims_str) // file%record_dim_name
        END IF

        WRITE(file%unit, "(A)") "#     " // TRIM(file%var_names(i)) //        &
                                "(" // TRIM(var_dims_str) // ") ;"

! Write info about each of the variables attributes
        DO j = 1,file%nattrs
          IF ( file%attr_var_ids(j) == i )                                    &
            WRITE(file%unit, "(A)") "#         " // TRIM(file%attr_values(j)) // " ;"
        END DO

        WRITE(file%unit, "(A)") "#"
      END DO

! Set the format string to use for output, now we know how many elements are in
! the buffer
      out_format_str = "(" // TRIM(to_string(buffer_size)) // "G15.7e2)"

  END SELECT

  RETURN

END SUBROUTINE file_ascii_enddef
#endif
