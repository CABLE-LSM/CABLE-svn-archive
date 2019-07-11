#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

SUBROUTINE get_var_levs_dims(var_id, ndims, dim_names_in, dim_names_out, dim_sizes)

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the levels dimensions
!   that it uses
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
  INTEGER, INTENT(OUT), OPTIONAL :: ndims
                                 ! The number of levels dimensions that the
                                 ! variable has
  CHARACTER(len=*), INTENT(OUT), OPTIONAL :: dim_names_in(:)
                                 ! The dimension names to use for the variable
                                 ! in input files
  CHARACTER(len=*), INTENT(OUT), OPTIONAL :: dim_names_out(:)
                                 ! The dimension names to use for the variable
                                 ! in output files
  INTEGER, INTENT(OUT), OPTIONAL :: dim_sizes(:)
                                 ! The dimension sizes to for the variable


!-----------------------------------------------------------------------------


! var_id is an index into the metadata array

! We decide what to return based on the 'variable type' of the requested variable
  SELECT CASE ( metadata(var_id)%var_type )
    CASE ( VAR_TYPE_SURFACE )
! Just indicate that there are no vertical dimensions
      IF ( PRESENT(ndims) ) ndims = 0

! Single extra dimension cases for different sizes
    CASE ( VAR_TYPE_PFT )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = pft_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = PFT_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = pft_dim_size

    CASE ( VAR_TYPE_NVG )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = nvg_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = NVG_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = nvg_dim_size

    CASE ( VAR_TYPE_TYPE )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = type_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = TYPE_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = type_dim_size

    CASE ( VAR_TYPE_TILE )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = tile_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = TILE_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = tile_dim_size

    CASE ( VAR_TYPE_SOIL )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = soil_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = SOIL_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = soil_dim_size

    CASE ( VAR_TYPE_SCPOOL )
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) ) dim_names_in(1) = scpool_dim_name
      IF ( PRESENT(dim_names_out) ) dim_names_out(1) = SCPOOL_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) ) dim_sizes(1) = scpool_dim_size

    CASE ( VAR_TYPE_SNOW )
! Snow variables have 2 levels dimensions
      IF ( PRESENT(ndims) ) ndims = 2
      IF ( PRESENT(dim_names_in) )                                            &
        dim_names_in(1:2) = (/ tile_dim_name, snow_dim_name /)
      IF ( PRESENT(dim_names_out) )                                           &
        dim_names_out(1:2) = (/ TILE_DIM_NAME_OUT, SNOW_DIM_NAME_OUT /)
      IF ( PRESENT(dim_sizes) )                                               &
        dim_sizes(1:2) = (/ tile_dim_size, snow_dim_size /)

    ! These CABLE var types can piggy back on existing types if so ordained
    CASE ( VAR_TYPE_CABLE_SOIL )
! CABLE soil variables have 2 levels dimensions
      IF ( PRESENT(ndims) ) ndims = 2
      IF ( PRESENT(dim_names_in) )                                            &
        dim_names_in(1:2) = (/ tile_dim_name, soil_dim_name /)
      IF ( PRESENT(dim_names_out) )                                           &
        dim_names_out(1:2) = (/ TILE_DIM_NAME_OUT, SOIL_DIM_NAME_OUT /)
      IF ( PRESENT(dim_sizes) )                                               &
        dim_sizes(1:2) = (/ tile_dim_size, soil_dim_size /)

    CASE ( VAR_TYPE_CABLE_SNOW )
! CABLE snow variables have 2 levels dimensions
      IF ( PRESENT(ndims) ) ndims = 2
      IF ( PRESENT(dim_names_in) )                                            &
        dim_names_in(1:2) = (/ tile_dim_name, snow_dim_name /)
      IF ( PRESENT(dim_names_out) )                                           &
        dim_names_out(1:2) = (/ TILE_DIM_NAME_OUT, SNOW_DIM_NAME_OUT /)
      IF ( PRESENT(dim_sizes) )                                               &
        dim_sizes(1:2) = (/ tile_dim_size, snow_dim_size /)

    CASE ( VAR_TYPE_CABLE1l_SNOW )
! CABLE snow variables per tile only have 1 level dimensions
      IF ( PRESENT(ndims) ) ndims = 1
      IF ( PRESENT(dim_names_in) )                                            &
        dim_names_in = tile_dim_name
      IF ( PRESENT(dim_names_out) )                                           &
        dim_names_out = TILE_DIM_NAME_OUT
      IF ( PRESENT(dim_sizes) )                                               &
        dim_sizes = tile_dim_size

    CASE DEFAULT
! Unrecognised variable type
      CALL log_fatal("get_var_levs_dims",                                     &
                     "Unrecognised variable type - check metadata")
  END SELECT

  RETURN

END SUBROUTINE get_var_levs_dims
#endif
