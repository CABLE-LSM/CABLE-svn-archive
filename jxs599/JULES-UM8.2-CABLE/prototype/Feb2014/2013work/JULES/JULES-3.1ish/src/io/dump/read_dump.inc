#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/dump/read_dump.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE read_dump(file_name, identifiers)

  USE io_constants, ONLY : MODE_READ

  USE dictionary_mod, ONLY : dict, dict_create, dict_set, dict_get, dict_free

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_introspect,               &
                       file_inquire_dim, file_inquire_var, file_read_var,     &
                       file_close

  USE ancil_info, ONLY : dim_cs1, frac, land_pts, nsmax, ntiles, sm_levels

  USE nstypes, ONLY : npft, ntype

  USE imogen_constants, ONLY : n_olevs, nfarray

  USE prognostics, ONLY : canht_ft, canopy, cs, gs, lai, nsnow, rgrain,       &
                          rgrainL, rho_snow_grnd, sice, sliq, snow_grnd,      &
                          snow_tile, snowdepth, t_soil, tsnow, tstar_tile

  USE imogen_progs, ONLY : co2_ppmv, co2_change_ppmv, dtemp_o, fa_ocean,      &
                           seed_rain

  USE trifctl, ONLY : cv

  USE snow_param, ONLY : ds

  USE p_s_parms, ONLY : sthu, sthf

  USE top_pdm, ONLY : sthzw, zw

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Check that the given file is a JULES dump compatible with the current
!   run, and read the given identifiers from it if they are there
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  CHARACTER(len=*) :: file_name  ! The dump file
  CHARACTER(len=*) :: identifiers(:)  ! The model identifiers for the variables
                                      ! to define

! Work variables
  TYPE(file_handle) :: file  ! The opened file

  INTEGER :: nvars  ! The number of variables we are processing

  INTEGER :: dim_size_file  ! The size of the dimension currently being
                            ! processed in the file


  INTEGER :: dim_id  ! The current dimension id
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name  ! The current dimension name
  INTEGER :: dim_size  ! The current dimension size

  TYPE(dict) :: dims  ! The required dimensions for the file - maps dimension
                      ! names to dimension sizes

  INTEGER :: var_id  ! The current variable id

  LOGICAL :: is_record_dim  ! Dummy input for dimension inquiry function

! Real versions of integer valued variables
  REAL :: nsnow_real(land_pts, ntiles)
  REAL :: seed_rain_real(4)

  INTEGER :: i  ! Loop counter


!-----------------------------------------------------------------------------

  nvars = SIZE(identifiers)

!-----------------------------------------------------------------------------
! We use the lowest level file API here, as we don't want to impose the input
! grid
!-----------------------------------------------------------------------------
  file = file_open(file_name, MODE_READ)

! We want to auto-detect the dimensions and variables in the file
  CALL file_introspect(file)

!-----------------------------------------------------------------------------
! Check the existence and size of dimensions of the requested variables
!-----------------------------------------------------------------------------
! First, collect information about the dimensions we need for the variables
! requested
! We collect this information in a dictionary that maps dimension name as key
! to dimension size as value (we can safely call dict_set with the same key
! more than once as the value will just be overwritten)
  dims = dict_create(MAX_DIM_DUMP, 1)

  DO i=1,nvars
    SELECT CASE ( identifiers(i) )
      CASE ( 'gs', 'sthzw', 'zw', 'cv' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)

      CASE ( 'canht', 'lai' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, PFT_DIM_NAME, npft)

      CASE ( 'cs' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, SC_POOL_DIM_NAME, dim_cs1)

      CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, TILE_DIM_NAME, ntiles)
        CALL dict_set(dims, SNOW_DIM_NAME, nsmax)

      CASE ( 'sthuf', 't_soil' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, SOIL_DIM_NAME, sm_levels)

      CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',            &
                                'snow_depth', 'snow_grnd', 'tstar_tile' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, TILE_DIM_NAME, ntiles)

      CASE ( 'frac' )
        CALL dict_set(dims, LAND_DIM_NAME, land_pts)
        CALL dict_set(dims, TYPE_DIM_NAME, ntype)

! Cases for IMOGEN variables
      CASE ( 'co2_ppmv', 'co2_change_ppmv' )
! scalar variables are represented by an array of dimension 1
        CALL dict_set(dims, SCALAR_DIM_NAME, 1)

      CASE ( 'dtemp_o' )
        CALL dict_set(dims, NOLEVS_DIM_NAME, n_olevs)

      CASE ( 'fa_ocean' )
        CALL dict_set(dims, NFARRAY_DIM_NAME, nfarray)

      CASE ( 'seed_rain' )
        CALL dict_set(dims, SEED_DIM_NAME, 4)

      CASE DEFAULT
        CALL log_fatal("read_dump",                                           &
                       "Unrecognised variable for dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! Check that the dimensions exist and have the correct size
  DO i = 1,dims%length
    dim_name = dims%keys(i)
    CALL dict_get(dims, dim_name, dim_size)

    CALL file_inquire_dim(file, dim_name, dim_id, dim_size_file, is_record_dim)

    IF ( dim_id < 0 )                                                         &
      CALL log_fatal("read_dump",                                             &
                     "Expected dimension '" // TRIM(dim_name) // "' in dump file")

    IF ( dim_size_file /= dim_size )                                          &
      CALL log_fatal("read_dump",                                             &
                     "Dimension '" // TRIM(dim_name) // "' has size not " //  &
                     "compatible with current run (required: " //             &
                     TRIM(to_string(dim_size)) // ", found: " //              &
                     TRIM(to_string(dim_size_file)) // ")")
  END DO

  CALL dict_free(dims)


!-----------------------------------------------------------------------------
! Set the requested variables from the file
! We assume that if the file passed all the checks on dimensions above, then
! it will be fine to fill variables here (i.e. we don't check the dimensions
! associated with the variables)
!-----------------------------------------------------------------------------
  DO i = 1,nvars
    CALL file_inquire_var(file, identifiers(i), var_id)

! Check that we managed to find the variable
    IF ( var_id < 1 )                                                         &
      CALL log_fatal("read_dump",                                             &
                     "Failed to find requested variable '" //                 &
                     TRIM(identifiers(i)) // "' in dump file")

! Since the dump variables have the exact same dimensions as their model
! counterparts, we can just pass the model variables into the read routines
    SELECT CASE ( identifiers(i) )
      CASE ( 'canopy' )
        CALL file_read_var(file, var_id, canopy)

      CASE ( 'cs' )
        CALL file_read_var(file, var_id, cs)

      CASE ( 'gs' )
        CALL file_read_var(file, var_id, gs)

      CASE ( 'snow_tile' )
        CALL file_read_var(file, var_id, snow_tile)

      CASE ( 'sthuf' )
! sthuf is held in sthu until it is processed
        CALL file_read_var(file, var_id, sthu)

      CASE ( 't_soil' )
        CALL file_read_var(file, var_id, t_soil)

      CASE ( 'tstar_tile' )
        CALL file_read_var(file, var_id, tstar_tile)

      CASE ( 'lai' )
        CALL file_read_var(file, var_id, lai)

      CASE ( 'canht' )
        CALL file_read_var(file, var_id, canht_ft)

      CASE ( 'frac' )
        CALL file_read_var(file, var_id, frac)

      CASE ( 'sthzw' )
        CALL file_read_var(file, var_id, sthzw)

      CASE ( 'zw' )
        CALL file_read_var(file, var_id, zw)

      CASE ( 'rgrain' )
        CALL file_read_var(file, var_id, rgrain)

      CASE ( 'cv' )
        CALL file_read_var(file, var_id, cv)

      CASE ( 'rho_snow' )
        CALL file_read_var(file, var_id, rho_snow_grnd)

      CASE ( 'snow_depth' )
        CALL file_read_var(file, var_id, snowDepth)

      CASE ( 'snow_grnd' )
        CALL file_read_var(file, var_id, snow_grnd)

      CASE ( 'nsnow' )
        CALL file_read_var(file, var_id, nsnow_real)
        nsnow = NINT(nsnow_real)

      CASE ( 'snow_ds' )
        CALL file_read_var(file, var_id, ds)

      CASE ( 'snow_ice' )
        CALL file_read_var(file, var_id, sice)

      CASE ( 'snow_liq' )
        CALL file_read_var(file, var_id, sliq)

      CASE ( 'tsnow' )
        CALL file_read_var(file, var_id, tsnow)

      CASE ( 'rgrainl' )
        CALL file_read_var(file, var_id, rgrainL)

! Cases for IMOGEN variables
      CASE ( 'co2_ppmv' )
        CALL file_read_var(file, var_id, co2_ppmv)

      CASE ( 'co2_change_ppmv' )
        CALL file_read_var(file, var_id, co2_change_ppmv)

      CASE ( 'dtemp_o' )
        CALL file_read_var(file, var_id, dtemp_o)

      CASE ( 'fa_ocean' )
        CALL file_read_var(file, var_id, fa_ocean)

      CASE ( 'seed_rain' )
        CALL file_read_var(file, var_id, seed_rain_real)
        seed_rain = NINT(seed_rain_real)

      CASE DEFAULT
        CALL log_fatal("read_dump",                                           &
                       "Unexpected variable in dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! We are done with the file
  CALL file_close(file)

  RETURN

END SUBROUTINE read_dump
#endif
