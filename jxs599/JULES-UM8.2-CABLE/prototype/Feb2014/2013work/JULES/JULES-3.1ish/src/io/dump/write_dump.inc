#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/dump/write_dump.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE write_dump()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_DIM_VAR, MODE_WRITE

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE dictionary_mod, ONLY : dict, dict_create, dict_set, dict_get, dict_free

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_def_dim, file_def_var,    &
                       file_enddef, file_write_var, file_close

  USE output_mod, ONLY : output_dir, run_id

  USE model_time_mod, ONLY : current_time, is_spinup, spinup_cycle

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
!   Writes a dump file for the current timestep
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file_name
                                      ! The filename to use for the dump file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: dt_string
                                      ! The datetime string to use in the file
                                      ! name
  CHARACTER(len=IDENTIFIER_LEN) :: identifiers(MAX_VAR_DUMP)
                                      ! The model identifiers for the variables
                                      ! to put in the dump

  TYPE(file_handle) :: file  ! The dump file

  INTEGER :: nvars  ! The number of variables we are processing

! Variables used when defining dimensions
  TYPE(dict) :: dims  ! The dimensions that we need to define
                      ! Indexed by dimension name (name => dim size)
  TYPE(dict) :: dim_ids  ! The ids of the the dimensions once they have been
                         ! defined
                         ! Indexed by dimension name (name => dim id)
  INTEGER :: dim_id  ! The current dimension id
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_name  ! The current dimension name
  INTEGER :: dim_size  ! The current dimension size

! Variables used when defining variables
  INTEGER :: var_ndims  ! The number of dimensions the current variable has
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_dims(MAX_DIM_VAR)
                        ! The dimension names for the current variable
  INTEGER :: var_dim_ids(MAX_DIM_VAR)
                        ! The id in file of each dimension for the current
                        ! variable

  TYPE(dict) :: var_ids  ! The ids of the variables once they have been defined
                         ! Indexed by identifier (identifier => var id)
  INTEGER :: var_id  ! The current variable id

  INTEGER :: i,j  ! Loop counters


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get the list of identifiers that we are going to output
!-----------------------------------------------------------------------------
  CALL required_vars_for_configuration(nvars, identifiers)

!-----------------------------------------------------------------------------
! Generate the file name that we want to use and open the file
!-----------------------------------------------------------------------------
! File name starts with run id + indicator of a dump file
  file_name = TRIM(run_id) // ".dump."

! Include the current spinup cycle if there is one
  IF ( is_spinup )                                                            &
    file_name = TRIM(file_name) // "spin" // TRIM(to_string(spinup_cycle)) // "."

! Then current date and time
  WRITE(dt_string, '(I4.4,I2.2,I2.2)') current_time%year, current_time%month, &
                                       current_time%day
  dt_string = TRIM(dt_string) // "." // TRIM(to_string(current_time%time))
  file_name = TRIM(file_name) // TRIM(dt_string)

! Add the extension based on dump format
  SELECT CASE ( DUMP_FORMAT )
    CASE ( FORMAT_ASCII )
      file_name = TRIM(file_name) // ".asc"

    CASE ( FORMAT_NCDF )
      file_name = TRIM(file_name) // ".nc"

    CASE DEFAULT
      CALL log_fatal("write_dump",                                            &
                     "Unrecognised file format - " // TRIM(DUMP_FORMAT))
  END SELECT

! Prepend the output directory
  file_name = TRIM(output_dir) // "/" // TRIM(file_name)

! We use the lowest level file API here, as we don't want to impose a grid
  file = file_open(file_name, MODE_WRITE)

!-----------------------------------------------------------------------------
! Check what dimensions we need to create based on the variables that will
! be in the dump
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
        CALL log_fatal("write_dump",                                          &
                       "Unrecognised variable for dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! Define the dimensions and gather up their ids
  dim_ids = dict_create(MAX_DIM_DUMP, 1)
  DO i = 1,dims%length
    dim_name = dims%keys(i)
    CALL dict_get(dims, dim_name, dim_size)

! None of the dimensions in dump files are record dimensions
    dim_id = file_def_dim(file, dim_name, dim_size)
    CALL dict_set(dim_ids, dim_name, dim_id)
  END DO


!-----------------------------------------------------------------------------
! Define the variables that will be in the dump
!-----------------------------------------------------------------------------
  var_ids = dict_create(MAX_VAR_DUMP, 1)
  DO i = 1,nvars
    SELECT CASE ( identifiers(i) )
      CASE ( 'gs', 'sthzw', 'zw', 'cv' )
        var_ndims = 1
        var_dims(1) = LAND_DIM_NAME

      CASE ( 'canht', 'lai' )
        var_ndims = 2
        var_dims(1:2) = (/ LAND_DIM_NAME, PFT_DIM_NAME /)

      CASE ( 'cs' )
        var_ndims = 2
        var_dims(1:2) = (/ LAND_DIM_NAME, SC_POOL_DIM_NAME /)

      CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
        var_ndims = 3
        var_dims(1:3) = (/ LAND_DIM_NAME, TILE_DIM_NAME, SNOW_DIM_NAME /)

      CASE ( 'sthuf', 't_soil' )
        var_ndims = 2
        var_dims(1:2) = (/ LAND_DIM_NAME, SOIL_DIM_NAME /)

      CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',            &
                                'snow_depth', 'snow_grnd', 'tstar_tile' )
        var_ndims = 2
        var_dims(1:2) = (/ LAND_DIM_NAME, TILE_DIM_NAME /)

      CASE ( 'frac' )
        var_ndims = 2
        var_dims(1:2) = (/ LAND_DIM_NAME, TYPE_DIM_NAME /)

! Cases for IMOGEN variables
      CASE ( 'co2_ppmv', 'co2_change_ppmv' )
! scalars are represented by an array of dimension 1
        var_ndims = 1
        var_dims(1) = SCALAR_DIM_NAME

      CASE ( 'dtemp_o' )
        var_ndims = 1
        var_dims(1) = NOLEVS_DIM_NAME

      CASE ( 'fa_ocean' )
        var_ndims = 1
        var_dims(1) = NFARRAY_DIM_NAME

      CASE ( 'seed_rain' )
        var_ndims = 1
        var_dims(1) = SEED_DIM_NAME

      CASE DEFAULT
        CALL log_fatal("write_dump",                                          &
                       "Unrecognised variable for dump - " // TRIM(identifiers(i)))
    END SELECT

! Convert the dimension names into dimension ids for file_def_var
    DO j = 1,var_ndims
      CALL dict_get(dim_ids, var_dims(j), var_dim_ids(j))
    END DO

! Actually define the variable and store its id
! No variables use the record dimension, since there isn't one
! We use the identifier as the variable name in the file
    var_id = file_def_var(file, identifiers(i), var_dim_ids(1:var_ndims), .FALSE.)
    CALL dict_set(var_ids, identifiers(i), var_id)
  END DO

! We have finished defining things
  CALL file_enddef(file)


!-----------------------------------------------------------------------------
! Write data to file
! We can pass the model variables to write var since the variables in file
! have the exact same dimensions
!-----------------------------------------------------------------------------
  DO i = 1,nvars
    CALL dict_get(var_ids, identifiers(i), var_id)

    SELECT CASE ( identifiers(i) )
      CASE ( 'canopy' )
        CALL file_write_var(file, var_id, canopy)

      CASE ( 'cs' )
        CALL file_write_var(file, var_id, cs)

      CASE ( 'gs' )
        CALL file_write_var(file, var_id, gs)

      CASE ( 'snow_tile' )
        CALL file_write_var(file, var_id, snow_tile)

      CASE ( 'sthuf' )
! We dump the total wetness
        CALL file_write_var(file, var_id, sthu + sthf)

      CASE ( 't_soil' )
        CALL file_write_var(file, var_id, t_soil)

      CASE ( 'tstar_tile' )
        CALL file_write_var(file, var_id, tstar_tile)

      CASE ( 'lai' )
        CALL file_write_var(file, var_id, lai)

      CASE ( 'canht' )
        CALL file_write_var(file, var_id, canht_ft)

      CASE ( 'frac' )
        CALL file_write_var(file, var_id, frac)

      CASE ( 'sthzw' )
        CALL file_write_var(file, var_id, sthzw)

      CASE ( 'zw' )
        CALL file_write_var(file, var_id, zw)

      CASE ( 'rgrain' )
        CALL file_write_var(file, var_id, rgrain)

      CASE ( 'cv' )
        CALL file_write_var(file, var_id, cv)

      CASE ( 'rho_snow' )
        CALL file_write_var(file, var_id, rho_snow_grnd)

      CASE ( 'snow_depth' )
        CALL file_write_var(file, var_id, snowDepth)

      CASE ( 'snow_grnd' )
        CALL file_write_var(file, var_id, snow_grnd)

      CASE ( 'nsnow' )
        CALL file_write_var(file, var_id, REAL(nsnow))

      CASE ( 'snow_ds' )
        CALL file_write_var(file, var_id, ds)

      CASE ( 'snow_ice' )
        CALL file_write_var(file, var_id, sice)

      CASE ( 'snow_liq' )
        CALL file_write_var(file, var_id, sliq)

      CASE ( 'tsnow' )
        CALL file_write_var(file, var_id, tsnow)

      CASE ( 'rgrainl' )
        CALL file_write_var(file, var_id, rgrainL)

! Cases for IMOGEN variables
      CASE ( 'co2_ppmv' )
        CALL file_write_var(file, var_id, co2_ppmv)

      CASE ( 'co2_change_ppmv' )
        CALL file_write_var(file, var_id, co2_change_ppmv)

      CASE ( 'dtemp_o' )
        CALL file_write_var(file, var_id, dtemp_o)

      CASE ( 'fa_ocean' )
        CALL file_write_var(file, var_id, fa_ocean)

      CASE ( 'seed_rain' )
        CALL file_write_var(file, var_id, REAL(seed_rain))

      CASE DEFAULT
        CALL log_fatal("write_dump",                                          &
                       "Unrecognised variable for dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! We are done with the file and dictionaries
  CALL file_close(file)
  CALL dict_free(dims)
  CALL dict_free(dim_ids)
  CALL dict_free(var_ids)

  RETURN

END SUBROUTINE write_dump
#endif
