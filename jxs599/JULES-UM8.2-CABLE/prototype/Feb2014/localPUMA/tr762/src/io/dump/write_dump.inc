#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE write_dump()

  USE io_constants, ONLY : MAX_FILE_NAME_LEN, MAX_DIM_VAR, MODE_WRITE

  USE parallel_mod, ONLY : MASTER_TASK_ID, is_master_task, gather_land_field

  USE model_interface_mod, ONLY : IDENTIFIER_LEN

  USE dictionary_mod, ONLY : dict, dict_create, dict_set, dict_get,           &
                             dict_has_key, dict_free

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_def_dim, file_def_var,    &
                       file_enddef, file_write_var, file_close

  USE output_mod, ONLY : output_dir, run_id

  USE model_time_mod, ONLY : current_time, is_spinup, spinup_cycle

  USE model_grid_mod, ONLY : global_land_pts

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
!   Note that the writing of the dump is done by the master task with the
!   values gathered from other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
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
  TYPE(dict) :: file_dim_ids  ! Dictionary of the dimensions that have been defined
                              ! Maps dim_name => dim_id

  INTEGER :: ndims  ! The number of levels dims for the current variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_VAR)
                    ! The levels dimension names for the current variable
  INTEGER :: dim_sizes(MAX_DIM_VAR)
                    ! The sizes of the levels dims for the current variable
  INTEGER :: dim_ids(MAX_DIM_VAR)
                   ! The ids in file of the levels dims for the current variable

  INTEGER :: var_ids(MAX_VAR_DUMP)
                        ! The ids of the variables in the dump file

  INTEGER :: i, j, m, n  ! Loop counters


! Arrays to hold global land points version of data gathered in master task
! before writing
  REAL, ALLOCATABLE :: global_data(:)  ! For data with no vertical levels
  REAL, ALLOCATABLE :: global_data_levs(:,:)   ! With one vertical level
  REAL, ALLOCATABLE :: global_data_snow(:,:,:) ! With two "vertical" levels
                                               ! I.E. snow variables


!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Get the list of identifiers that we are going to output
!-----------------------------------------------------------------------------
  CALL required_vars_for_configuration(nvars, identifiers)

!-----------------------------------------------------------------------------
! In the master task only, we open a new file and define the required
! dimensions and variables
!-----------------------------------------------------------------------------
  IF ( is_master_task() ) THEN
!-----------------------------------------------------------------------------
! Generate the file name that we want to use and open the file
!-----------------------------------------------------------------------------
! File name starts with run id + indicator of a dump file
    file_name = TRIM(run_id) // ".dump."

! Include the current spinup cycle if there is one
    IF ( is_spinup )                                                          &
      file_name = TRIM(file_name) // "spin" // TRIM(to_string(spinup_cycle)) // "."

! Then current date and time
    WRITE(dt_string, '(I4.4,I2.2,I2.2)') current_time%year,                   &
                                         current_time%month,                  &
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
        CALL log_fatal("write_dump",                                          &
                       "Unrecognised file format - " // TRIM(DUMP_FORMAT))
    END SELECT

! Prepend the output directory
    file_name = TRIM(output_dir) // "/" // TRIM(file_name)

! We use the lowest level file API here, as we don't want to impose a grid
    file = file_open(file_name, MODE_WRITE)

!-----------------------------------------------------------------------------
! Create the dimensions and variables
!-----------------------------------------------------------------------------
    file_dim_ids = dict_create(MAX_DIM_DUMP, INT(1))

    DO i = 1,nvars

!-----------------------------------------------------------------------------
! Get information about the dimensions used by the variable
!-----------------------------------------------------------------------------
      SELECT CASE ( identifiers(i) )
        CASE ( 'gs', 'sthzw', 'zw', 'cv' )
          ndims = 1
          dim_names(1) = LAND_DIM_NAME
          dim_sizes(1) = global_land_pts

        CASE ( 'canht', 'lai' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, PFT_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, npft /)

        CASE ( 'cs' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, SC_POOL_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, dim_cs1 /)

        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          ndims = 3
          dim_names(1:3) = (/ LAND_DIM_NAME, TILE_DIM_NAME, SNOW_DIM_NAME /)
          dim_sizes(1:3) = (/ global_land_pts, ntiles, nsmax /)

        CASE ( 'sthuf', 't_soil' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, SOIL_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, sm_levels /)

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
                                  'snow_depth', 'snow_grnd', 'tstar_tile' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, TILE_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, ntiles /)

        CASE ( 'frac' )
          ndims = 2
          dim_names(1:2) = (/ LAND_DIM_NAME, TYPE_DIM_NAME /)
          dim_sizes(1:2) = (/ global_land_pts, ntype /)

! Cases for IMOGEN variables
        CASE ( 'co2_ppmv', 'co2_change_ppmv' )
! scalar variables are represented by an array of dimension 1
          ndims = 1
          dim_names(1) = SCALAR_DIM_NAME
          dim_sizes(1) = 1

        CASE ( 'dtemp_o' )
          ndims = 1
          dim_names(1) = NOLEVS_DIM_NAME
          dim_sizes(1) = n_olevs

        CASE ( 'fa_ocean' )
          ndims = 1
          dim_names(1) = NFARRAY_DIM_NAME
          dim_sizes(1) = nfarray

        CASE ( 'seed_rain' )
          ndims = 1
          dim_names(1) = SEED_DIM_NAME
          dim_sizes(1) = 4

        CASE DEFAULT
          CALL log_fatal("write_dump",                                        &
                         "Unrecognised variable for dump - " // TRIM(identifiers(i)))
      END SELECT

!-----------------------------------------------------------------------------
! Define the dimensions if they have not already been defined
! We use a dictionary to keep track of defined dimension ids
!
! At the same time, gather up the dimension ids needed by the current variable
!-----------------------------------------------------------------------------
      DO j = 1,ndims
! If it has not yet been defined, define the dimension, storing its id
        IF ( .NOT. dict_has_key(file_dim_ids, dim_names(j)) )                 &
          CALL dict_set(                                                      &
            file_dim_ids, dim_names(j),                                       &
            file_def_dim(file, dim_names(j), dim_sizes(j))                    &
          )

! Get the dimension id from the dict and add it to the list for this variable
        CALL dict_get(file_dim_ids, dim_names(j), dim_ids(j))
      END DO

!-----------------------------------------------------------------------------
! Define the variable, saving the id in the file for later
!-----------------------------------------------------------------------------
      var_ids(i) = file_def_var(file, identifiers(i), dim_ids(1:ndims), .FALSE.)

    END DO

!-----------------------------------------------------------------------------
! We have finished defining things
!-----------------------------------------------------------------------------
    CALL file_enddef(file)
    CALL dict_free(file_dim_ids)

  END IF  ! MASTER TASK


!-----------------------------------------------------------------------------
! Gather data from other tasks and write it to file
!-----------------------------------------------------------------------------
! Allocate the global data arrays
  IF ( is_master_task() ) THEN
    ALLOCATE(global_data(global_land_pts))
    ALLOCATE(global_data_levs(global_land_pts, MAX(npft,dim_cs1,sm_levels,ntiles,ntype)))
    ALLOCATE(global_data_snow(global_land_pts, ntiles, nsmax))
  ELSE
    ALLOCATE(global_data(1))
    ALLOCATE(global_data_levs(1,1))
    ALLOCATE(global_data_snow(1,1,1))
  END IF

  DO i = 1,nvars
! Gather the variables into a global array to write to file
! Note that gather_land_field can only gather one land_pts array at a time,
! so to gather variables with multiple levels we must loop
    SELECT CASE ( identifiers(i) )
      CASE ( 'gs' )
        CALL gather_land_field(gs, global_data)

      CASE ( 'sthzw' )
        CALL gather_land_field(sthzw, global_data)

      CASE ( 'zw' )
        CALL gather_land_field(zw, global_data)

      CASE ( 'cv' )
        CALL gather_land_field(cv, global_data)

      CASE ( 'canht' )
        DO n = 1,npft
          CALL gather_land_field(canht_ft(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'lai' )
        DO n = 1,npft
          CALL gather_land_field(lai(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'cs' )
        DO n = 1,dim_cs1
          CALL gather_land_field(cs(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'sthuf' )
! sthuf is held in sthu until it is processed
        DO n = 1,sm_levels
          CALL gather_land_field(sthu(:,n) + sthf(:,n), global_data_levs(:,n))
        END DO

      CASE ( 't_soil' )
        DO n = 1,sm_levels
          CALL gather_land_field(t_soil(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'canopy' )
        DO n = 1,ntiles
          CALL gather_land_field(canopy(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'nsnow' )
        DO n = 1,ntiles
          CALL gather_land_field(REAL(nsnow(:,n)), global_data_levs(:,n))
        END DO

      CASE ( 'rgrain' )
        DO n = 1,ntiles
          CALL gather_land_field(rgrain(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'rho_snow' )
        DO n = 1,ntiles
          CALL gather_land_field(rho_snow_grnd(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'snow_tile' )
        DO n = 1,ntiles
          CALL gather_land_field(snow_tile(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'snow_depth' )
        DO n = 1,ntiles
          CALL gather_land_field(snowDepth(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'snow_grnd' )
        DO n = 1,ntiles
          CALL gather_land_field(snow_grnd(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'tstar_tile' )
        DO n = 1,ntiles
          CALL gather_land_field(tstar_tile(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'frac' )
        DO n = 1,ntype
          CALL gather_land_field(frac(:,n), global_data_levs(:,n))
        END DO

      CASE ( 'rgrainl' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL gather_land_field(rgrainL(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ds' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL gather_land_field(ds(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ice' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL gather_land_field(sice(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO


      CASE ( 'snow_liq' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL gather_land_field(sliq(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO


      CASE ( 'tsnow' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL gather_land_field(tsnow(:,m,n), global_data_snow(:,m,n))
          END DO
        END DO

! Since each task runs its own version of IMOGEN, we just use the values from
! the master task
      CASE ( 'co2_ppmv', 'co2_change_ppmv', 'dtemp_o', 'fa_ocean', 'seed_rain' )
! Nothing to do

      CASE DEFAULT
        CALL log_fatal("write_dump",                                          &
                       "No code to gather variable for dump - " // TRIM(identifiers(i)))
    END SELECT

!-----------------------------------------------------------------------------
! In the master task, write the global data to file
!-----------------------------------------------------------------------------
    IF ( is_master_task() ) THEN

      SELECT CASE ( identifiers(i) )
! If it is a land_pts array with no levels associated, write the global_data array
        CASE ( 'gs', 'sthzw', 'zw', 'cv' )
          CALL file_write_var(file, var_ids(i), global_data)

! If it is a variable with one vertical level, write the appropriate number of
! levels from global_data_levs
        CASE ( 'canht', 'lai' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:npft))

        CASE ( 'cs' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:dim_cs1))

        CASE ( 'sthuf', 't_soil' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:sm_levels))

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
                                  'snow_depth', 'snow_grnd', 'tstar_tile' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:ntiles))

        CASE ( 'frac' )
          CALL file_write_var(file, var_ids(i), global_data_levs(:,1:ntype))

! If it is a snow variable, first allocate the array if required, then read it
        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          CALL file_write_var(file, var_ids(i), global_data_snow)

! Cases for IMOGEN variables
! Each task runs its own version of IMOGEN - we just write the master task's
! versions
        CASE ( 'co2_ppmv' )
          CALL file_write_var(file, var_ids(i), co2_ppmv)

        CASE ( 'co2_change_ppmv' )
          CALL file_write_var(file, var_ids(i), co2_change_ppmv)

        CASE ( 'dtemp_o' )
          CALL file_write_var(file, var_ids(i), dtemp_o)

        CASE ( 'fa_ocean' )
          CALL file_write_var(file, var_ids(i), fa_ocean)

        CASE ( 'seed_rain' )
          CALL file_write_var(file, var_ids(i), REAL(seed_rain))

        CASE DEFAULT
          CALL log_fatal("write_dump",                                        &
                         "Unrecognised variable for dump - " // TRIM(identifiers(i)))
      END SELECT
    END IF  ! MASTER TASK
  END DO

! We are done with the file and dictionaries
  IF ( is_master_task() ) CALL file_close(file)

  IF ( ALLOCATED(global_data) ) DEALLOCATE(global_data)
  IF ( ALLOCATED(global_data_levs) ) DEALLOCATE(global_data_levs)
  IF ( ALLOCATED(global_data_snow) ) DEALLOCATE(global_data_snow)

  RETURN

END SUBROUTINE write_dump
#endif
