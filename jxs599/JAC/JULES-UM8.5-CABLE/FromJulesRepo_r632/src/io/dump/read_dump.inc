#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE read_dump(file_name, identifiers)

  USE io_constants, ONLY : MODE_READ, MAX_DIM_VAR

  USE parallel_mod, ONLY : MASTER_TASK_ID, is_master_task, scatter_land_field

  USE string_utils_mod, ONLY : to_string

  USE file_mod, ONLY : file_handle, file_open, file_introspect,               &
                       file_inquire_dim, file_inquire_var, file_read_var,     &
                       file_close

  USE model_grid_mod, ONLY : global_land_pts

  USE ancil_info, ONLY : land_pts, dim_cs1, frac, nsmax, ntiles, sm_levels

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

  INCLUDE "mpif.h"

!-----------------------------------------------------------------------------
! Description:
!   Check that the given file is a JULES dump compatible with the current
!   run, and read the given identifiers from it.
!   Note that the reading of the dump is done by the master task and the
!   results scattered to other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
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

! Used when defining dimensions and variables
  INTEGER :: ndims
                            ! The number of dimensions the current variable has
  CHARACTER(len=MAX_SDF_NAME_LEN) :: dim_names(MAX_DIM_VAR)
                            ! The dimension names the current variable should use
  INTEGER :: dim_sizes(MAX_DIM_VAR)
                            ! The dimension sizes for the current variable

  INTEGER :: dim_ids(MAX_DIM_VAR)  ! The dimension ids for the current variable as
                                   ! calculated from file_inquire_dim
  LOGICAL :: is_record_dim  ! Detects if the current dimension is a record dim

  INTEGER :: ndims_file  ! The number of dimensions the variable has in file
                         ! Compared to ndims above for each variable
  INTEGER :: dim_ids_file(MAX_DIM_VAR)
                         ! The ids of the dimensions the variable has in the file
                         ! Compared to dim_ids above to verify the variable has the
                         ! correct dimensions

  INTEGER :: var_ids(SIZE(identifiers))  ! The ids of the variables in the dump file

  LOGICAL :: is_record_var  ! Indicates if a variable uses the record dimension


! Real versions of integer valued variables
  REAL :: nsnow_real(land_pts, ntiles)
  REAL :: seed_rain_real(4)

  INTEGER :: i, j, m, n  ! Loop counters

  INTEGER :: error  ! Variable to collect MPI errors - most MPI implementations
                    ! bail on error, so this is not checked


! Arrays to hold global land points version of data read in master task
! before scattering
  REAL, ALLOCATABLE :: global_data(:)  ! For data with no vertical levels
  REAL, ALLOCATABLE :: global_data_levs(:,:)   ! With one vertical level
  REAL, ALLOCATABLE :: global_data_snow(:,:,:) ! With two "vertical" levels
                                               ! I.E. snow variables


!-----------------------------------------------------------------------------


  nvars = SIZE(identifiers)

!-----------------------------------------------------------------------------
! In the master task only, we open the file and check that the correct
! dimensions exist and are of a size compatible with this run
!-----------------------------------------------------------------------------
  IF ( is_master_task() ) THEN
!-----------------------------------------------------------------------------
! We use the lowest level file API here, as we don't want to impose the input
! grid
!-----------------------------------------------------------------------------
    file = file_open(file_name, MODE_READ)

! We want to auto-detect the dimensions and variables in the file
    CALL file_introspect(file)

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
          CALL log_fatal("read_dump",                                         &
                         "Unrecognised variable for dump - " // TRIM(identifiers(i)))
      END SELECT

!-----------------------------------------------------------------------------
! Check the dimensions exist and have the correct size
!-----------------------------------------------------------------------------
      DO j = 1,ndims
! Retrive information about the dimension from the file we store the id
! for use outside this loop
        CALL file_inquire_dim(                                                &
          file, dim_names(j), dim_ids(j), dim_size_file, is_record_dim        &
        )

! Check that we found a dimension
        IF ( dim_ids(j) < 0 )                                                 &
          CALL log_fatal("read_dump",                                         &
                         "Could not find expected dimension '" //             &
                         TRIM(dim_names(j)) // "' in dump file")

! Check that the dimension is not a record dimension (there shouldn't be one
! in dump files)
        IF ( is_record_dim )                                                  &
          CALL log_fatal("read_dump",                                         &
                         "Dimension '" // TRIM(dim_names(j)) // "' is a " //  &
                         "record dimension - should not exist in dump file")

! Check that the dimension has the correct size
        IF ( dim_size_file /= dim_sizes(j) )                                  &
          CALL log_fatal("read_dump",                                         &
                         "Dimension '" // TRIM(dim_names(j)) // "' has " //   &
                         "size incompatible with current run (required: " //  &
                         TRIM(to_string(dim_sizes(j))) // ", found: " //      &
                         TRIM(to_string(dim_size_file)) // ")")
      END DO  ! dims

!-----------------------------------------------------------------------------
! Check that the variable exists and has the correct dimensions
!-----------------------------------------------------------------------------
! Retrieve information about the variable from the file
      CALL file_inquire_var(                                                  &
        file, identifiers(i), var_ids(i), ndims_file, dim_ids_file,           &
        is_record_var                                                         &
      )

! Check that we found a variable
      IF ( var_ids(i) < 1 )                                                   &
        CALL log_fatal("read_dump",                                           &
                       "Failed to find requested variable '" //               &
                       TRIM(identifiers(i)) // "' in dump file")

! Check that the number of dimensions match
      IF ( ndims_file /= ndims )                                              &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect number of dimensions in dump file (" //     &
                       "expected: " // TRIM(to_string(ndims)) // ", " //      &
                       "found: " // TRIM(to_string(ndims_file)) // ")")

! Check that the dimension ids match
      IF ( .NOT. ALL(dim_ids(1:ndims) == dim_ids_file(1:ndims)) )             &
        CALL log_fatal("read_dump",                                           &
                       "Variable '" // TRIM(identifiers(i)) // "' has " //    &
                       "incorrect dimensions in dump file")

    END DO  ! vars

  END IF  ! MASTER TASK


!-----------------------------------------------------------------------------
! Set the requested variables from the file
!
! This is done by reading the value of the variable on global land points
! in the master task, then scattering it to the other tasks
!
! We assume that if the file passed all the checks on dimensions above, then
! it will be fine to fill variables here (i.e. we don't check the dimensions
! associated with the variables)
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
!-----------------------------------------------------------------------------
! In the master task, read the global data
!-----------------------------------------------------------------------------
    IF ( is_master_task() ) THEN
      SELECT CASE ( identifiers(i) )
! If it is a land_pts array with no levels associated, read into the
! global_data array
        CASE ( 'gs', 'sthzw', 'zw', 'cv' )
          CALL file_read_var(file, var_ids(i), global_data)

! If it is a variable with one vertical level, read the appropriate number of
! levels into global_data_levs
        CASE ( 'canht', 'lai' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:npft))

        CASE ( 'cs' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:dim_cs1))

        CASE ( 'sthuf', 't_soil' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:sm_levels))

        CASE ( 'canopy', 'nsnow', 'rgrain', 'rho_snow', 'snow_tile',          &
                                  'snow_depth', 'snow_grnd', 'tstar_tile' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:ntiles))

        CASE ( 'frac' )
          CALL file_read_var(file, var_ids(i), global_data_levs(:,1:ntype))

! If it is a snow variable, first allocate the array if required, then read it
        CASE ( 'rgrainl', 'snow_ds', 'snow_ice', 'snow_liq', 'tsnow' )
          CALL file_read_var(file, var_ids(i), global_data_snow)

! Cases for IMOGEN variables
! Each task runs its own version of IMOGEN - these variables are broadcast to
! all tasks below
        CASE ( 'co2_ppmv' )
          CALL file_read_var(file, var_ids(i), co2_ppmv)

        CASE ( 'co2_change_ppmv' )
          CALL file_read_var(file, var_ids(i), co2_change_ppmv)

        CASE ( 'dtemp_o' )
          CALL file_read_var(file, var_ids(i), dtemp_o)

        CASE ( 'fa_ocean' )
          CALL file_read_var(file, var_ids(i), fa_ocean)

        CASE ( 'seed_rain' )
          CALL file_read_var(file, var_ids(i), seed_rain_real)

        CASE DEFAULT
          CALL log_fatal("read_dump",                                         &
                         "Unexpected variable in dump - " // TRIM(identifiers(i)))
      END SELECT
    END IF  ! MASTER TASK

! Now scatter the variables into their final destinations
! Note that scatter_land_field can only scatter one land_pts array at a time,
! so to scatter variables with multiple levels we must loop
    SELECT CASE ( identifiers(i) )
      CASE ( 'gs' )
        CALL scatter_land_field(global_data, gs)

      CASE ( 'sthzw' )
        CALL scatter_land_field(global_data, sthzw)

      CASE ( 'zw' )
        CALL scatter_land_field(global_data, zw)

      CASE ( 'cv' )
        CALL scatter_land_field(global_data, cv)

      CASE ( 'canht' )
        DO n = 1,npft
          CALL scatter_land_field(global_data_levs(:,n), canht_ft(:,n))
        END DO

      CASE ( 'lai' )
        DO n = 1,npft
          CALL scatter_land_field(global_data_levs(:,n), lai(:,n))
        END DO

      CASE ( 'cs' )
        DO n = 1,dim_cs1
          CALL scatter_land_field(global_data_levs(:,n), cs(:,n))
        END DO

      CASE ( 'sthuf' )
! sthuf is held in sthu until it is processed
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_levs(:,n), sthu(:,n))
        END DO

      CASE ( 't_soil' )
        DO n = 1,sm_levels
          CALL scatter_land_field(global_data_levs(:,n), t_soil(:,n))
        END DO

      CASE ( 'canopy' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), canopy(:,n))
        END DO

      CASE ( 'nsnow' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), nsnow_real(:,n))
        END DO
        nsnow = NINT(nsnow_real)

      CASE ( 'rgrain' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), rgrain(:,n))
        END DO

      CASE ( 'rho_snow' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), rho_snow_grnd(:,n))
        END DO

      CASE ( 'snow_tile' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), snow_tile(:,n))
        END DO

      CASE ( 'snow_depth' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), snowDepth(:,n))
        END DO

      CASE ( 'snow_grnd' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), snow_grnd(:,n))
        END DO

      CASE ( 'tstar_tile' )
        DO n = 1,ntiles
          CALL scatter_land_field(global_data_levs(:,n), tstar_tile(:,n))
        END DO

      CASE ( 'frac' )
        DO n = 1,ntype
          CALL scatter_land_field(global_data_levs(:,n), frac(:,n))
        END DO

      CASE ( 'rgrainl' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL scatter_land_field(global_data_snow(:,m,n), rgrainL(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ds' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL scatter_land_field(global_data_snow(:,m,n), ds(:,m,n))
          END DO
        END DO

      CASE ( 'snow_ice' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL scatter_land_field(global_data_snow(:,m,n), sice(:,m,n))
          END DO
        END DO


      CASE ( 'snow_liq' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL scatter_land_field(global_data_snow(:,m,n), sliq(:,m,n))
          END DO
        END DO


      CASE ( 'tsnow' )
        DO n = 1,nsmax
          DO m = 1,ntiles
            CALL scatter_land_field(global_data_snow(:,m,n), tsnow(:,m,n))
          END DO
        END DO

! IMOGEN variables are just broadcast to all tasks
      CASE ( 'co2_ppmv' )
        CALL MPI_BCAST(co2_ppmv, 1, MPI_REAL,                                 &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'co2_change_ppmv' )
        CALL MPI_BCAST(co2_change_ppmv, 1, MPI_REAL,                          &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'dtemp_o' )
        CALL MPI_BCAST(dtemp_o, n_olevs, MPI_REAL,                            &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'fa_ocean' )
        CALL MPI_BCAST(fa_ocean, nfarray, MPI_REAL,                           &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)

      CASE ( 'seed_rain' )
        CALL MPI_BCAST(seed_rain_real, 4, MPI_REAL,                           &
                       MASTER_TASK_ID, MPI_COMM_WORLD, error)
        seed_rain = NINT(seed_rain_real)

      CASE DEFAULT
        CALL log_fatal("read_dump",                                           &
                       "Unexpected variable in dump - " // TRIM(identifiers(i)))
    END SELECT
  END DO

! We are done with the file
  IF ( is_master_task() ) CALL file_close(file)

  IF ( ALLOCATED(global_data) ) DEALLOCATE(global_data)
  IF ( ALLOCATED(global_data_levs) ) DEALLOCATE(global_data_levs)
  IF ( ALLOCATED(global_data_snow) ) DEALLOCATE(global_data_snow)

  RETURN

END SUBROUTINE read_dump
#endif
