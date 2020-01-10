#if !defined(UM_JULES)

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads in JULES_IRRIG namelist if necessary

SUBROUTINE init_irrig

USE io_constants, ONLY: max_file_name_len, max_sdf_name_len, namelist_unit

USE missing_data_mod, ONLY:                                                   &
!  imported scalar parameters
     imdi, rmdi

USE model_interface_mod, ONLY: identifier_len, populate_var, get_var_id

USE string_utils_mod, ONLY: to_string

USE jules_vegetation_mod, ONLY: l_irrig_dmd

USE crop_vars_mod, ONLY: irrtiles, frac_irr_all_tiles, nirrtile

USE max_dimensions, ONLY: npft_max

USE dump_mod, ONLY: ancil_dump_read

USE jules_surface_types_mod, ONLY: npft

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!    Reads in irrigation fraction
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
INTEGER :: error  ! Error indicator

INTEGER :: i  ! Loop counter

!-----------------------------------------------------------------------------
! Definition of the jules_irrig namelist
!-----------------------------------------------------------------------------
LOGICAL :: read_from_dump
CHARACTER(LEN=max_file_name_len) :: FILE
                      ! The name of the file to use for irrigation fraction
                      ! if read_file = .false.
LOGICAL :: read_file
                      !   T - read from file
                      !   F - use a constant value for all points
CHARACTER(LEN=max_sdf_name_len) :: var_name
                      ! The name of the variable in the file
REAL(KIND=real_jlslsm) :: const_frac_irr
                      ! The constant value to use for each variable if
                      ! read_file = F for that variable
LOGICAL :: frac_irrig_all_tiles
                      ! T - assign irr fraction to all tiles
                      ! F - assign irr fraction to specfied tiles only
INTEGER :: irrigtiles(npft_max)


!work
INTEGER :: year, month, day, time

NAMELIST  / jules_irrig/ read_from_dump, frac_irrig_all_tiles, nirrtile,      &
                       irrigtiles, read_file, FILE, var_name, const_frac_irr


!-----------------------------------------------------------------------------

! Nothing to do if irrigation demand model is not on
IF ( .NOT. l_irrig_dmd ) RETURN

!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
read_from_dump = .FALSE.
read_file = .TRUE.  ! Default value
FILE=''      ! Empty file name
var_name  = ''      ! Empty variable name.
const_frac_irr       = rmdi     ! Missing data value.
frac_irrig_all_tiles = .TRUE.
irrigtiles(:)        = imdi     ! Missing data value.

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
CALL log_info("init_irrig", "Reading JULES_IRRIG namelist...")

READ(namelist_unit, NML = jules_irrig, IOSTAT = error)

IF ( error /= 0 )                                                             &
  CALL log_fatal("init_irrig",                                                &
                "Error reading namelist JULES_IRRIG " //                      &
                "(IOSTAT=" // TRIM(to_string(error)) // ")")

ancil_dump_read%irrig = read_from_dump

IF ( .NOT. ancil_dump_read%irrig) THEN !we read from the ancil file

  !---------------------------------------------------------------------------
  ! Determine if variable to be read from file
  !---------------------------------------------------------------------------
  IF ( read_file ) THEN
    CALL log_info("init_irrig", "'frac_irrig' will be read from file")

    ! If the variable will be filled from file, register it here
    IF ( LEN_TRIM(FILE) == 0 )                                                &
      CALL log_fatal("init_irrig", "No file name provided")

    CALL fill_variables_from_file(FILE, (/ 'frac_irrig' /), (/ var_name /))
  ELSE
    ! If the variable is being set as a constant, populate it here.
    ! First check that a value has been provided.
    IF ( ABS( const_frac_irr - rmdi ) < EPSILON(1.0) )                        &
      CALL log_fatal("init_irrig",                                            &
                     "No constant value provided for frac_irrig.")

    CALL log_info("init_irrig",                                               &
                  "'frac_irrig' will be set to a constant = " //              &
                  to_string(const_frac_irr))

    CALL populate_var(get_var_id('frac_irrig'), const_val = const_frac_irr)
  END IF

  !---------------------------------------------------------------------------
  ! Set frac_irr_all_tiles from namelist value
  !---------------------------------------------------------------------------
  frac_irr_all_tiles = frac_irrig_all_tiles

  !---------------------------------------------------------------------------
  ! Process pft names to be assigned irr fraction (if not all)
  !---------------------------------------------------------------------------
  IF ( .NOT. frac_irr_all_tiles ) THEN
    IF ( ( nirrtile > npft ) .OR. ( nirrtile < 1 ) ) THEN
      CALL log_info("init_irrig",                                             &
                    "nirrtile should be greater than 1 and less than or" //   &
                    "equal to the number of pfts.")
    END IF
    ! Copy pft indices from namelist into allocated array.
    ! First check that values were given.
    IF ( ANY( irrigtiles(1:nirrtile) == imdi ) )                              &
      CALL log_fatal("init_irrig",                                            &
                     "Insufficient values provided for irrigtiles.")
    irrtiles(1:nirrtile) = irrigtiles(1:nirrtile)
  END IF

ELSE !We read from the dump file
  CALL log_info("init_frac",                                                  &
                "frac will be read from the dump file.  " //                  &
                "Namelist values ignored")

END IF !.NOT. ancil_dump_read%frac

END SUBROUTINE init_irrig
#endif
