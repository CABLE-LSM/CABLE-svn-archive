#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/initial_conditions/init_ic.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_ic()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE water_constants_mod, ONLY : rho_water

  USE string_utils_mod, ONLY : to_string

  USE dump_mod, ONLY : MAX_VAR_DUMP, required_vars_for_configuration, read_dump

  USE input_mod, ONLY : fill_variables_from_file

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE switches, ONLY : l_top, l_triffid

  USE switches_urban, ONLY : l_urban2t

  USE nstypes, ONLY : npft, ice, urban, urban_canyon, urban_roof

  USE seed, ONLY : frac_min

  USE ancil_info, ONLY : land_pts, lice_pts, soil_pts, sm_levels, frac,       &
                         lice_index, soil_index, tile_pts, tile_index

  USE prognostics, ONLY : smcl, t_soil

  USE p_s_parms, ONLY : sthu, sthf, smvcst, b, sathh

  USE soil_param, ONLY : dzsoil

  USE urban_param, ONLY : wrr

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Sets up the initial conditions for the run
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER :: nvars_required  ! The number of variables that are
                             ! required in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_VAR_DUMP)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file   ! The number of variables that will be set
                          ! from the given file (template?)

  LOGICAL :: reset_done  ! Indicates if a reset of frac to frac_min was
                         ! performed

  REAL :: urban_fraction  ! Used in partitioning of urban fraction into
                          ! canyon and roof

  INTEGER :: i,j,l  ! Loop counters

  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_initial namelist
!-----------------------------------------------------------------------------
  LOGICAL :: dump_file  ! T - the given file is a dump file
                        ! F - the given file is not a dump file

  LOGICAL :: total_snow
                        ! Switch indicating how the snow model is initialised
                        !   T - only snow_tile is needed
                        !       If nsmax>0, the layer values are determined
                        !       by the model
                        !   F - all snow variables are supplied directly (what
                        !       these are depends on nsmax)
                        ! total_snow = F is only allowed for dump files

  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_VAR_DUMP)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_VAR_DUMP)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
                        ! Defaults to T for every variable
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_VAR_DUMP)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_VAR_DUMP)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(MAX_VAR_DUMP)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /jules_initial/ total_snow, dump_file, file,                       &
                           nvars, var, use_file, var_name, tpl_name, const_val


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nvars_required = 0
  nvars_file     = 0
  dump_file      = .FALSE.
  total_snow     = .TRUE.
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is to set every variable from file

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_ic", "Reading JULES_INITIAL namelist...")

  OPEN(NAMELIST_UNIT, FILE='initial_conditions.nml',                          &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ic",                                                 &
                   "Error opening namelist file initial_conditions.nml " //   &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  READ(NAMELIST_UNIT, nml=jules_initial, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ic",                                                 &
                   "Error reading namelist JULES_INITIAL " //                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_ic",                                                 &
                   "Error closing namelist file initial_conditions.nml " //   &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up initial conditions using namelist values
!-----------------------------------------------------------------------------
! First check the combination of options
  IF ( .NOT. dump_file .AND. .NOT. total_snow )                               &
    CALL log_fatal("init_ic",                                                 &
                   "total_snow = F is only allowed when initialising " //     &
                   "from a dump file")

! Set up the required variables - we get the list by calling a procedure in
! dump_mod that tells us the required prognostic variables for the current
! model configuration
! We indicate whether or not we are using total_snow, as this affects the
! required variables
! We also indicate that we do not want IMOGEN prognostics in the list, even
! if they are required, since they are initialised in init_imogen
  CALL required_vars_for_configuration(                                       &
    nvars_required, required_vars, total_snow, .FALSE.                        &
  )

! If we are initialising from a dump and no variables were specified, then
! we assume that all variables will be initialised from the dump file
  IF ( dump_file .AND. nvars < 1 ) THEN
    CALL log_info("init_ic",                                                  &
                  "No variables given - will attempt to initialise all " //   &
                  "required variables from specified dump file")
    nvars = nvars_required
    var(1:nvars) = required_vars(1:nvars)
! Every variable will use the file
    use_file(:) = .TRUE.
! We don't need to set var_name, tpl_name or const_val since they are never
! used in the case of a dump file anyway
  END IF

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_ic",                                               &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

!-----------------------------------------------------------------------------
! Check which variables we will be using and partition them into variables
! set to constant values and variables set from file
!-----------------------------------------------------------------------------
  DO i = 1,nvars
!-----------------------------------------------------------------------------
! If the variable is one of the required vars, then we will be using it
!-----------------------------------------------------------------------------
    IF ( ANY(required_vars(1:nvars_required) == var(i)) ) THEN
      IF ( use_file(i) ) THEN
        CALL log_info("init_ic",                                              &
                      "'" // TRIM(var(i)) // "' will be read from file")

! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
! Since nvars_file <= i (so we will not overwrite unprocessed values)
! and we do not need the values from these arrays for any non-file variables
! from now on, we can just compress them down onto variables that are in the file
        var(nvars_file) = var(i)
        var_name(nvars_file) = var_name(i)
        tpl_name(nvars_file) = tpl_name(i)
      ELSE
! If the variable is being set as a constant, just populate it here
        CALL log_info("init_ic",                                              &
                      "'" // TRIM(var(i)) // "' will be set to a constant")

        CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
      END IF
    ELSE
! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_ic",                                                &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
  END DO

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    IF ( dump_file ) THEN
! If we are using a dump file, use read_dump to fill the variables
      CALL read_dump(file, var(1:nvars_file))
    ELSE IF ( tpl_has_var_name(file) ) THEN
! If we are using a non-dump file with a variable name template, loop through
! the variables setting one from each file
      DO i = 1,nvars_file
        CALL fill_variables_from_file(                                        &
          tpl_substitute_var(file, tpl_name(i)),                              &
          (/ var(i) /), (/ var_name(i) /)                                     &
        )
      END DO
    ELSE
! We are not using a file name template, so set all variables from the same
! file
      CALL fill_variables_from_file(                                          &
        file, var(1:nvars_file), var_name(1:nvars_file)                       &
      )
    END IF
  END IF


!*****************************************************************************
! Further processing depending on options specified
!*****************************************************************************
! Calculate soil moisture content from wetness
! NOTE that at this point, sthu contains total soil wetness (sthu + sthf)
  DO i = 1,sm_levels
    smcl(:,i) = rho_water * dzsoil(i) * sthu(:,i) * smvcst(:,i)
  ENDDO


!------------------------------------------------------------------------------
!   If using two-tile urban schemes then split urban fraction here between
!   the canyon and roof. This has to be done instead of in init_urban to be
!   consistent with triffid
!------------------------------------------------------------------------------
! The ice tile is no longer hijacked for standalone JULES, like it is currently
! in the UM, so the check for the presence of ice does not need to be done
! here. However, it has been left it in to be consistent.
  IF ( l_urban2t ) THEN
    CALL log_info("init_ic",                                                  &
                  "Either URBAN-2T or MORUSES is in use - splitting urban " //&
                  "tile into canyon/roof")

    DO l = 1, land_pts
      IF (frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
        urban_fraction = frac(l,urban)
        frac(l,urban_canyon) = urban_fraction * wrr(l)
        frac(l,urban_roof)   = urban_fraction - frac(l,urban_canyon)
      ELSE IF ( frac(l,urban) > 0.0 .AND. frac(l,ice) > 0.0 ) THEN
        CALL log_warn("init_ic",                                              &
                      "WARNING: ice and urban co-exist for land point " //    &
                      TRIM(to_string(l)))
      END IF
    END DO
  END IF

!-------------------------------------------------------------------------------
!   If using TRIFFID (with or without competing veg), ensure that fractions of
!   PFTs are not below minimum.
!   Only do this over soil points - land ice points should have zero fractions.
!-------------------------------------------------------------------------------
  IF ( l_triffid ) THEN
! Set up a flag to see if any points were reset to frac_min
    reset_done = .FALSE.
    DO j = 1,soil_pts
      i = soil_index(j)
      IF ( ANY( frac(i,:) < frac_min ) ) THEN
! Reset all small values. Renormalisation is done later, but will fail if
! frac_min is sufficiently large. We only reset PFT tiles
        WHERE ( frac(i,1:npft) < frac_min ) frac(i,1:npft) = frac_min
        reset_done = .TRUE.
      ENDIF
    ENDDO

    IF ( reset_done )                                                         &
      CALL log_warn("init_ic",                                                &
                    "frac < frac_min at one or more points - reset to " //    &
                    "frac_min at those points")
  ENDIF

!-------------------------------------------------------------------------------
! Check that frac sums to 1.0 (with a bit of leeway).
!-------------------------------------------------------------------------------
  DO i = 1,land_pts
    IF ( ABS( SUM(frac(i,:))-1.0 ) >= 1.0e-2 ) THEN
! If the discrepancy is big enough, bail
      CALL log_fatal("init_ic",                                               &
                     "frac does not sum to 1 at point " //                    &
                     TRIM(to_string(i)) //                                    &
                     " and the discrepancy is too big to be removed")
    ELSE IF ( ABS( SUM(frac(i,:)) - 1.0 ) > 1.0e-4 ) THEN
! Correct a small discrepancy
      CALL log_warn("init_ic",                                                &
                    "frac does not sum to 1 at point " //                     &
                     TRIM(to_string(i)) //                                    &
                     " - removing small discrepancy")
! Ignore small discrepancies and (re)normalise
      frac(i,:) = frac(i,:) / SUM(frac(i,:))
    ENDIF

! Ignore any discrepancy below the threshold completely
  ENDDO

!-------------------------------------------------------------------------------
! Process the ice fraction field.
! Identify land ice points, only if the ice surface type is specified.
!-------------------------------------------------------------------------------
  lice_pts = 0
  lice_index(:) = 0
  IF ( ice > 0 ) THEN
    DO l = 1,land_pts
      IF ( frac(l,ice) > 0.0 ) THEN
! This is a land ice point.
        lice_pts = lice_pts + 1
        lice_index(lice_pts) = l
! At present, land ice and soil points are mutually exclusive.
! Check this is not a soil point
        IF ( ANY(soil_index == l) )                                           &
          CALL log_fatal("init_ic",                                           &
                         "Land ice points and soil points are mutually exclusive")

! Check that ice fraction is one (cannot have partial ice coverage).
        IF ( ABS(frac(l,ice) - 1.0) > EPSILON(1.0) )                          &
          CALL log_error("init_ic",                                           &
                         "Ice fraction must be 1 at an ice point")
      ENDIF
    ENDDO

    CALL log_info("init_ic",                                                  &
                  "Number of land ice points = " // TRIM(to_string(lice_pts)))
  ENDIF

!-----------------------------------------------------------------------------
! Check that all land points have been identified as either soil or ice
!-----------------------------------------------------------------------------
  IF ( soil_pts + lice_pts /= land_pts )                                      &
    CALL log_fatal("init_ic",                                                 &
                   "All points should be either soil or land ice points - " //&
                   "have land_pts = " // TRIM(to_string(land_pts)) //         &
                   "and soil_pts + lice_pts = " //                            &
                   TRIM(to_string(soil_pts + lice_pts)))

!-----------------------------------------------------------------------------
! Set up tile index
!-----------------------------------------------------------------------------
  CALL TILEPTS( LAND_PTS, FRAC, TILE_PTS, TILE_INDEX )

! For URBAN-2T or MORUSES: Check that urban canyons also have roofs
  IF ( l_urban2t ) THEN
    IF ( tile_pts(urban_canyon) /= tile_pts(urban_roof) )                     &
      CALL log_fatal("init_ic",                                               &
                     "URBAN-2T or MORUSES - # canyons /= # roofs")
  END IF

!-----------------------------------------------------------------------------
! Deal with "simple" initialisation of snow variables.
!-----------------------------------------------------------------------------
  IF ( total_snow ) CALL total_snow_init()

!-----------------------------------------------------------------------------
! Calculate frozen and unfrozen fractions of soil moisture.
!-----------------------------------------------------------------------------
  CALL FREEZE_SOIL (LAND_PTS, SM_LEVELS,                                      &
                           B, DZSOIL, SATHH, SMCL, T_SOIL, SMVCST, STHU, STHF)

!-----------------------------------------------------------------------------
! Finish initialising TOPMODEL
!-----------------------------------------------------------------------------
  IF ( l_top ) CALL topmodel_init()

  RETURN

END SUBROUTINE init_ic
#endif
