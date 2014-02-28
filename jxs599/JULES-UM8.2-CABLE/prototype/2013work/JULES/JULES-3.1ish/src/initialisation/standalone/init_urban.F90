#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/init_urban.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

SUBROUTINE init_urban()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE input_mod, ONLY : fill_variables_from_file

  USE switches, ONLY : l_cosz

  USE ancil_info, ONLY : land_pts, tile_pts, tile_index, frac

  USE nstypes, ONLY : npft, ice, urban, urban_canyon, urban_roof

  USE switches_urban

  USE urban_param

  USE nvegparm

  USE c_z0h_z0m, ONLY : z0h_z0m

  USE logging_mod, ONLY : log_info, log_debug, log_warn, log_error, log_fatal

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises urban parameters and properties
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_URBAN_VARS = 9
         ! The maximum possible number of TOPMODEL variables that can be given

  INTEGER :: nvars_required      ! The number of variables that are
                                 ! required in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_URBAN_VARS)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)

  REAL :: sc_hwr(land_pts), d_h(land_pts)  ! Work variables
  REAL :: lambdaf, lambdap

  INTEGER :: i,l  ! Index variables

  INTEGER :: error, error_sum  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the urban_properties namelist - this specifies how urban
! properties are set
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_URBAN_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_URBAN_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_URBAN_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_URBAN_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(MAX_URBAN_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /urban_properties/ file, nvars, var, use_file, var_name, tpl_name, &
                              const_val


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nvars_required = 0
  nvars_file     = 0
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is to read every variable from file

!-----------------------------------------------------------------------------
! Do checks on urban tiles to see if an urban-2T scheme is enabled
!-----------------------------------------------------------------------------
  IF ( urban > 0 .AND. urban_canyon > 0 )                                     &
    CALL log_fatal("init_urban",                                              &
                   "Cannot have both urban and urban_canyon tiles")

! The urban_roof has to be specified to use the two-tile urban schemes. Check
! that the roof has the correct counterparts
  IF ( urban_roof > 0 ) THEN
    IF ( urban < 0 .AND. urban_canyon < 0 )                                   &
      CALL log_fatal("init_urban",                                            &
                     "Cannot have urban_roof without urban_canyon - to " //   &
                     "use URBAN-1T scheme set urban_roof tile to -1")
  ELSE
! Check that if not present then neither is the canyon or MORUSES
    IF ( urban_canyon > 0 )                                                   &
      CALL log_fatal("init_urban",                                            &
                     "Cannot have urban_canyon without urban_roof - to " //   &
                     "use URBAN-1T scheme set urban_canyon tile to -1")
    IF ( l_moruses )                                                          &
      CALL log_fatal("init_urban",                                            &
                     "MORUSES has no urban roof surface type")
  END IF

! Either "urban_canyon" or "urban" can be specified as
! urban fraction as a whole is specified and then is split between canyon &
! roof depending on canyon fraction in init_urban.
  IF ( urban_roof > 0 ) THEN
    l_urban2t = .TRUE.
    IF ( urban_canyon < 0 ) urban_canyon = urban
    IF ( l_moruses ) THEN
      CALL log_info("init_urban",                                             &
                    "Using urban surface scheme - MORUSES")
    ELSE
      CALL log_info("init_urban",                                             &
                    "Using urban surface scheme - URBAN-2T")
    END IF
  END IF

! If the two-tile urban schemes are not used, make sure the appropriate
! switches are set to .FALSE. to avoid unnecessary calculations and leave.
! The urban arrays are not allocated yet so there is no need to initialise
! them to zero
  IF ( .NOT. l_urban2t ) THEN
    l_urban_empirical      = .FALSE.
    l_moruses_albedo       = .FALSE.
    l_moruses_emissivity   = .FALSE.
    l_moruses_rough        = .FALSE.
    l_moruses_storage      = .FALSE.
    l_moruses_storage_thin = .FALSE.
    l_moruses_macdonald    = .FALSE.
    RETURN
  END IF

!-----------------------------------------------------------------------------
! Read the namelists, now we know they are required
!-----------------------------------------------------------------------------
! Open the urban namelist file
  OPEN(NAMELIST_UNIT, FILE='urban.nml',                                       &
                 STATUS='old', POSITION='rewind', ACTION='read', IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_urban",                                              &
                   "Error opening namelist file urban.nml " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! There are three namelists to read from this file
  CALL log_info("init_urban", "Reading URBAN_SWITCHES namelist...")
  READ(NAMELIST_UNIT, nml=urban_switches, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_urban",                                              &
                   "Error reading namelist URBAN_SWITCHES " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_urban", "Reading URBAN2T_PARAM namelist...")
  READ(NAMELIST_UNIT, nml=urban2t_param, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_urban",                                              &
                   "Error reading namelist URBAN2T_PARAM " //                 &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

  CALL log_info("init_urban", "Reading URBAN_PROPERTIES namelist...")
  READ(NAMELIST_UNIT, nml=urban_properties, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_urban",                                              &
                   "Error reading namelist URBAN_PROPERTIES " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


! Close the namelist file
  CLOSE(NAMELIST_UNIT, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_urban",                                              &
                   "Error closing namelist file urban.nml " //                &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")


!-----------------------------------------------------------------------------
! Process the namelist values and set derived variables
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! Check that the run actually has urban and that urban schemes are not run
! in error
!-----------------------------------------------------------------------------
  CALL tilepts(land_pts, frac, tile_pts, tile_index)
! urban_canyon tile number is set above so should be present if two-tiles
! schemes are used
  IF ( tile_pts(urban_canyon) == 0 )                                          &
    CALL log_warn("init_urban",                                               &
                  "URBAN-2T or MORUSES is selected but there are no " //      &
                  "urban land points - extra calculations may be being " //   &
                  "performed that will not impact on results")

!-----------------------------------------------------------------------------
! Check urban switch logic is sensible
!-----------------------------------------------------------------------------
  IF ( l_moruses ) THEN
! Check that MORUSES has some parametrisations turned on
    IF ( .NOT. l_moruses_albedo .AND. .NOT. l_moruses_emissivity &
       .AND. .NOT. l_moruses_rough .AND. .NOT. l_moruses_storage )            &
      CALL log_fatal("init_urban", "MORUSES has no parametrisations turned on")
  ELSE
! If MORUSES is not used then ALL MORUSES switches MUST be .FALSE.
    l_moruses_albedo       = .FALSE.
    l_moruses_emissivity   = .FALSE.
    l_moruses_rough        = .FALSE.
    l_moruses_storage      = .FALSE.
    l_moruses_storage_thin = .FALSE.
    l_moruses_macdonald    = .FALSE.
  END IF

! Check MORUSES switch logic
  IF ( l_moruses .AND. l_urban_empirical .AND. .NOT. l_moruses_macdonald ) THEN
    CALL log_info("init_urban",                                               &
                  "Using l_urban_empirical with MORUSES - MacDonald (1998) "//&
                  "formulation for roughness length for momentum and " //     &
                  "displacement height must be used for consistency")
    l_moruses_macdonald = .TRUE.
  END IF

  IF ( l_moruses_albedo .AND. .NOT. l_cosz ) THEN
    CALL log_info("init_urban", "Using l_moruses_albedo - l_cosz must be T")
    l_cosz = .TRUE.
  END IF

  CALL log_info("init_urban",                                                 &
                "l_urban2T = " // TRIM(to_string(l_urban2T)) // "; " //       &
                "l_urban_empirical = " // TRIM(to_string(l_urban_empirical)) // &
                "; l_moruses = " // TRIM(to_string(l_moruses)))


!-----------------------------------------------------------------------------
! Copy urban-specific values of nveg properties into convenience variables
!-----------------------------------------------------------------------------
  albsnc_c   = albsnc_nvg(urban_canyon-npft)
  albsnc_rf  = albsnc_nvg(urban_roof-npft)

  albsnf_c   = albsnf_nvg(urban_canyon-npft)
  albsnf_rf  = albsnf_nvg(urban_roof-npft)

  catch_c    = catch_nvg(urban_canyon-npft)
  catch_rf   = catch_nvg(urban_roof-npft)

  gs_c       = gs_nvg(urban_canyon-npft)
  gs_rf      = gs_nvg(urban_roof-npft)

  infil_c    = infil_nvg(urban_canyon-npft)
  infil_rf   = infil_nvg(urban_roof-npft)

  z0_c       = z0_nvg(urban_canyon-npft)
  z0_rf      = z0_nvg(urban_roof-npft)

  z0h_z0m_c  = z0h_z0m(urban_canyon)
  z0h_z0m_rf = z0h_z0m(urban_roof)

  ch_c       = ch_nvg(urban_canyon-npft)
  ch_rf      = ch_nvg(urban_roof-npft)

  vf_c       = vf_nvg(urban_canyon-npft)
  vf_rf      = vf_nvg(urban_roof-npft)

  emis_c     = emis_nvg(urban_canyon-npft)
  emis_rf    = emis_nvg(urban_roof-npft)


!-----------------------------------------------------------------------------
! Allocate and initialise arrays
!-----------------------------------------------------------------------------
  error_sum = 0
  ALLOCATE( wrr(land_pts), STAT=error )
  error_sum = error_sum + error
  IF ( error_sum == 0 ) wrr(:) = 0.0

  IF ( l_moruses ) THEN
    ALLOCATE( hgt(land_pts)  , STAT=error )
    error_sum = error_sum + error
    ALLOCATE( hwr(land_pts)  , STAT=error )
    error_sum = error_sum + error
    ALLOCATE( disp(land_pts) , STAT=error )
    error_sum = error_sum + error
    ALLOCATE( ztm(land_pts)  , STAT=error )
    error_sum = error_sum + error
    ALLOCATE( albwl(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( albrd(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( emisw(land_pts), STAT=error )
    error_sum = error_sum + error
    ALLOCATE( emisr(land_pts), STAT=error )
    error_sum = error_sum + error

    IF  ( error_sum == 0 ) THEN
      hgt(:)   = 0.0
      hwr(:)   = 0.0
      albwl(:) = 0.0
      albrd(:) = 0.0
      emisw(:) = 0.0
      emisr(:) = 0.0
      ztm(:)   = 0.0
      disp(:)  = 0.0
    ENDIF
  ENDIF

! Check for error.
  IF ( error_sum /= 0 )                                                       &
    CALL log_fatal("init_urban", "Error allocating urban arrays")


!-----------------------------------------------------------------------------
! Read values for urban properties
!-----------------------------------------------------------------------------
! Set up the required variables
  IF ( l_moruses ) THEN
! All variables are required for MORUSES
    nvars_required = MAX_URBAN_VARS
    required_vars(:) = (/ 'wrr  ', 'hwr  ', 'hgt  ', 'ztm  ', 'disp ',        &
                          'albwl', 'albrd', 'emisw', 'emisr' /)
  ELSE
! For urban2t, we only need wrr
    nvars_required = 1
    required_vars(1) = 'wrr'
  END IF

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_urban",                                            &
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
        CALL log_info("init_urban",                                           &
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
        CALL log_info("init_urban",                                           &
                      "'" // TRIM(var(i)) // "' will be set to a constant")

        CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
      END IF
    ELSE
! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_urban",                                             &
                    "Provided variable '" // TRIM(var(i)) //                  &
                    "' is not required, so will be ignored")
    END IF
  END DO

!-----------------------------------------------------------------------------
! Set variables from file
!-----------------------------------------------------------------------------
  IF ( nvars_file > 0 ) THEN
    IF ( tpl_has_var_name(file) ) THEN
! We are using a file name template, so loop through the variables setting
! one from each file
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


!-----------------------------------------------------------------------------
! Empirical relationships derived from correlating CEH urban fraction and
! LUCID urban geometry data for London. Obtained from collaboration with the
! University of Reading. See:
!     Bohnenstengel, S.I., Evans, S., Clark, P., Belcher, S.E. (2010);
!     Simulations of the London urban heat island, Q.J.R.Meteorol. Soc., to
!     be submitted.
! for more information

! Check for ice has been left in to be consistent with UM, but is not actually
! required here
!-----------------------------------------------------------------------------
  IF ( l_urban_empirical ) THEN
    CALL log_info("init_urban",                                               &
                  "Using empirical relationships for urban geometry: wrr")

    IF ( l_moruses )                                                          &
      CALL log_info("init_urban",                                             &
                    "Using empirical relationships for urban geometry: hwr")

    DO l = 1, land_pts
      IF ( frac(l,urban) > 0.0 .AND. ABS(frac(l,ice)) < EPSILON(1.0) ) THEN
        lambdap = 22.878*frac(l,urban)**6 - 59.473*frac(l,urban)**5           &
                + 57.749*frac(l,urban)**4 - 25.108*frac(l,urban)**3           &
                + 4.3337*frac(l,urban)**2 + 0.1926*frac(l,urban)              &
                + 0.036
        wrr(l) = 1.0 - lambdap

        IF ( l_moruses ) THEN
          lambdaf = 16.412*frac(l,urban)**6 - 41.855*frac(l,urban)**5         &
                  + 40.387*frac(l,urban)**4 - 17.759*frac(l,urban)**3         &
                  + 3.2399*frac(l,urban)**2 + 0.0626*frac(l,urban)            &
                  + 0.0271
          hwr(l) = 4.0 * ATAN(1.0)/2.0 * lambdaf / ( 1.0 - lambdap )
        END IF
      END IF
    END DO
  END IF

  IF ( l_moruses ) THEN
    IF ( l_urban_empirical ) THEN
      CALL log_info("init_urban",                                             &
                    "Using empirical relationships for urban geometry: hgt")

      DO l = 1,land_pts
        IF (frac(l,urban) > 0.0 .AND. ABS(frac(l,ice)) < EPSILON(1.0) ) THEN
          hgt(l) = 167.409  * frac(l,urban)**5 - 337.853  * frac(l,urban)**4  &
                 + 247.813  * frac(l,urban)**3 -  76.3678 * frac(l,urban)**2  &
                 +  11.4832 * frac(l,urban)    +   4.48226
        END IF
      END DO
    END IF

    IF ( l_moruses_macdonald ) THEN
! Macdonald Formulation
      CALL log_info("init_urban", "Using MacDonald formulation")

      sc_hwr(:) = 0.5 * ( hwr(:) / (2.0 * ATAN(1.0)) )
      d_h(:)    = 1.0 - wrr(:) * ( a**(wrr(:) - 1.0) )
      disp(:)   = d_h(:) * hgt(:)

      DO l = 1,land_pts
        IF ( wrr(l) > 0.0 ) THEN
          ztm(l) = (cdz * (1.0 - d_h(l)) * sc_hwr(l) * wrr(l) / kappa2)**(-0.5)
          ztm(l) = (1.0 - d_h(l)) * EXP(-ztm(l))
          ztm(l) = ztm(l) * hgt(l)
          ztm(l) = MAX(ztm(l), z0m_mat)
        END IF
      END DO
    END IF
  END IF

  RETURN

END SUBROUTINE init_urban
#endif
