#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/initialisation/standalone/ancillaries/init_soil.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE init_soil()

  USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

  USE string_utils_mod, ONLY : to_string

  USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

  USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

  USE ancil_info, ONLY : land_pts, soil_pts, soil_index, sm_levels

  USE soil_param, ONLY : jules_soil_param, dzsoil, dzsoil_io, zsmc, zst

  USE p_s_parms, ONLY : sathh, satcon, hcon, smvcst

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in the soil properties and checks them for consistency
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: MAX_SOIL_VARS = 9  ! The maximum possible number of
                                           ! soil variables that can be given

  INTEGER :: nvars_required      ! The number of soil variables that are
                                 ! required in this configuration
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(MAX_SOIL_VARS)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)
! Variables passed to fill_variables_from_file
  CHARACTER(len=IDENTIFIER_LEN) :: file_var(MAX_SOIL_VARS)
                        ! The variable identifiers of the variables to set
                        ! from file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_var_name(MAX_SOIL_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_tpl_name(MAX_SOIL_VARS)
                        ! The name to substitute in a template for each
                        ! variable

  INTEGER :: i,l  ! Index variables

  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the jules_soil_props namelist
!-----------------------------------------------------------------------------
  LOGICAL :: const_z            ! T - the same properties are used for each
                                !     soil layer
                                ! F - properties for each layer are read from
                                !     file
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section
  CHARACTER(len=IDENTIFIER_LEN) :: var(MAX_SOIL_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(MAX_SOIL_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(MAX_SOIL_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(MAX_SOIL_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(MAX_SOIL_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable
  NAMELIST /jules_soil_props/ const_z, file, nvars, var, use_file, var_name,  &
                                             tpl_name, const_val


!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nvars_required = 0
  nvars_file     = 0
  const_z        = .FALSE. ! Default is to read a value for each soil level
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is for every variable to be read from file

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_soil", "Reading JULES_SOIL_PARAM namelist...")

! First, we read the soil_param namelist
  READ(NAMELIST_UNIT, nml=jules_soil_param, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error reading namelist JULES_SOIL_PARAM " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

! Initialise dzsoil by copying from fixed length IO counterpart
  dzsoil(:) = dzsoil_io(1:sm_levels)

! Check that zsmc and zst are within soil depth.
  IF ( SUM(dzsoil(:)) < zsmc )                                                &
    CALL log_fatal("init_soil",                                               &
                   "zsmc is below bottom of soil column")

  IF ( SUM(dzsoil(:)) < zst )                                                 &
    CALL log_fatal("init_soil",                                               &
                   "zst is below bottom of soil column")


  CALL log_info("init_soil", "Reading JULES_SOIL_PROPS namelist...")

! Then read the soil properties namelist
  READ(NAMELIST_UNIT, nml=jules_soil_props, IOSTAT=error)
  IF ( error /= 0 )                                                           &
    CALL log_fatal("init_soil",                                               &
                   "Error reading namelist JULES_SOIL_PROPS " //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up soil properties using namelist values
!-----------------------------------------------------------------------------
! Set up the required variables
! All the soil variables are always required
  nvars_required = MAX_SOIL_VARS
  required_vars(:) = (/ 'b      ', 'sathh  ', 'satcon ', 'sm_sat ', 'sm_crit',&
                        'sm_wilt', 'hcap   ', 'hcon   ', 'albsoil' /)


!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------
  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == required_vars(i)) )                        &
      CALL log_fatal("init_soil",                                             &
                     "No value given for required variable '" //              &
                     TRIM(required_vars(i)) // "'")
  END DO

!-----------------------------------------------------------------------------
! Constant Z (i.e. spacially varying but constant through vertical levels) is
! implemented by having a separate input variable in model_interface_mod
! called <var>_const_z that has no vertical levels
!
! Hence, once the previous check is done, we add _const_z to both required
! and provided variable identifiers if asked for
!-----------------------------------------------------------------------------
  IF ( const_z ) THEN
    DO i = 1,nvars
! albsoil is the only soil variable that does not have multiple levels
      IF ( var(i) /= 'albsoil' )                                              &
        var(i) = TRIM(var(i)) // "_const_z"
    END DO

    DO i = 1,nvars_required
! albsoil is the only soil variable that does not have multiple levels
      IF ( required_vars(i) /= 'albsoil' )                                    &
        required_vars(i) = TRIM(required_vars(i)) // "_const_z"
    END DO
  END IF

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
        CALL log_info("init_soil",                                            &
                      "'" // TRIM(var(i)) // "' will be read from file")

! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        file_var(nvars_file) = var(i)
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = tpl_name(i)
      ELSE
! If the variable is being set as a constant, just populate it here
        CALL log_info("init_soil",                                            &
                      "'" // TRIM(var(i)) // "' will be set to a constant")

        CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
      END IF
    ELSE
! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_soil",                                              &
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
          tpl_substitute_var(file, file_tpl_name(i)),                         &
          (/ file_var(i) /), (/ file_var_name(i) /)                           &
        )
      END DO
    ELSE
! We are not using a file name template, so set all variables from the same
! file
      CALL fill_variables_from_file(                                          &
        file,file_var(1:nvars_file), file_var_name(1:nvars_file)              &
      )
    END IF
  END IF

!-----------------------------------------------------------------------------
! Set up derived values
!-----------------------------------------------------------------------------
! Set surface values.
  hcon(:,0)   = hcon(:,1)
  satcon(:,0) = satcon(:,1)

! Check that sathh>=0 - a common error!
  IF ( ANY( sathh(:,:) < 0.0 ) )                                              &
    CALL log_fatal("init_soil",                                               &
                   "sathh < 0.0 detected - for JULES, sathh is " //           &
                   "abs(saturated head)")

!-------------------------------------------------------------------------------
! Detect soil points.
! If top layer saturation moisture content > 0, this is a soil point.
! Note that land ice points are no longer assigned here.
!---------------------------------smvcst----------------------------------------------
  soil_pts = 0
  soil_index(:) = 0
  DO l = 1,land_pts
    IF( smvcst(l,1) > EPSILON(smvcst(l,1)) ) THEN
      soil_pts = soil_pts + 1
      soil_index(soil_pts) = l
    ENDIF
  ENDDO

  CALL log_info("init_soil",                                                  &
                "Number of soil points = " // TRIM(to_string(soil_pts)))

  IF ( soil_pts <= 0 )                                                        &
    CALL log_warn("init_soil",                                                &
                  "There are no soil points - any land points are land ice")

  RETURN

END SUBROUTINE init_soil
#endif
