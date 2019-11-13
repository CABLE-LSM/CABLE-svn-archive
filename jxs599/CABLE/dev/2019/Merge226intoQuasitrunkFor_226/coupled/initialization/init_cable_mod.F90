#if !defined(UM_JULES) 
!jh:#this will need to be removed 

MODULE init_cable_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_cable_grid()

USE cable_types_mod,          ONLY : mp
USE cable_types_mod,          ONLY : l_tile_pts
USE ancil_info,               ONLY : surft_pts, frac_surft, land_pts
USE jules_surface_types_mod,  ONLY : ntype

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

INTEGER :: i, j

! Determine active tiles map
l_tile_pts(:,:) = .FALSE.

DO j = 1, ntype
  DO i = 1, land_pts
    IF ( frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO

RETURN

END SUBROUTINE init_cable_grid

SUBROUTINE init_cable_veg()
USE cable_types_mod,            ONLY : mp
USE cable_other_constants_mod,  ONLY : nrb
USE cable_types_mod,            ONLY : SurfaceTypeID => SurfaceTypeID_cbl
USE cable_types_mod,            ONLY : l_tile_pts
USE cable_types_mod,            ONLY : VegIN 
USE cable_types_mod,            ONLY : VegXfang, VegTaul, VegRefl
USE ancil_info,                 ONLY : nsurft, land_pts, frac_surft
USE init_params_mod,            ONLY : init_veg_from_vegin_JAC
IMPLICIT NONE
integer :: JSurfaceTypeID(land_pts,nsurft)  
integer :: i
!-----------------------------------------------------------------------------
! Description:
!   init_cables veg parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

!pack surface type: L_tile_pts mask set FROM surf_couple_rad
SurfaceTypeID = 0
JSurfaceTypeID = 0
do i=1,nsurft
  if( frac_surft(1,i) > 0 ) JSurfaceTypeID(:,i) = i
end do

SurfaceTypeID = PACK( JSurfaceTypeID, L_tile_pts)

call init_veg_from_vegin_JAC(VegXfang, VegTaul, VegRefl, mp, nrb,   &
                  SurfaceTypeID, Vegin%Xfang, Vegin%Taul, Vegin%Refl )

END SUBROUTINE init_cable_veg



SUBROUTINE init_cable_progs()

USE io_constants, ONLY : MAX_SDF_NAME_LEN, MAX_FILE_NAME_LEN, NAMELIST_UNIT

USE string_utils_mod, ONLY : to_string

USE templating_mod, ONLY : tpl_has_var_name, tpl_substitute_var

USE model_interface_mod, ONLY : IDENTIFIER_LEN, populate_var, get_var_id

USE ancil_info, ONLY : land_pts!, soil_pts, soil_index, sm_levels

USE input_mod, ONLY: fill_variables_from_file

USE logging_mod, ONLY: log_info, log_debug, log_warn, log_error, log_fatal

!  USE soil_param, ONLY : jules_soil_param, dzsoil, dzsoil_io, zsmc, zst
!  USE p_s_parms, ONLY : sathh, satcon, hcon, smvcst

!USE cable_data_mod

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads in CABLE prognostic variables and checks them for consistency
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------
! Work variables
  INTEGER, PARAMETER :: nCABLE_VARS  = 10

  INTEGER :: nvars_required      ! The number of soil variables that are
                                 ! required in this configuration
  
  CHARACTER(len=IDENTIFIER_LEN) :: required_vars(nCABLE_VARS)
                                 ! The variable identifiers of the required
                                 ! variables

  INTEGER :: nvars_file       ! The number of variables that will be set
                              ! from the given file (template?)

! Variables passed to fill_variables_from_file
  CHARACTER(len=IDENTIFIER_LEN) :: file_var(nCABLE_VARS)
                        ! The variable identifiers of the variables to set
                        ! from file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_var_name(nCABLE_VARS)
                        ! The name of each variable in the file

  CHARACTER(len=MAX_SDF_NAME_LEN) :: file_tpl_name(nCABLE_VARS)
                        ! The name to substitute in a template for each
                        ! variable

  INTEGER :: i,l  ! Index variables

  INTEGER :: error  ! Error indicator

!-----------------------------------------------------------------------------
! Definition of the cable_progs namelist
!-----------------------------------------------------------------------------
  CHARACTER(len=MAX_FILE_NAME_LEN) :: file
                        ! The name of the file (or variable name template) to
                        ! use for variables that need to be filled from file

  INTEGER :: nvars      ! The number of variables in this section

  CHARACTER(len=IDENTIFIER_LEN) :: var(nCABLE_VARS)
                        ! The variable identifiers of the variables
  LOGICAL :: use_file(nCABLE_VARS)
                        !   T - the variable uses the file
                        !   F - the variable is set using a constant value
  CHARACTER(len=MAX_SDF_NAME_LEN) :: var_name(nCABLE_VARS)
                        ! The name of each variable in the file
  CHARACTER(len=MAX_SDF_NAME_LEN) :: tpl_name(nCABLE_VARS)
                        ! The name to substitute in a template for each
                        ! variable
  REAL :: const_val(nCABLE_VARS)
  INTEGER:: iconst_val(nCABLE_VARS)
                        ! The constant value to use for each variable if
                        ! use_file = F for that variable

  NAMELIST /cable_progs/ file, nvars, use_file, var, var_name

!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Initialise
!-----------------------------------------------------------------------------
  nvars_required = 0
  nvars_file     = 0
  nvars          = 0
  use_file(:)    = .TRUE.  ! Default is for every variable to be read from file

!-----------------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------------
  CALL log_info("init_cable", "Reading CABLE_VARS namelist...")

! First, we read the soil_param namelist
  READ(NAMELIST_UNIT, nml=cable_progs, IOSTAT=error)

 IF ( error /= 0 )                                                           &
    CALL log_fatal("init_cable",                                               &
                   "Error reading namelist CABLE" //              &
                   "(IOSTAT=" // TRIM(to_string(error)) // ")")

!-----------------------------------------------------------------------------
! Set up soil properties using namelist values
!-----------------------------------------------------------------------------
! Set up the required variables
! All the CABLE variables are always required for CABLE runs
  nvars_required = nCABLE_VARS
  required_vars(:) = (/ &
                          'ThreeLayerSnowFlag_CABLE',                          &
                          'OneLyrSnowDensity_CABLE ',                          &
                          'SnowAge_CABLE           ',                          &
                          'SnowDensity_CABLE       ',                          &
                          'SnowMass_CABLE          ',                          &
                          'SnowDepth_CABLE         ',                          &
                          'SnowTemp_CABLE          ',                          &
                          'FrozenSoilFrac_CABLE    ',                          &
                          'SoilMoisture_CABLE      ',                          &
                          'SoilTemp_CABLE          '                           &
                          /)

!-----------------------------------------------------------------------------
! Check that all the required variables are there
!-----------------------------------------------------------------------------

  DO i = 1,nvars_required
    IF ( .NOT. ANY(var(1:nvars) == TRIM(required_vars(i))) )                   &
    !IF ( trim( var(1) ) /= TRIM( required_vars(1) ) )                         &
      CALL log_fatal("init_cable",                                             &
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
        CALL log_info("init_cable",                                            &
                      "'" // TRIM(var(i)) // "' will be read from file")

! If the variable will be filled from file, register it here
        nvars_file = nvars_file + 1
        file_var(nvars_file) = var(i)
        file_var_name(nvars_file) = var_name(i)
        file_tpl_name(nvars_file) = tpl_name(i)
      ELSE
! If the variable is being set as a constant, just populate it here
        CALL log_info("init_cable",                                            &
                      "'" // TRIM(var(i)) // "' will be set to a constant")

        CALL populate_var(get_var_id(var(i)), CONST_VAL=const_val(i))
      END IF
    ELSE
! If the variable is not a required variable, warn about not using it
      CALL log_warn("init_cable",                                              &
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

  RETURN

END SUBROUTINE init_cable_progs

END MODULE init_cable_mod

#endif
