#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************


SUBROUTINE get_var_attrs(var_id, int_attributes, real_attributes, char_attributes)

  USE io_constants, ONLY : MAX_ATTR_VAR, MDI

  USE dictionary_mod, ONLY : dict, dict_create, dict_set

  USE string_utils_mod, ONLY : to_string

  USE switches, ONLY : l_triffid

  USE soil_param, ONLY : zsmc

  USE model_grid_mod, ONLY : latitude, longitude

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given an identifer for a model variable, returns the attributes that
!   should be defined on output files
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  INTEGER, INTENT(IN) :: var_id  ! Identifies the model variable
  TYPE(dict), INTENT(OUT) :: int_attributes  ! Attributes with integer values
  TYPE(dict), INTENT(OUT) :: real_attributes  ! Attributes with real values
  TYPE(dict), INTENT(OUT) :: char_attributes  ! Attributes with character values


! Work variables
  CHARACTER(len=IDENTIFIER_LEN) :: identifier  ! The string identifier for the
                                               ! variable


!-----------------------------------------------------------------------------


! Create the dictionary objects
  int_attributes = dict_create(MAX_ATTR_VAR, INT(1))
  real_attributes = dict_create(MAX_ATTR_VAR, REAL(1.0))
  char_attributes = dict_create(MAX_ATTR_VAR, "char")


! Get the string identifier
  identifier = get_string_identifier(var_id)


! All variables except lat and lon will have a missing data attribute and
! a coordinate attribute
! lat and lon do not need these since they are interpreted as coordinate variables
! as per the CF convention
  IF ( identifier /= 'latitude' .AND. identifier /= 'longitude' ) THEN
    CALL dict_set(real_attributes, "missing_value", MDI)
    CALL dict_set(real_attributes, "_FillValue", MDI)
    CALL dict_set(char_attributes, "coordinates", "latitude longitude")
  END IF


! Add the standard "long_name" and "units" attributes from the metadata array
! var_id is an index into the metadata array
  CALL dict_set(char_attributes, "long_name", metadata(var_id)%long_name)
  CALL dict_set(char_attributes, "units", metadata(var_id)%units)


! Use a select statement to add extra attributes/override current values if required
! The SELECT statement uses the string identifier to select a CASE to avoid being
! dependent on the implementation of integer ids
  SELECT CASE ( identifier )
    CASE ( 'latitude' )
      CALL dict_set(char_attributes, "standard_name", "latitude")

    CASE ( 'longitude' )
      CALL dict_set(char_attributes, "standard_name", "longitude")

    CASE ( 'frac' )
      CALL dict_set(real_attributes, "valid_min", 0.0)
      CALL dict_set(real_attributes, "valid_max", 1.0)

    CASE ( 'smc_avail_top' )
! Overwrite long_name to use a value that is calculated using zsmc
      CALL dict_set(char_attributes, "long_name",                             &
                    "Gridbox available moisture in top " // TRIM(to_string(zsmc)) // "m of soil")

    CASE ( 'cs' )
! Overwrite long_name to use a value that specifies the soil carbon pools
      IF ( l_triffid ) THEN
        CALL dict_set(char_attributes, "long_name", "Gridbox soil carbon in each pool (DPM,RPM,bio,hum)")
      ELSE
        CALL dict_set(char_attributes, "long_name", "Gridbox soil carbon (single pool)")
      END IF

    CASE ( 'resp_s' )
! Overwrite long_name to use a value that specifies the soil carbon pools
      IF ( l_triffid ) THEN
        CALL dict_set(char_attributes, "long_name", "Gridbox soil respiration from each pool (DPM,RPM,bio,hum)")
      ELSE
        CALL dict_set(char_attributes, "long_name", "Gridbox soil respiration (single pool)")
      END IF
  END SELECT

  RETURN

END SUBROUTINE get_var_attrs
#endif
