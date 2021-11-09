
MODULE cable_fields_mod

! module containing instances of the data types for CABLE in jules standalone
! analagous to jules_fields_mod. This is the central place where ubiquitous
! fields are distributed and USEd from throughout the model are instantiated.
! Currently only prognostic fields are here but many more are on the way

!type definitions
USE progs_cbl_vars_mod,       ONLY: progs_cbl_vars_type,                      &
                                    progs_cbl_vars_data_type
USE veg_params_type_mod_cbl,  ONLY: veg_params_data_type_cbl ,                &
                                    veg_params_type_cbl
USE soil_params_type_mod_cbl, ONLY: soil_params_data_type_cbl,                &
                                    soil_params_type_cbl

PUBLIC

!TYPES to hold the data
TYPE(progs_cbl_vars_data_type),   TARGET :: progs_cbl_vars_data
TYPE(veg_params_data_type_cbl),   TARGET :: veg_pars_data_cbl
TYPE(soil_params_data_type_cbl),  TARGET :: soil_pars_data_cbl

!TYPES we pass around. These happen to be pointers to the data types above
!but this should be transparent
TYPE(progs_cbl_vars_type) :: progs_cbl_vars
TYPE(veg_params_type_cbl)  :: veg_pars_cbl
TYPE(soil_params_type_cbl) :: soil_pars_cbl

END MODULE cable_fields_mod
