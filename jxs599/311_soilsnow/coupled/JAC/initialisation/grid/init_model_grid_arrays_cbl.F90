MODULE init_model_grid_arrays_mod_cbl

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='init_model_grid_arrays_mod_cbl'

CONTAINS

SUBROUTINE init_model_grid_arrays_cbl( mp, nsurft, land_pts, surft_pts,       &
                                       frac_surft, l_tile_pts, veg, veg_ptr,  &
                                       vegin, soil, soil_ptr, soilin )
!Subrs CALLed
USE allocate_cable_arrays_mod,  ONLY: allocate_cable_arrays
USE veg_params_type_mod_cbl,    ONLY: init_veg_cbl
USE soil_params_type_mod_cbl,   ONLY: init_soil_cbl

!Type declarations
USE vegin_pars_mod_cbl,         ONLY: vegin_type
USE soilin_pars_mod_cbl,        ONLY: soilin_type
USE veg_params_type_mod_cbl,    ONLY: veg_params_data_type_cbl
USE veg_params_type_mod_cbl,    ONLY: veg_params_type_cbl
USE soil_params_type_mod_cbl,   ONLY: soil_params_data_type_cbl
USE soil_params_type_mod_cbl,   ONLY: soil_params_type_cbl

IMPLICIT  NONE
!Redeclare input args:
INTEGER, INTENT(IN) :: mp
INTEGER, INTENT(IN) :: land_pts
INTEGER, INTENT(IN) :: nsurft
INTEGER, INTENT(IN) :: surft_pts(nsurft)  
REAL,    INTENT(IN) :: frac_surft(land_pts,nsurft) 
LOGICAL, INTENT(IN) :: L_tile_pts(land_pts, nsurft)

TYPE(vegin_type)                :: vegin
TYPE(soilin_type)               :: soilin 
TYPE(veg_params_data_type_cbl)  :: veg
TYPE(veg_params_type_cbl)       :: veg_ptr
TYPE(soil_params_data_type_cbl) :: soil
TYPE(soil_params_type_cbl)      :: soil_ptr
! End HEader

! Allocate gridded arrays to kept beteen timesteps - analagous to "state" var
! NB: Many of These arrays in current implementation are NOT required here and
! in future will be treated as local/automatic vars in CABLE 
CALL allocate_cable_arrays( mp, veg, veg_ptr, soil, soil_ptr)

! Here we initialize vars, specifically ID-ed to CALL CABLE albedo() at t=1
! in future there will be further t=1 inits migrated from the interface
CALL init_veg_cbl( mp, frac_surft, veg, vegin, L_tile_pts )
CALL init_soil_cbl( mp, soil, soilin )

return
END SUBROUTINE init_model_grid_arrays_cbl

End MODULE init_model_grid_arrays_mod_cbl
