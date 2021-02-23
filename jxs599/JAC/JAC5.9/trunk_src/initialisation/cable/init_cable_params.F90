MODULE init_cable_params_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: init_cable_params

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='init_cable_params_mod'

CONTAINS

SUBROUTINE init_cable_params(mp, frac_surft )
USE allocate_veg_params_mod,      ONLY: allocate_veg_parameter_type
USE allocate_soil_params_mod,      ONLY: allocate_soil_parameter_type
USE cable_params_mod,         ONLY: soil => soil_cbl
USE cable_params_mod,         ONLY: veg => veg_cbl
USE init_cable_pftparms_mod,  ONLY: init_cable_veg
USE init_cable_soilparms_mod, ONLY: init_cable_soil
USE ancil_info,        ONLY: nsurft, land_pts

IMPLICIT NONE

INTEGER :: mp
REAL:: frac_surft(land_pts,nsurft)  
!-----------------------------------------------------------------------------
! Description:
! Set CABLE params as are used in CABLE (per active tile) - 
! - from namelist values read (per PFT/soil type)
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE TECHNICAL
!-----------------------------------------------------------------------------

CALL allocate_soil_parameter_type(soil, mp)
CALL allocate_veg_parameter_type(veg, mp)
CALL init_cable_veg( frac_surft)
CALL init_cable_soil()

END SUBROUTINE init_cable_params

END MODULE init_cable_params_mod
