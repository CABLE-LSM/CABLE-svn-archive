MODULE allocate_cable_arrays_mod

!Common Non-science modules
USE yomhook,                  ONLY: lhook, dr_hook
USE ereport_mod,              ONLY: ereport
USE parkind1,                 ONLY: jprb, jpim

IMPLICIT NONE

PRIVATE

PUBLIC :: allocate_cable_arrays

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOCATE_CABLE_ARRAYS_MOD'

CONTAINS

SUBROUTINE allocate_cable_arrays( mp, veg, veg_ptr, soil, soil_ptr )

USE veg_params_type_mod_cbl,  ONLY: alloc_veg_params_type_cbl,                &
                                    assoc_veg_params_type_cbl,                &
                                    veg_params_data_type_cbl,                 &
                                    veg_params_type_cbl
USE soil_params_type_mod_cbl, ONLY: alloc_soil_params_type_cbl,               &
                                    assoc_soil_params_type_cbl,               &
                                    soil_params_data_type_cbl,                &
                                    soil_params_type_cbl

implicit  none
integer :: mp
TYPE(veg_params_data_type_cbl)  :: veg
TYPE(soil_params_data_type_cbl) :: soil
TYPE(veg_params_type_cbl)  :: veg_ptr
TYPE(soil_params_type_cbl) :: soil_ptr
!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

CALL alloc_veg_params_type_cbl( mp, veg )
CALL assoc_veg_params_type_cbl( veg, veg_ptr )
CALL alloc_soil_params_type_cbl( mp, soil )
CALL assoc_soil_params_type_cbl( soil, soil_ptr )

END SUBROUTINE allocate_cable_arrays

END MODULE allocate_cable_arrays_mod
