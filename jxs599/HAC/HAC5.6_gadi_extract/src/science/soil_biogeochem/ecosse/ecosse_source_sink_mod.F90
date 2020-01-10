#if !defined(UM_JULES)
!******************************COPYRIGHT**************************************
! (c) Centre for Ecology and Hydrology.
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms
! and conditions set out therein.
!
! [Met Office Ref SC0237] 
!******************************COPYRIGHT**************************************

MODULE ecosse_source_sink_mod

!-----------------------------------------------------------------------------
! Description:
!   Deals with sources (e.g. deposition) and sinks (e.g. plant uptake) of
!   soil C and N.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in BIOGEOCHEMISTRY
!
! Code Description:
!   Language: Fortran 90.
!-----------------------------------------------------------------------------

USE ancil_info, ONLY:                                                         &
  ! imported scalars
  dim_cs1, nz_soilc => dim_cslayer

USE ecosse_utils_mod, ONLY:                                                   &
  ! imported parameters
  nt

USE jules_soil_ecosse_mod, ONLY:                                              &
  ! imported scalars
  l_soil_N

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

PRIVATE
PUBLIC ecosse_source_sink

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                       &
  ModuleName = 'ECOSSE_SOURCE_SINK_MOD'

CONTAINS

!#############################################################################
!#############################################################################
  
SUBROUTINE ecosse_source_sink( land_pts, vs_pts, convert_veg_input,           &
                         l_add_plant_inputs, vs_index,                        &
                         plant_input_c_dpm, plant_input_c_rpm,                &
                         plant_input_n_dpm, plant_input_n_rpm,                &
                         c_dpm, c_rpm, n_dpm, n_rpm )

! Description:
!   Deals with sources (e.g. deposition) and sinks (e.g. plant uptake) of
!   soil C and N.

USE jules_surface_types_mod, ONLY:                                            &
  ! imported scalars
  npft, ntype

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! The number of land points.
  vs_pts
    ! The number of points with veg and/or soil.

REAL, INTENT(IN) ::                                                           &
  convert_veg_input
    ! Timestep for plant inputs as a fraction of veg model timestep.

LOGICAL, INTENT(IN) ::                                                        &
  l_add_plant_inputs
    ! Flag indicating if plant inputs are to be added this timestep.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  vs_index(land_pts)
    ! Indices of veg/soil points.

REAL, INTENT(IN) ::                                                           &
  plant_input_c_dpm(land_pts,nz_soilc),                                       &
    ! Carbon added to DPM pool by plant inputs (kg m-2).
  plant_input_c_rpm(land_pts,nz_soilc),                                       &
    ! Carbon added to RPM pool by plant inputs (kg m-2).
  plant_input_n_dpm(land_pts,nz_soilc),                                       &
    ! Nitrogen added to DPM pool by plant inputs (kg m-2).
  plant_input_n_rpm(land_pts,nz_soilc)
    ! Nitrogen added to RPM pool by plant inputs (kg m-2).

!-----------------------------------------------------------------------------
! Array arguments with intent(inout).
! Note these are dimensioned with nt but at present only a single time level
! is required. All nt levels are provided for future use.
!-----------------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                        &
  c_dpm(land_pts,nz_soilc,nt),                                                &
    ! C in decomposable plant material (kg m-2).
  c_rpm(land_pts,nz_soilc,nt),                                                &
    ! C in resistant plant material (kg m-2).
  n_dpm(land_pts,nz_soilc,nt),                                                &
    ! N in decomposable plant material (kg m-2).
  n_rpm(land_pts,nz_soilc,nt)
    ! N in resistant plant material (kg m-2).

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ECOSSE_SOURCE_SINK'

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Remove plant uptake of N.
! Not done yet!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Add fixed N to soil, if calculated by TRIFFID.
! Not done yet!
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Add plant litterfall inputs when required.
!-----------------------------------------------------------------------------
IF ( l_add_plant_inputs ) THEN

  CALL add_plant_inputs( land_pts, vs_pts, vs_index,                          &
                         plant_input_c_dpm, plant_input_c_rpm,                &
                         plant_input_n_dpm, plant_input_n_rpm,                &
                         c_dpm(:,:,nt), c_rpm(:,:,nt),                        &
                         n_dpm(:,:,nt), n_rpm(:,:,nt) )

  ! If nt=1, our copy of the prognostic variables can potentially now
  ! include values <= 0 (if the litter flux into a small pool was < 0). This
  ! would likely need to be dealt with here.
END IF

!-----------------------------------------------------------------------------
! Add inputs from atmospheric deposition.
! Not done yet!
!-----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ecosse_source_sink
    
!#############################################################################
!#############################################################################

SUBROUTINE add_plant_inputs( land_pts, vs_pts,vs_index,                       &
               plant_input_c_dpm, plant_input_c_rpm,                          &
               plant_input_n_dpm, plant_input_n_rpm,                          &
               c_dpm, c_rpm, n_dpm, n_rpm )

IMPLICIT NONE

! Description:
!   Adds plant inputs of organic matter to soil.
!   Note that this subroutine can, at least in theory, result in a soil store
!   <0 if a store is depleted because the veg model has calculated
!   "litterfall"<0. That will be dealt with by future developments.

!-----------------------------------------------------------------------------
! Scalar arguments with intent(in).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                        &
  land_pts,                                                                   &
    ! The number of land points.
  vs_pts
    ! The number of points with veg and/or soil.
    ! The number of land points.

!-----------------------------------------------------------------------------
! Array arguments with intent(in)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN)  :: vs_index(land_pts)
       ! Indices of veg/soil points.

REAL, INTENT(IN) ::                                                           &
  plant_input_c_dpm(land_pts,nz_soilc),                                       &
    ! Carbon added to DPM pool by plant inputs (kg m-2).
  plant_input_c_rpm(land_pts,nz_soilc),                                       &
    ! Carbon added to RPM pool by plant inputs (kg m-2).
  plant_input_n_dpm(land_pts,nz_soilc),                                       &
    ! Nitrogen added to DPM pool by plant inputs (kg m-2).
  plant_input_n_rpm(land_pts,nz_soilc)
    ! Nitrogen added to RPM pool by plant inputs (kg m-2).

!-----------------------------------------------------------------------------
! Array arguments with intent(inout).
!-----------------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                        &
  c_dpm(land_pts,nz_soilc),                                                   &
    ! C in DPM (kg m-2).
  c_rpm(land_pts,nz_soilc),                                                   &
    ! C in RPM (kg m-2).
  n_dpm(land_pts,nz_soilc),                                                   &
    ! N in DPM (kg m-2).
  n_rpm(land_pts,nz_soilc) 
    ! N in RPM (kg m-2).

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ADD_PLANT_INPUTS'

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER :: i, j   ! work

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Add plant material to soil pools.
!-----------------------------------------------------------------------------
! Carbon.
DO j = 1,vs_pts
  i = vs_index(j)
  c_dpm(i,:) = c_dpm(i,:) + plant_input_c_dpm(i,:)
  c_rpm(i,:) = c_rpm(i,:) + plant_input_c_rpm(i,:)
END DO

! Nitrogen.
IF ( l_soil_N ) THEN
  DO j = 1,vs_pts
    i = vs_index(j)
    n_dpm(i,:) = n_dpm(i,:) + plant_input_n_dpm(i,:)
    n_rpm(i,:) = n_rpm(i,:) + plant_input_n_rpm(i,:)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE add_plant_inputs
  
!#############################################################################
!#############################################################################

END MODULE ecosse_source_sink_mod
#endif
