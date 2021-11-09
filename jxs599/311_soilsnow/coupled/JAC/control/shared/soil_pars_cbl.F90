MODULE soil_params_type_mod_cbl

IMPLICIT NONE
!Public subroutines
PUBLIC :: alloc_soil_params_type_cbl
PUBLIC :: dealloc_soil_params_type_cbl
PUBLIC :: assoc_soil_params_type_cbl
PUBLIC :: nullify_soil_params_type_cbl
PUBLIC :: init_soil_cbl
!Public data 
PUBLIC :: soil_params_data_type_cbl
PUBLIC :: soil_params_type_cbl
PRIVATE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOIL_PARAMS_CBL_VARS_MOD'
!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------
! Soil parameters used:
TYPE soil_params_data_type_cbl

  INTEGER, ALLOCATABLE ::                                                     &
    isoilm(:)     ! integer soil type

  REAL, ALLOCATABLE ::                                                        &
    bch(:),     & ! parameter b in Campbell equation
    c3(:),      & ! c3 drainage coeff (fraction)
    clay(:),    & ! fraction of soil which is clay
    css(:),     & ! soil specific heat capacity [kJ/kg/K]
    hsbh(:),    & ! difsat * etasat (=hyds*abs(sucs)*bch)
    hyds(:),    & ! hydraulic conductivity @ saturation [m/s], Ksat
    i2bp3(:),   & ! par. one in K vis suction (=nint(bch)+2)
    ibp2(:),    & ! par. two in K vis suction (fn of pbch)
    rhosoil(:), & ! soil density [kg/m3]
    sand(:),    & ! fraction of soil which is sand
    sfc(:),     & ! vol H2O @ field capacity
    silt(:),    & ! fraction of soil which is silt
    ssat(:),    & ! vol H2O @ saturation
    sucs(:),    & ! suction at saturation (m)
    swilt(:),   & ! vol H2O @ wilting
    zse(:),     & ! thickness of each soil layer (1=top) in m
    zshh(:),    & ! distance between consecutive layer midpoints (m)
                  ! vars intro for Ticket #27
    soilcol(:), & ! keep color for all patches/tiles
    albsoilf(:),& ! soil reflectance
    cnsd(:),    & ! thermal conductivity of dry soil [W/m/K]
    pwb_min(:)    ! working variable (swilt/ssat)**ibp2

  REAL, ALLOCATABLE ::                                                        &
    heat_cap_lower_limit(:,:),                                                &
    zse_vec(:,:),                                                             &
    css_vec(:,:),                                                             &
    cnsd_vec(:,:),                                                            &
    albsoil(:,:)    ! soil reflectance (2nd dim. BP 21Oct2009)
     
  !parameters for GW module that vary with soil layer
  REAL, ALLOCATABLE ::                                                        &
    sucs_vec(:,:),    & !psi at saturation in [mm]
    hyds_vec(:,:),    & !saturated hydraulic conductivity  [mm/s]
    bch_vec(:,:),     & !C and H B [none]
    clay_vec(:,:),    & !fraction of soil that is clay [frac]
    sand_vec(:,:),    & !fraction of soil that is sand [frac]
    silt_vec(:,:),    & !fraction of soil that is silt [frac]
    org_vec(:,:),     & !fration of soil made of organic soils [frac]
    rhosoil_vec(:,:), & !soil density  [kg/m3]
    ssat_vec(:,:),    & !volumetric water content at saturation [mm3/mm3]
    watr(:,:),        & !residual water content of the soil [mm3/mm3]
    sfc_vec(:,:),     & !field capcacity (hk = 1 mm/day)
    swilt_vec(:,:)      ! wilting point (hk = 0.02 mm/day)

  REAL, ALLOCATABLE ::                                                        &
    drain_dens(:),  & !  drainage density ( mean dist to rivers/streams )
    elev(:),        & !elevation above sea level
    elev_std(:),    & !elevation above sea level
    slope(:),       & !mean slope of grid cell
    slope_std(:)      !stddev of grid cell slope

  !parameters for GW module for the aquifer
  REAL, ALLOCATABLE ::                                                        &
    GWsucs_vec(:),  &  !head in the aquifer [mm]
    GWhyds_vec(:),  &  !saturated hydraulic conductivity of the aquifer [mm/s]
    GWbch_vec(:),   & !clapp and horn b of the aquifer   [none]
    GWssat_vec(:),  & !saturated water content of the aquifer [mm3/mm3]
    GWwatr(:),      & !residual water content of the aquifer [mm3/mm3]
    GWz(:),         & !node depth of the aquifer    [m]
    GWdz(:),        & !thickness of the aquifer   [m]
    GWrhosoil_vec(:)    !density of the aquifer substrate [kg/m3]

  ! Additional SLI parameters
  INTEGER, ALLOCATABLE :: nhorizons(:) ! number of soil horizons
  INTEGER, ALLOCATABLE :: ishorizon(:,:) ! horizon number 1:nhorizons
  REAL, ALLOCATABLE ::                                                        &
    clitt(:),      & ! litter (tC/ha)
    zeta(:),       & ! macropore parameter
    fsatmax(:)       ! variably saturated area parameter

END TYPE soil_params_data_type_cbl

! Pointers to Soil parameters used:
TYPE soil_params_type_cbl

  INTEGER, POINTER  ::                                                        &
    isoilm(:)     ! integer soil type

  REAL, POINTER ::                                                            &
    bch(:),     & ! parameter b in Campbell equation
    c3(:),      & ! c3 drainage coeff (fraction)
    clay(:),    & ! fraction of soil which is clay
    css(:),     & ! soil specific heat capacity [kJ/kg/K]
    hsbh(:),    & ! difsat * etasat (=hyds*abs(sucs)*bch)
    hyds(:),    & ! hydraulic conductivity @ saturation [m/s], Ksat
    i2bp3(:),   & ! par. one in K vis suction (=nint(bch)+2)
    ibp2(:),    & ! par. two in K vis suction (fn of pbch)
    rhosoil(:), & ! soil density [kg/m3]
    sand(:),    & ! fraction of soil which is sand
    sfc(:),     & ! vol H2O @ field capacity
    silt(:),    & ! fraction of soil which is silt
    ssat(:),    & ! vol H2O @ saturation
    sucs(:),    & ! suction at saturation (m)
    swilt(:),   & ! vol H2O @ wilting
    zse(:),     & ! thickness of each soil layer (1=top) in m
    zshh(:),    & ! distance between consecutive layer midpoints (m)
                  ! vars intro for Ticket #27
    soilcol(:), & ! keep color for all patches/tiles
    albsoilf(:),& ! soil reflectance
    cnsd(:),    & ! thermal conductivity of dry soil [W/m/K]
    pwb_min(:)    ! working variable (swilt/ssat)**ibp2

  REAL, POINTER ::                                                            &
    heat_cap_lower_limit(:,:),                                                &
    zse_vec(:,:),                                                             &
    css_vec(:,:),                                                             &
    cnsd_vec(:,:),                                                            &
    albsoil(:,:)    ! soil reflectance (2nd dim. BP 21Oct2009)
     
  !parameters for GW module that vary with soil layer
  REAL, POINTER ::                                                            &
    sucs_vec(:,:),    & !psi at saturation in [mm]
    hyds_vec(:,:),    & !saturated hydraulic conductivity  [mm/s]
    bch_vec(:,:),     & !C and H B [none]
    clay_vec(:,:),    & !fraction of soil that is clay [frac]
    sand_vec(:,:),    & !fraction of soil that is sand [frac]
    silt_vec(:,:),    & !fraction of soil that is silt [frac]
    org_vec(:,:),     & !fration of soil made of organic soils [frac]
    rhosoil_vec(:,:), & !soil density  [kg/m3]
    ssat_vec(:,:),    & !volumetric water content at saturation [mm3/mm3]
    watr(:,:),        & !residual water content of the soil [mm3/mm3]
    sfc_vec(:,:),     & !field capcacity (hk = 1 mm/day)
    swilt_vec(:,:)      ! wilting point (hk = 0.02 mm/day)

  REAL, POINTER ::                                                            &
    drain_dens(:),  & !  drainage density ( mean dist to rivers/streams )
    elev(:),        & !elevation above sea level
    elev_std(:),    & !elevation above sea level
    slope(:),       & !mean slope of grid cell
    slope_std(:)      !stddev of grid cell slope

  !parameters for GW module for the aquifer
  REAL, POINTER ::                                                            &
    GWsucs_vec(:),  &  !head in the aquifer [mm]
    GWhyds_vec(:),  &  !saturated hydraulic conductivity of the aquifer [mm/s]
    GWbch_vec(:),   & !clapp and horn b of the aquifer   [none]
    GWssat_vec(:),  & !saturated water content of the aquifer [mm3/mm3]
    GWwatr(:),      & !residual water content of the aquifer [mm3/mm3]
    GWz(:),         & !node depth of the aquifer    [m]
    GWdz(:),        & !thickness of the aquifer   [m]
    GWrhosoil_vec(:)    !density of the aquifer substrate [kg/m3]

  ! Additional SLI parameters
  INTEGER, POINTER :: nhorizons(:) ! number of soil horizons
  INTEGER, POINTER :: ishorizon(:,:) ! horizon number 1:nhorizons
  REAL, POINTER  ::                                                           &
    clitt(:),      & ! litter (tC/ha)
    zeta(:),       & ! macropore parameter
    fsatmax(:)       ! variably saturated area parameter

END TYPE soil_params_type_cbl

CONTAINS

SUBROUTINE alloc_soil_params_type_cbl( mp, var)

USE ereport_mod,              ONLY: ereport

USE jules_model_environment_mod,  ONLY: lsm_id, cable
USE grid_constants_mod_cbl, ONLY: nrb,  & !total # rad. "bands"
                                  nsl     ! # soil layers

!-----------------------------------------------------------------------------
! Description:
!   Allocates soil parameter arrays.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) :: mp
TYPE(soil_params_data_type_cbl), INTENT(INOUT) :: var

INTEGER ::                                                                    &
  error        = 0,                                                           &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                           &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode      = 101
                       ! Variable to use in error report

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_PARAMETER_TYPE'
!End of header

IF ( lsm_id == cable ) THEN
  ALLOCATE( var% isoilm(mp),        stat = error )
  ALLOCATE( var% albsoil(mp, nrb),  stat = error )
  ALLOCATE( var% albsoilf(mp),      stat = error )
  ALLOCATE( var% soilcol(mp),       stat = error )
  ALLOCATE( var% bch(mp),           stat = error )    
  ALLOCATE( var% c3(mp),            stat = error )
  ALLOCATE( var% clay(mp),          stat = error )
  ALLOCATE( var% css(mp),           stat = error )
  ALLOCATE( var% hsbh(mp),          stat = error )
  ALLOCATE( var% hyds(mp),          stat = error )
  ALLOCATE( var% i2bp3(mp),         stat = error )
  ALLOCATE( var% ibp2(mp),          stat = error )
  ALLOCATE( var% rhosoil(mp),       stat = error )
  ALLOCATE( var% sand(mp),          stat = error )
  ALLOCATE( var% sfc(mp),           stat = error )
  ALLOCATE( var% silt(mp),          stat = error )
  ALLOCATE( var% ssat(mp),          stat = error )
  ALLOCATE( var% sucs(mp),          stat = error )
  ALLOCATE( var% swilt(mp),         stat = error )
  ALLOCATE( var% zse(nsl),          stat = error )
  ALLOCATE( var% zshh(nsl+1),       stat = error )
  ALLOCATE( var% cnsd(mp),          stat = error )
  ALLOCATE( var% pwb_min(mp),       stat = error )
  !Aquifer properties               
  ALLOCATE( var%GWhyds_vec(mp),     stat = error )
  ALLOCATE( var%GWsucs_vec(mp),     stat = error )
  ALLOCATE( var%GWbch_vec(mp),      stat = error )
  ALLOCATE( var%GWssat_vec(mp),     stat = error )
  ALLOCATE( var%GWwatr(mp),         stat = error )
  ALLOCATE( var%GWz(mp),            stat = error )
  ALLOCATE( var%GWdz(mp),           stat = error )
  ALLOCATE( var%GWrhosoil_vec(mp),  stat = error )
  !soil properties (vary by layer)
  ALLOCATE( var% zse_vec(mp,nsl),   stat = error )
  ALLOCATE( var% css_vec(mp,nsl),   stat = error )
  ALLOCATE( var% cnsd_vec(mp,nsl),  stat = error )
  ALLOCATE( var%hyds_vec(mp,nsl),   stat = error )
  ALLOCATE( var%sucs_vec(mp,nsl),   stat = error )
  ALLOCATE( var%bch_vec(mp,nsl),    stat = error )
  ALLOCATE( var%ssat_vec(mp,nsl),   stat = error )
  ALLOCATE( var%watr(mp,nsl),       stat = error )
  ALLOCATE( var%sfc_vec(mp,nsl),    stat = error )
  ALLOCATE( var%swilt_vec(mp,nsl),  stat = error )
  ALLOCATE( var%sand_vec(mp,nsl),   stat = error )
  ALLOCATE( var%clay_vec(mp,nsl),   stat = error )
  ALLOCATE( var%silt_vec(mp,nsl),   stat = error )
  ALLOCATE( var%org_vec(mp,nsl),    stat = error )
  ALLOCATE( var%rhosoil_vec(mp,nsl),stat = error )
  ALLOCATE( var% heat_cap_lower_limit(mp,nsl),   stat = error )
                                                 
  ALLOCATE( var%drain_dens(mp),     stat = error )
  ALLOCATE( var%elev(mp),           stat = error )
  ALLOCATE( var%elev_std(mp),       stat = error )
  ALLOCATE( var%slope(mp),          stat = error )
  ALLOCATE( var%slope_std(mp),      stat = error )
  
  ! Allocate variables for SLI soil model:        
  ALLOCATE ( var % nhorizons(mp),     stat = error )
  ALLOCATE ( var % ishorizon(mp,nsl), stat = error )  
  ALLOCATE ( var % clitt(mp),         stat = error )
  ALLOCATE ( var % zeta(mp),          stat = error )
  ALLOCATE ( var % fsatmax(mp),       stat = error )
  ALLOCATE ( var % swilt_vec(mp,nsl), stat = error )
  ALLOCATE ( var % ssat_vec(mp,nsl),  stat = error )
  ALLOCATE ( var % sfc_vec(mp,nsl),   stat = error )
ELSE ! lsm_id==JULES
  ALLOCATE( var% isoilm(1),        stat = error )
  ALLOCATE( var% albsoil(1, 1),  stat = error )
  ALLOCATE( var% albsoilf(1),      stat = error )
  ALLOCATE( var% soilcol(1),       stat = error )
  ALLOCATE( var% bch(1),           stat = error )    
  ALLOCATE( var% c3(1),            stat = error )
  ALLOCATE( var% clay(1),          stat = error )
  ALLOCATE( var% css(1),           stat = error )
  ALLOCATE( var% hsbh(1),          stat = error )
  ALLOCATE( var% hyds(1),          stat = error )
  ALLOCATE( var% i2bp3(1),         stat = error )
  ALLOCATE( var% ibp2(1),          stat = error )
  ALLOCATE( var% rhosoil(1),       stat = error )
  ALLOCATE( var% sand(1),          stat = error )
  ALLOCATE( var% sfc(1),           stat = error )
  ALLOCATE( var% silt(1),          stat = error )
  ALLOCATE( var% ssat(1),          stat = error )
  ALLOCATE( var% sucs(1),          stat = error )
  ALLOCATE( var% swilt(1),         stat = error )
  ALLOCATE( var% zse(1),          stat = error )
  ALLOCATE( var% zshh(1),       stat = error )
  ALLOCATE( var% cnsd(1),          stat = error )
  ALLOCATE( var% pwb_min(1),       stat = error )
  !Aquifer properties               
  ALLOCATE( var%GWhyds_vec(1),     stat = error )
  ALLOCATE( var%GWsucs_vec(1),     stat = error )
  ALLOCATE( var%GWbch_vec(1),      stat = error )
  ALLOCATE( var%GWssat_vec(1),     stat = error )
  ALLOCATE( var%GWwatr(1),         stat = error )
  ALLOCATE( var%GWz(1),            stat = error )
  ALLOCATE( var%GWdz(1),           stat = error )
  ALLOCATE( var%GWrhosoil_vec(1),  stat = error )
  !soil properties (vary by layer)
  ALLOCATE( var% zse_vec,   stat = error )
  ALLOCATE( var% css_vec(1,1),   stat = error )
  ALLOCATE( var% cnsd_vec(1,1),  stat = error )
  ALLOCATE( var%hyds_vec(1,1),   stat = error )
  ALLOCATE( var%sucs_vec(1,1),   stat = error )
  ALLOCATE( var%bch_vec(1,1),    stat = error )
  ALLOCATE( var%ssat_vec(1,1),   stat = error )
  ALLOCATE( var%watr(1,1),       stat = error )
  ALLOCATE( var%sfc_vec(1,1),    stat = error )
  ALLOCATE( var%swilt_vec(1,1),  stat = error )
  ALLOCATE( var%sand_vec(1,1),   stat = error )
  ALLOCATE( var%clay_vec(1,1),   stat = error )
  ALLOCATE( var%silt_vec(1,1),   stat = error )
  ALLOCATE( var%org_vec(1,1),    stat = error )
  ALLOCATE( var%rhosoil_vec(1,1),stat = error )
  ALLOCATE( var% heat_cap_lower_limit(1,1),   stat = error )
                                                 
  ALLOCATE( var%drain_dens(1),     stat = error )
  ALLOCATE( var%elev(1),           stat = error )
  ALLOCATE( var%elev_std(1),       stat = error )
  ALLOCATE( var%slope(1),          stat = error )
  ALLOCATE( var%slope_std(1),      stat = error )
  
  ! Allocate variables for SLI soil model:        
  ALLOCATE ( var % nhorizons(1),     stat = error )
  ALLOCATE ( var % ishorizon(1,1), stat = error )  
  ALLOCATE ( var % clitt(1),         stat = error )
  ALLOCATE ( var % zeta(1),          stat = error )
  ALLOCATE ( var % fsatmax(1),       stat = error )
  ALLOCATE ( var % swilt_vec(1,1), stat = error )
  ALLOCATE ( var % ssat_vec(1,1),  stat = error )
  ALLOCATE ( var % sfc_vec(1,1),   stat = error )
END IF

!something needs to be done with these Zeroes 
IF (error_sum == 0) THEN
  var%isoilm(:) = 0
  var%albsoil(:,:) = 0.0
  var%albsoilf(:) = 0.0
  var%soilcol(:) = 0.0
ELSE
  CALL ereport("allocate_cable_arrays", errcode,                            &
               "Error allocating CABLE model soil parameter arrays")
END IF
!something needs to be done with these blind clobbers
var%GWwatr(:) = 0.05
var%watr(:,:) = 0.05
!ALLOCATE ( var % swilt_vec(mp,nsl) )
!ALLOCATE ( var % ssat_vec(mp,nsl) )
!ALLOCATE ( var % sfc_vec(mp,nsl) )

RETURN

END SUBROUTINE alloc_soil_params_type_cbl



SUBROUTINE dealloc_soil_params_type_cbl( var )
  !Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(soil_params_data_type_cbl), INTENT(INOUT) :: var
  
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='dealloc_soil_params_type_cbl'
  
!End of header
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE(   var% isoilm       )
DEALLOCATE(   var% albsoil      )
DEALLOCATE(   var% albsoilf     )
DEALLOCATE(   var% soilcol      )
DEALLOCATE(   var% bch          )
DEALLOCATE(   var% c3           )
DEALLOCATE(   var% clay         )
DEALLOCATE(   var% css          )
DEALLOCATE(   var% hsbh         )
DEALLOCATE(   var% hyds         )
DEALLOCATE(   var% i2bp3        )
DEALLOCATE(   var% ibp2         )
DEALLOCATE(   var% rhosoil      )
DEALLOCATE(   var% sand         )
DEALLOCATE(   var% sfc          )
DEALLOCATE(   var% silt         )
DEALLOCATE(   var% ssat         )
DEALLOCATE(   var% sucs         )
DEALLOCATE(   var% swilt        )
DEALLOCATE(   var% zse          )
DEALLOCATE(   var% zshh         )
DEALLOCATE(   var% cnsd         )
DEALLOCATE(   var% pwb_min      )
DEALLOCATE(   var%GWhyds_vec    )
DEALLOCATE(   var%GWsucs_vec    )
DEALLOCATE(   var%GWbch_vec     )
DEALLOCATE(   var%GWssat_vec    )
DEALLOCATE(   var%GWwatr        )
DEALLOCATE(   var%GWz           )
DEALLOCATE(   var%GWdz          )
DEALLOCATE(   var%GWrhosoil_vec )
DEALLOCATE(   var% zse_vec      )
DEALLOCATE(   var% css_vec      )
DEALLOCATE(   var% cnsd_vec     )
DEALLOCATE(   var%hyds_vec      )
DEALLOCATE(   var%sucs_vec      )
DEALLOCATE(   var%bch_vec       )
DEALLOCATE(   var%ssat_vec      )
DEALLOCATE(   var%watr          )
DEALLOCATE(   var%sfc_vec       )
DEALLOCATE(   var%swilt_vec     )
DEALLOCATE(   var%sand_vec      )
DEALLOCATE(   var%clay_vec      )
DEALLOCATE(   var%silt_vec      )
DEALLOCATE(   var%org_vec       )
DEALLOCATE(   var%rhosoil_vec   )
DEALLOCATE(   var%drain_dens    )
DEALLOCATE(   var%elev          )
DEALLOCATE(   var%elev_std      )
DEALLOCATE(   var%slope         )
DEALLOCATE(   var%slope_std     )
DEALLOCATE(   var % nhorizons   )
DEALLOCATE(   var % ishorizon   )
DEALLOCATE(   var % clitt       )
DEALLOCATE(   var % zeta        )
DEALLOCATE(   var % fsatmax     )
DEALLOCATE(   var % swilt_vec   )
DEALLOCATE(   var % ssat_vec    )
DEALLOCATE(   var % sfc_vec     )
DEALLOCATE(   var% heat_cap_lower_limit )
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE dealloc_soil_params_type_cbl


SUBROUTINE assoc_soil_params_type_cbl( var, var_ptr )

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(soil_params_data_type_cbl), INTENT(INOUT), TARGET  :: var
TYPE(soil_params_type_cbl), INTENT(INOUT) :: var_ptr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='assoc_soil_params_type_cbl'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL nullify_soil_params_type_cbl(var_ptr)

  var_ptr% isoilm       => var% isoilm
  var_ptr% albsoil      => var% albsoil
  var_ptr% albsoilf     => var% albsoilf
  var_ptr% soilcol      => var% soilcol
  var_ptr% bch          => var% bch
  var_ptr% c3           => var% c3 
  var_ptr% clay         => var% clay 
  var_ptr% css          => var% css 
  var_ptr% hsbh         => var% hsbh 
  var_ptr% hyds         => var% hyds 
  var_ptr% i2bp3        => var% i2bp3 
  var_ptr% ibp2         => var% ibp2 
  var_ptr% rhosoil      => var% rhosoil 
  var_ptr% sand         => var% sand 
  var_ptr% sfc          => var% sfc 
  var_ptr% silt         => var% silt 
  var_ptr% ssat         => var% ssat 
  var_ptr% sucs         => var% sucs 
  var_ptr% swilt        => var% swilt 
  var_ptr% zse          => var% zse
  var_ptr% zshh         => var% zshh
  var_ptr% cnsd         => var% cnsd 
  var_ptr% pwb_min      => var% pwb_min 
  var_ptr%GWhyds_vec    => var%GWhyds_vec 
  var_ptr%GWsucs_vec    => var%GWsucs_vec 
  var_ptr%GWbch_vec     => var%GWbch_vec 
  var_ptr%GWssat_vec    => var%GWssat_vec 
  var_ptr%GWwatr        => var%GWwatr 
  var_ptr%GWz           => var%GWz 
  var_ptr%GWdz          => var%GWdz 
  var_ptr%GWrhosoil_vec => var%GWrhosoil_vec 
  var_ptr% zse_vec      => var% zse_vec 
  var_ptr% css_vec      => var% css_vec
  var_ptr% cnsd_vec     => var% cnsd_vec
  var_ptr%hyds_vec      => var%hyds_vec 
  var_ptr%sucs_vec      => var%sucs_vec 
  var_ptr%bch_vec       => var%bch_vec 
  var_ptr%ssat_vec      => var%ssat_vec 
  var_ptr%watr          => var%watr 
  var_ptr%sfc_vec       => var%sfc_vec 
  var_ptr%swilt_vec     => var%swilt_vec 
  var_ptr%sand_vec      => var%sand_vec 
  var_ptr%clay_vec      => var%clay_vec 
  var_ptr%silt_vec      => var%silt_vec 
  var_ptr%org_vec       => var%org_vec 
  var_ptr%rhosoil_vec   => var%rhosoil_vec
  var_ptr%drain_dens    => var%drain_dens 
  var_ptr%elev          => var%elev 
  var_ptr%elev_std      => var%elev_std 
  var_ptr%slope         => var%slope 
  var_ptr%slope_std     => var%slope_std 
  var_ptr % nhorizons   => var % nhorizons 
  var_ptr % ishorizon   => var % ishorizon 
  var_ptr % clitt       => var % clitt
  var_ptr % zeta        => var % zeta 
  var_ptr % fsatmax     => var % fsatmax 
  var_ptr % swilt_vec   => var % swilt_vec 
  var_ptr % ssat_vec    => var % ssat_vec 
  var_ptr % sfc_vec     => var % sfc_vec 
  var_ptr% heat_cap_lower_limit => var% heat_cap_lower_limit 

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE assoc_soil_params_type_cbl

SUBROUTINE nullify_soil_params_type_cbl(var_ptr)

USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(soil_params_type_cbl), INTENT(INOUT) :: var_ptr
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='nullify_soil_params_type_cbl'
  
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY( var_ptr% isoilm       )
NULLIFY( var_ptr% albsoil      )
NULLIFY( var_ptr% albsoilf     )
NULLIFY( var_ptr% soilcol      ) 
NULLIFY( var_ptr% bch          ) 
NULLIFY( var_ptr% c3           )
NULLIFY( var_ptr% clay         )
NULLIFY( var_ptr% css          )
NULLIFY( var_ptr% hsbh         )
NULLIFY( var_ptr% hyds         )
NULLIFY( var_ptr% i2bp3        )
NULLIFY( var_ptr% ibp2         )
NULLIFY( var_ptr% rhosoil      )
NULLIFY( var_ptr% sand         )
NULLIFY( var_ptr% sfc          )
NULLIFY( var_ptr% silt         )
NULLIFY( var_ptr% ssat         )
NULLIFY( var_ptr% sucs         )
NULLIFY( var_ptr% swilt        )
NULLIFY( var_ptr% zse          )
NULLIFY( var_ptr% zshh         )
NULLIFY( var_ptr% cnsd         )
NULLIFY( var_ptr% pwb_min      )
NULLIFY( var_ptr%GWhyds_vec    )
NULLIFY( var_ptr%GWsucs_vec    )
NULLIFY( var_ptr%GWbch_vec     )
NULLIFY( var_ptr%GWssat_vec    )
NULLIFY( var_ptr%GWwatr        )
NULLIFY( var_ptr%GWz           )
NULLIFY( var_ptr%GWdz          )
NULLIFY( var_ptr%GWrhosoil_vec )
NULLIFY( var_ptr% zse_vec      )
NULLIFY( var_ptr% css_vec      )
NULLIFY( var_ptr% cnsd_vec     )
NULLIFY( var_ptr%hyds_vec      )
NULLIFY( var_ptr%sucs_vec      )
NULLIFY( var_ptr%bch_vec       )
NULLIFY( var_ptr%ssat_vec      )
NULLIFY( var_ptr%watr          )
NULLIFY( var_ptr%sfc_vec       )
NULLIFY( var_ptr%swilt_vec     )
NULLIFY( var_ptr%sand_vec      )
NULLIFY( var_ptr%clay_vec      )
NULLIFY( var_ptr%silt_vec      )
NULLIFY( var_ptr%org_vec       )       
NULLIFY( var_ptr%rhosoil_vec   )       
NULLIFY( var_ptr%drain_dens    )       
NULLIFY( var_ptr%elev          )       
NULLIFY( var_ptr%elev_std      )       
NULLIFY( var_ptr%slope         )       
NULLIFY( var_ptr%slope_std     )       
NULLIFY( var_ptr % nhorizons   )       
NULLIFY( var_ptr % ishorizon   )       
NULLIFY( var_ptr % clitt       )       
NULLIFY( var_ptr % zeta        )       
NULLIFY( var_ptr % fsatmax     )       
NULLIFY( var_ptr % swilt_vec   )
NULLIFY( var_ptr % ssat_vec    )
NULLIFY( var_ptr % sfc_vec     )
NULLIFY( var_ptr% heat_cap_lower_limit )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nullify_soil_params_type_cbl

SUBROUTINE init_soil_cbl( mp, soil, soilin )
USE soilin_pars_mod_cbl,  ONLY: soilin_type
USE ancil_info,           ONLY: nsurft, land_pts

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   init_cables soil panrameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER ::mp
TYPE(soil_params_data_type_cbl), TARGET :: soil
TYPE(soilin_type),               TARGET :: soilin
INTEGER, PARAMETER :: soil_isoilm=2

!jhan:we can bring soil in here and discriminate based on veg%iveg to set PF
!im sure the others are available? bexp etc

!Arbitrarily we SET soil_isoilm=2. Why Eva did this I dont know
!The bulk of soil pars we get from spatially defined vars thru the UM
!We only actually use Permafrost parameters but this is fixed in CABLE 
!interface. The important thing is to send the mp allocated arrays!
soil%css      = soilin%css(soil_isoilm)      
soil%clay     = soilin%clay(soil_isoilm)         
soil%silt     = soilin%silt(soil_isoilm)         
soil%sand     = soilin%sand(soil_isoilm)         
!NO JULES analog. Only ever used as a product, for which there is a JULES var 
soil%rhosoil  = soilin%rhosoil(soil_isoilm) 

RETURN

END SUBROUTINE init_soil_cbl

END MODULE soil_params_type_mod_cbl


