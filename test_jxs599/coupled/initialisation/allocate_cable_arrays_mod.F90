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

SUBROUTINE allocate_cable_arrays()

USE cable_types_mod,      ONLY: mp
USE cable_types_mod,      ONLY: L_tile_pts
USE cable_params_mod,     ONLY: soil, veg
USE ancil_info,           ONLY: surft_pts

!-----------------------------------------------------------------------------
! Description:
!   Allocates the CABLE model arrays using sizes determined during
!   initialisation
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Determine the number of active tiles
mp = SUM(surft_pts)

CALL allocate_soil_parameter_type(soil, mp)
CALL allocate_veg_parameter_type(veg, mp)

END SUBROUTINE allocate_cable_arrays

SUBROUTINE allocate_veg_parameter_type(var, mp)

USE cable_params_mod,            ONLY: veg_parameter_type
USE cable_other_constants_mod,  ONLY: nrb
USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for vegetation parameters.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(veg_parameter_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_VEG_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% iveg(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% hc(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% vlai(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% xfang(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var%refl(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var%taul(mp,nrb), stat = error )
error_sum = error_sum + error

ALLOCATE( var% canst1(mp) )
ALLOCATE( var% dleaf(mp) )
ALLOCATE( var% ejmax(mp) )
!H!    ALLOCATE( var% iLU(mp) )
!H!    ALLOCATE( var% meth(mp) )
ALLOCATE( var% frac4(mp) )
ALLOCATE( var% xalbnir(mp) )
ALLOCATE( var% rp20(mp) )
ALLOCATE( var% rpcoef(mp) )
ALLOCATE( var% rs20(mp) )
ALLOCATE( var% shelrb(mp) )
ALLOCATE( var% vegcf(mp) )
ALLOCATE( var% tminvj(mp) )
!H!    ALLOCATE( var% toptvj(mp) )
ALLOCATE( var% tmaxvj(mp) )
ALLOCATE( var% vbeta(mp) )
ALLOCATE( var% vcmax(mp) )
ALLOCATE( var%extkn(mp) )
ALLOCATE( var%wai(mp) )
!H!    ALLOCATE( var%deciduous(mp) )
    !was nrb(=3), but never uses (:,3) in model
!H!    ALLOCATE( var%vlaimax(mp) )
ALLOCATE( var%a1gs(mp) )
ALLOCATE( var%d0gs(mp) )
ALLOCATE( var%alpha(mp) )
ALLOCATE( var%convex(mp) )
ALLOCATE( var%cfrd(mp) )
ALLOCATE( var%gswmin(mp) )
ALLOCATE( var%conkc0(mp) )
ALLOCATE( var%conko0(mp) )
ALLOCATE( var%ekc(mp) )
ALLOCATE( var%eko(mp) )
ALLOCATE( var% g0(mp) )   ! Ticket #56.
ALLOCATE( var% g1(mp) )   ! Ticket #56.

ALLOCATE( var%froot(mp,sm_levels) )

ALLOCATE ( var % rootbeta(mp) )
!H!    ALLOCATE ( var % gamma(mp) )
!H!    ALLOCATE ( var % F10(mp) )
ALLOCATE ( var % zr(mp) )
ALLOCATE ( var % clitt(mp) )

!H!    ALLOCATE ( var % disturbance_interval(mp,2) )
!H!    ALLOCATE ( var % disturbance_intensity(mp,2) )




IF (error_sum == 0) THEN
  var%iveg(:) = 0
  var%hc(:) = 0.0
  var%vlai(:) = 0.0
  var%refl(:,:) = 0.0
  var%taul(:,:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model veg parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_veg_parameter_type


SUBROUTINE allocate_soil_parameter_type(var, mp)

USE cable_params_mod,            ONLY: soil_parameter_type
USE cable_other_constants_mod,  ONLY: nrb
USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates soil parameter arrays.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(soil_parameter_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_PARAMETER_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% isoilm(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoil(mp, nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% albsoilf(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% soilcol(mp), stat = error )
error_sum = error_sum + error
!H! from cable_define_types
!H!    ALLOCATE( var% bch(mp) )
!H!    ALLOCATE( var% c3(mp) )
!H!    ALLOCATE( var% clay(mp) )
!H!    ALLOCATE( var% css(mp) )
!H!    ALLOCATE( var% hsbh(mp) )
!H!    ALLOCATE( var% hyds(mp) )
!H!    ALLOCATE( var% i2bp3(mp) )
!H!    ALLOCATE( var% ibp2(mp) )
!H!    ALLOCATE( var% rhosoil(mp) )
!H!    ALLOCATE( var% sand(mp) )
!H!    ALLOCATE( var% sfc(mp) )
!H!    ALLOCATE( var% silt(mp) )
!H!    ALLOCATE( var% ssat(mp) )
!H!    ALLOCATE( var% sucs(mp) )
!H!    ALLOCATE( var% swilt(mp) )
!H!    ALLOCATE( var% zse(ms) )
!H!    ALLOCATE( var% zshh(ms+1) )
!H!    ALLOCATE( var% cnsd(mp) )
!H!    ALLOCATE( var% pwb_min(mp) )
!H!    !mrd561
!H!    !MD
!H!    !Aquifer properties
!H!    ALLOCATE( var%GWhyds_vec(mp) )
!H!    ALLOCATE( var%GWsucs_vec(mp) )
!H!    ALLOCATE( var%GWbch_vec(mp) )
!H!    ALLOCATE( var%GWssat_vec(mp) )
!H!    ALLOCATE( var%GWwatr(mp) )
!H!    var%GWwatr(:) = 0.05
!H!    ALLOCATE( var%GWz(mp) )
!H!    ALLOCATE( var%GWdz(mp) )
!H!    ALLOCATE( var%GWrhosoil_vec(mp) )
!H!    !soil properties (vary by layer)
!H!    ALLOCATE( var% zse_vec(mp,ms) )
!H!    ALLOCATE( var% heat_cap_lower_limit(mp,ms) )
!H!    ALLOCATE( var% css_vec(mp,ms) )
!H!    ALLOCATE( var% cnsd_vec(mp,ms) )
!H!    ALLOCATE( var%hyds_vec(mp,ms) )
!H!    ALLOCATE( var%sucs_vec(mp,ms) )
!H!    ALLOCATE( var%bch_vec(mp,ms) )
!H!    ALLOCATE( var%ssat_vec(mp,ms) )
!H!    ALLOCATE( var%watr(mp,ms) )
!H!    var%watr(:,:) = 0.05
!H!    ALLOCATE( var%sfc_vec(mp,ms) )
!H!    ALLOCATE( var%swilt_vec(mp,ms) )
!H!    ALLOCATE( var%sand_vec(mp,ms) )
!H!    ALLOCATE( var%clay_vec(mp,ms) )
!H!    ALLOCATE( var%silt_vec(mp,ms) )
!H!    ALLOCATE( var%org_vec(mp,ms) )
!H!    ALLOCATE( var%rhosoil_vec(mp,ms) )
!H!
!H!    ALLOCATE( var%drain_dens(mp) )
!H!    ALLOCATE( var%elev(mp) )
!H!    ALLOCATE( var%elev_std(mp) )
!H!    ALLOCATE( var%slope(mp) )
!H!    ALLOCATE( var%slope_std(mp) )
!H!
!H!    ! Allocate variables for SLI soil model:
!H!    ALLOCATE ( var % nhorizons(mp) )
!H!    ALLOCATE ( var % ishorizon(mp,ms) )
!H!    ALLOCATE ( var % clitt(mp) )
!H!    ALLOCATE ( var % zeta(mp) )
!H!    ALLOCATE ( var % fsatmax(mp) )
!H!    !ALLOCATE ( var % swilt_vec(mp,ms) )
!H!    !ALLOCATE ( var % ssat_vec(mp,ms) )
!H!    !ALLOCATE ( var % sfc_vec(mp,ms) )
!H!    IF(.NOT.(ASSOCIATED(var % swilt_vec))) ALLOCATE ( var % swilt_vec(mp,ms) )
!H!    IF(.NOT.(ASSOCIATED(var % ssat_vec))) ALLOCATE ( var % ssat_vec(mp,ms) )
!H!    IF(.NOT.(ASSOCIATED(var % sfc_vec))) ALLOCATE ( var % sfc_vec(mp,ms) )



IF (error_sum == 0) THEN
  var%isoilm(:) = 0
  var%albsoil(:,:) = 0.0
  var%albsoilf(:) = 0.0
  var%soilcol(:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model soil parameter arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_soil_parameter_type

END MODULE allocate_cable_arrays_mod
