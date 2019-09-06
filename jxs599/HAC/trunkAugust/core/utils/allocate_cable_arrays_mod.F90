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

USE cable_types_mod,      ONLY: mp, met, rad, soil, ssnow, canopy, veg
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

CALL allocate_met_type(met, mp)
CALL allocate_radiation_type(rad, mp)
CALL allocate_soil_parameter_type(soil, mp)
CALL allocate_soil_snow_type(ssnow, mp)
CALL allocate_canopy_type(canopy, mp)
CALL allocate_veg_parameter_type(veg, mp)

END SUBROUTINE allocate_cable_arrays


SUBROUTINE allocate_radiation_type(var, mp)

USE cable_types_mod,            ONLY: radiation_type
USE cable_other_constants_mod,  ONLY: nrb, mf, swb

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for radiation.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(radiation_type), INTENT(INOUT) :: var
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_RADIATION_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% albedo(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% extkb(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% extkd(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% rhocdf(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% reffdf(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% reffbm(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% extkbm(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% extkdm(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% cexpkbm(mp,swb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% cexpkdm(mp,swb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% fbeam(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% rhocbm(mp,nrb), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%albedo(:,:) = 0.0
  var%extkb(:) = 0.0
  var%rhocdf(:,:) = 0.0
  var%reffdf(:,:) = 0.0
  var%reffbm(:,:) = 0.0
  var%extkbm(:,:) = 0.0
  var%extkdm(:,:) = 0.0
  var%cexpkbm(:,:) = 0.0
  var%fbeam(:,:) = 0.0
  var%rhocbm(:,:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model radiation arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_radiation_type


SUBROUTINE allocate_veg_parameter_type(var, mp)

USE cable_types_mod,            ONLY: veg_parameter_type
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
ALLOCATE( var%refl(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var%taul(mp,nrb), stat = error )
error_sum = error_sum + error

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

USE cable_types_mod,            ONLY: soil_parameter_type
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


SUBROUTINE allocate_soil_snow_type(var, mp)

USE cable_types_mod,              ONLY: soil_snow_type, snow_rho1l,           &
                                        isnow_flg3l, snow_age
USE cable_other_constants_mod,    ONLY: nrb, msn
USE jules_surface_types_mod,      ONLY: ntype
USE ancil_info,                   ONLY: land_pts
USE jules_soil_mod,               ONLY: sm_levels

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for the soil snow type.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(soil_snow_type), INTENT(INOUT) :: var
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_SOIL_SNOW_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% albsoilsn(mp,nrb), stat = error )
error_sum = error_sum + error
ALLOCATE( var% isflag(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% osnowd(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% snage(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% snowd(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% ssdnn(mp), stat = error )
error_sum = error_sum + error
ALLOCATE( var% tgg(mp,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( var% tggsn(mp,msn), stat = error )
error_sum = error_sum + error
ALLOCATE( var% wb(mp,sm_levels), stat = error )
error_sum = error_sum + error
ALLOCATE( var% wbice(mp,sm_levels), stat = error )
error_sum = error_sum + error

! Allocate snow variables normally obtained from UM in coupled runs
ALLOCATE (snow_rho1l(land_pts, ntype), stat = error )
error_sum = error_sum + error
ALLOCATE (isnow_flg3l(land_pts, ntype), stat = error )
error_sum = error_sum + error
ALLOCATE (snow_age(land_pts, ntype), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%albsoilsn(:,:) = 0.0
  var%isflag(:) = 0
  var%osnowd(:) = 0.0
  var%snage(:) = 0.0
  var%snowd(:) = 0.0
  var%ssdnn(:) = 0.0
  var%tgg(:,:) = 0.0
  var%tggsn(:,:) = 0.0
  var%wb(:,:) = 0.0
  var%wbice(:,:) = 0.0
  snow_rho1l(:,:) = 0.0
  isnow_flg3l(:,:) = 0.0
  snow_age(:,:) = 0.0

ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model soil snow arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_soil_snow_type


SUBROUTINE allocate_canopy_type(var, mp)

USE jules_soil_mod,             ONLY: sm_levels   ! number of soil levels
USE cable_types_mod,            ONLY: canopy_type
USE cable_other_constants_mod,  ONLY: mf, niter, n_assim_rates

!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for the canopy type.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(canopy_type), INTENT(INOUT) :: var
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_CANOPY_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( var% vlaiw(mp), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%vlaiw(:) = 0.0

ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model canopy arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_canopy_type


SUBROUTINE allocate_met_type(var, mp)

USE cable_types_mod,              ONLY: met_type
USE cable_other_constants_mod,    ONLY: swb

!-----------------------------------------------------------------------------
! Description:
!   Allocates the forcing variables.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

TYPE(met_type), INTENT(INOUT) :: var
INTEGER, INTENT(IN) :: mp

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_MET_TYPE'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE ( var % fsd(mp,swb), stat = error )
error_sum = error_sum + error
ALLOCATE ( var % tk(mp), stat = error )
error_sum = error_sum + error
ALLOCATE ( var % coszen(mp), stat = error )
error_sum = error_sum + error

IF (error_sum == 0) THEN
  var%fsd(:,:) = 0.0
  var%tk(:) = 0.0
  var%coszen(:) = 0.0
ELSE
  IF ( error_sum /= 0 )                                                       &
    CALL ereport("allocate_cable_arrays", errcode,                            &
                 "Error allocating CABLE model met arrays")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE allocate_met_type

END MODULE allocate_cable_arrays_mod
