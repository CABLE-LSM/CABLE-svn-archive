MODULE veg_params_type_mod_cbl

IMPLICIT NONE

!Public subroutines
PUBLIC :: alloc_veg_params_type_cbl
PUBLIC :: dealloc_veg_params_type_cbl
PUBLIC :: assoc_veg_params_type_cbl
PUBLIC :: nullify_veg_params_type_cbl
PUBLIC :: init_veg_cbl
!Public data 
PUBLIC :: veg_params_data_type_cbl
PUBLIC :: veg_params_type_cbl

PRIVATE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='VEG_PARAMS_CBL_VARS_MOD'

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Vegetation parameters used:
TYPE veg_params_data_type_cbl

  !jhan: surely %meth doesn't need to be 'mp' here
  INTEGER, ALLOCATABLE, PUBLIC ::                                             &
    meth(:),    & ! method for calculation of canopy fluxes and temp.
    iveg(:),    & ! vegetation type
    iLU(:)        ! land use type

  REAL, ALLOCATABLE, PUBLIC ::                                                &
    canst1(:),  & ! max intercepted water by canopy (mm/LAI)
    dleaf(:),   & ! chararacteristc legnth of leaf (m)
    ejmax(:),   & ! max pot. electron transp rate top leaf(mol/m2/s)
    frac4(:),   & ! fraction of c4 plants
    hc(:),      & ! roughness height of canopy (veg - snow)
    vlai(:),    & ! leaf area index
    xalbnir(:), & !
    rp20(:),    & ! plant respiration coefficient at 20 C
    rpcoef(:),  & ! temperature coef nonleaf plant respiration (1/C)
    rs20(:),    & ! soil respiration at 20 C [mol m-2 s-1]
    shelrb(:),  & ! sheltering factor (dimensionless)
    vegcf(:),   & ! kdcorbin, 08/10
    tminvj(:),  & ! min temperature of the start of photosynthesis
    toptvj(:),  & ! opt temperature of the start of photosynthesis
    tmaxvj(:),  & ! max temperature of the start of photosynthesis
    vbeta(:),   & !
    vcmax(:),   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
    xfang(:),   & ! leaf angle PARAMETER
    extkn(:),   & ! extinction coef for vertical
    vlaimax(:), & ! extinction coef for vertical
    wai(:),     & ! wood area index (stem+branches+twigs)
    a1gs(:),    & ! a1 parameter in stomatal conductance model
    d0gs(:),    & ! d0 in stomatal conductance model
    alpha(:),   & ! initial slope of J-Q response curve
    convex(:),  & ! convexity of J-Q response curve
    cfrd(:),    & ! ratio of day respiration to vcmax
    gswmin(:),  & ! minimal stomatal conductance
    conkc0(:),  & ! Michaelis-menton constant for carboxylase
    conko0(:),  & ! Michaelis-menton constant for oxygenase
    ekc(:),     & ! activation energy for caroxylagse
    eko(:),     & ! acvtivation enegery for oxygenase
    g0(:),      & ! Belinda's stomatal model intercept, Ticket #56.
    g1(:)         ! Belinda's stomatal model slope, Ticket #56.

  REAL, ALLOCATABLE, PUBLIC ::                                                &
    refl(:,:),    & !
    taul(:,:),    & !
    froot(:,:)      ! fraction of root in each soil layer

  ! Additional  veg parameters:
  REAL, ALLOCATABLE, PUBLIC ::                                                &
    rootbeta(:),  & ! param estimating vertical root mass distribution (froot)
    gamma(:),     & ! parameter in root efficiency function (Lai and Katul 2000)
    ZR(:),        & ! maximum rooting depth (cm)
    F10(:),       & ! fraction of roots in top 10 cm
    clitt(:)        !

  REAL, ALLOCATABLE, PUBLIC ::                                                &
    soil(:,:),                                                                &
    csoil(:,:),                                                               &
    cplant(:,:),                                                              &
    ratecs(:,:),                                                              &
    ratecp(:,:)

  ! Additional POP veg param
  INTEGER, ALLOCATABLE, PUBLIC :: disturbance_interval(:,:)
  REAL, ALLOCATABLE,    PUBLIC :: disturbance_intensity(:,:)

END TYPE veg_params_data_type_cbl

! Pointer version of Vegetation parameters used:
TYPE veg_params_type_cbl

  !jhan: surely %meth doesn't need to be 'mp' here
  INTEGER, POINTER, PUBLIC ::                                             &
    meth(:),    & ! method for calculation of canopy fluxes and temp.
    iveg(:),    & ! vegetation type
    iLU(:)        ! land use type

  REAL, POINTER, PUBLIC ::                                                &
    canst1(:),  & ! max intercepted water by canopy (mm/LAI)
    dleaf(:),   & ! chararacteristc legnth of leaf (m)
    ejmax(:),   & ! max pot. electron transp rate top leaf(mol/m2/s)
    frac4(:),   & ! fraction of c4 plants
    hc(:),      & ! roughness height of canopy (veg - snow)
    vlai(:),    & ! leaf area index
    xalbnir(:), & !
    rp20(:),    & ! plant respiration coefficient at 20 C
    rpcoef(:),  & ! temperature coef nonleaf plant respiration (1/C)
    rs20(:),    & ! soil respiration at 20 C [mol m-2 s-1]
    shelrb(:),  & ! sheltering factor (dimensionless)
    vegcf(:),   & ! kdcorbin, 08/10
    tminvj(:),  & ! min temperature of the start of photosynthesis
    toptvj(:),  & ! opt temperature of the start of photosynthesis
    tmaxvj(:),  & ! max temperature of the start of photosynthesis
    vbeta(:),   & !
    vcmax(:),   & ! max RuBP carboxylation rate top leaf (mol/m2/s)
    xfang(:),   & ! leaf angle PARAMETER
    extkn(:),   & ! extinction coef for vertical
    vlaimax(:), & ! extinction coef for vertical
    wai(:),     & ! wood area index (stem+branches+twigs)
    a1gs(:),    & ! a1 parameter in stomatal conductance model
    d0gs(:),    & ! d0 in stomatal conductance model
    alpha(:),   & ! initial slope of J-Q response curve
    convex(:),  & ! convexity of J-Q response curve
    cfrd(:),    & ! ratio of day respiration to vcmax
    gswmin(:),  & ! minimal stomatal conductance
    conkc0(:),  & ! Michaelis-menton constant for carboxylase
    conko0(:),  & ! Michaelis-menton constant for oxygenase
    ekc(:),     & ! activation energy for caroxylagse
    eko(:),     & ! acvtivation enegery for oxygenase
    g0(:),      & ! Belinda's stomatal model intercept, Ticket #56.
    g1(:)         ! Belinda's stomatal model slope, Ticket #56.

  REAL, POINTER, PUBLIC ::                                                &
    refl(:,:),    & !
    taul(:,:),    & !
    froot(:,:)      ! fraction of root in each soil layer

  ! Additional  veg parameters:
  REAL, POINTER, PUBLIC ::                                                &
    rootbeta(:),  & ! param estimating vertical root mass distribution (froot)
    gamma(:),     & ! parameter in root efficiency function (Lai and Katul 2000)
    ZR(:),        & ! maximum rooting depth (cm)
    F10(:),       & ! fraction of roots in top 10 cm
    clitt(:)        !

  REAL, POINTER, PUBLIC ::                                                &
    soil(:,:),                                                                &
    csoil(:,:),                                                               &
    cplant(:,:),                                                              &
    ratecs(:,:),                                                              &
    ratecp(:,:)

  ! Additional POP veg param
  INTEGER, POINTER, PUBLIC :: disturbance_interval(:,:)
  REAL, POINTER,    PUBLIC :: disturbance_intensity(:,:)

END TYPE veg_params_type_cbl

CONTAINS

SUBROUTINE alloc_veg_params_type_cbl( mp, var )
!-----------------------------------------------------------------------------
! Description:
!   Allocates variable arrays for vegetation parameters.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

USE ereport_mod,              ONLY: ereport

USE jules_model_environment_mod,  ONLY: lsm_id, cable
USE grid_constants_mod_cbl, ONLY:                                             & 
                            nrb,  & !total # rad. "bands"
                            nscs, & ! number of soil carbon stores
                            nvcs, & ! number of vegetation carbon stores
                            nsl     ! # soil layers 

IMPLICIT NONE

INTEGER, INTENT(IN) :: mp
TYPE(veg_params_data_type_cbl), INTENT(INOUT) :: var

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_VEG_PARAMETER_TYPE_CBL'

!End of header

IF ( lsm_id == cable ) THEN
  ALLOCATE( var% iveg(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% hc(mp),      stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vlai(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% xfang(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var%refl(mp,nrb), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var%taul(mp,nrb), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% canst1(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% dleaf(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% ejmax(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% meth(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% frac4(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% xalbnir(mp), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rp20(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rpcoef(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rs20(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% shelrb(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vegcf(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% tminvj(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% tmaxvj(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vbeta(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vcmax(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% extkn(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% wai(mp),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% a1gs(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% d0gs(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% alpha(mp),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% convex(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% cfrd(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% gswmin(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% conkc0(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% conko0(mp),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% ekc(mp),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% eko(mp),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% g0(mp),      stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% g1(mp),      stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%froot(mp,nsl), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%rootbeta(mp),  stat = error  )
  error_sum = error_sum + error
  ALLOCATE( var % GAMMA(mp),  stat = error  )
  error_sum = error_sum + error
  ALLOCATE( var % f10(mp),    stat = error  )
  error_sum = error_sum + error
  ALLOCATE ( var % zr(mp),    stat = error )
  error_sum = error_sum + error
  ALLOCATE ( var % clitt(mp), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%csoil(mp,nscs),  stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%cplant(mp,nvcs), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%ratecp(mp,nvcs), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%ratecs(mp,nscs), stat = error )
  error_sum = error_sum + error
  !H!    ALLOCATE( var% iLU(mp) )
  !H!    ALLOCATE( var% toptvj(mp) )
  !H!    ALLOCATE( var%deciduous(mp) )
  !H!    ALLOCATE( var%vlaimax(mp) )
  !H!    ALLOCATE ( var % disturbance_interval(mp,2) )
  !H!    ALLOCATE ( var % disturbance_intensity(mp,2) )
ELSE ! lsm_id==JULES
  ALLOCATE( var% iveg(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% hc(1),      stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vlai(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% xfang(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var%refl(1,1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var%taul(1,1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% canst1(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% dleaf(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% ejmax(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% meth(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% frac4(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% xalbnir(1), stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rp20(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rpcoef(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% rs20(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% shelrb(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vegcf(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% tminvj(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% tmaxvj(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vbeta(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% vcmax(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% extkn(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% wai(1),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% a1gs(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% d0gs(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% alpha(1),   stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% convex(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% cfrd(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% gswmin(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% conkc0(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% conko0(1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% ekc(1),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% eko(1),     stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% g0(1),      stat = error )
  error_sum = error_sum + error
  ALLOCATE( var% g1(1),      stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%froot(1,1), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%rootbeta(1),  stat = error  )
  error_sum = error_sum + error
  ALLOCATE( var % GAMMA(1),  stat = error  )
  error_sum = error_sum + error
  ALLOCATE( var % f10(1),    stat = error  )
  error_sum = error_sum + error
  ALLOCATE ( var % zr(1),    stat = error )
  error_sum = error_sum + error
  ALLOCATE ( var % clitt(1), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%csoil(1,1),  stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%cplant(1,1), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%ratecp(1,1), stat = error )
  error_sum = error_sum + error
  ALLOCATE(var%ratecs(1,1), stat = error )
  error_sum = error_sum + error
ENDIF

!init method to zero
IF (error_sum == 0) THEN
  var%meth(:)   = 0
ELSE
  CALL ereport("alloc_veg_cbl_params_type", errcode,                            &
                "Error allocating CABLE model veg parameter arrays")
END IF

RETURN

END SUBROUTINE alloc_veg_params_type_cbl

SUBROUTINE dealloc_veg_params_type_cbl( var )
  !Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(veg_params_data_type_cbl), INTENT(INOUT) :: var
  
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='dealloc_veg_params_type_cbl'
  
!End of header
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE(var%iveg   )
DEALLOCATE(var%hc     )
DEALLOCATE(var%vlai   )
DEALLOCATE(var%xfang  )
DEALLOCATE(var%refl   )
DEALLOCATE(var%taul   )
DEALLOCATE(var%canst1 )
DEALLOCATE(var%dleaf  )
DEALLOCATE(var%ejmax  )
DEALLOCATE(var%meth   )
DEALLOCATE(var%frac4  )
DEALLOCATE(var%xalbnir)
DEALLOCATE(var%rp20   )
DEALLOCATE(var%rpcoef )
DEALLOCATE(var%rs20   )
DEALLOCATE(var%shelrb )
DEALLOCATE(var%vegcf  )
DEALLOCATE(var%tminvj )
DEALLOCATE(var%tmaxvj )
DEALLOCATE(var%vbeta  )
DEALLOCATE(var%vcmax  )
DEALLOCATE(var%extkn  )
DEALLOCATE(var%wai    )
DEALLOCATE(var%a1gs   )
DEALLOCATE(var%d0gs   )
DEALLOCATE(var%alpha  )
DEALLOCATE(var%convex )
DEALLOCATE(var%cfrd   )
DEALLOCATE(var%gswmin )
DEALLOCATE(var%conkc0 )
DEALLOCATE(var%conko0 )
DEALLOCATE(var%ekc    )
DEALLOCATE(var%eko    )
DEALLOCATE(var%g0     )
DEALLOCATE(var%g1     )
DEALLOCATE(var%froot  )
DEALLOCATE(var%rootbeta )
DEALLOCATE(var% GAMMA  )
DEALLOCATE(var% f10    )
DEALLOCATE(var% zr     )
DEALLOCATE(var% clitt  )
DEALLOCATE(var%csoil  )
DEALLOCATE(var%cplant )
DEALLOCATE(var%ratecp )
DEALLOCATE(var%ratecs )
  
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN


END SUBROUTINE dealloc_veg_params_type_cbl


SUBROUTINE assoc_veg_params_type_cbl( var, var_ptr )

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(veg_params_data_type_cbl), INTENT(INOUT), TARGET  :: var
TYPE(veg_params_type_cbl), INTENT(INOUT) :: var_ptr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='assoc_veg_params_type_cbl'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL nullify_veg_params_type_cbl(var_ptr)

     var_ptr%iveg   =>     var%iveg   
     var_ptr%hc     =>     var%hc     
     var_ptr%vlai   =>     var%vlai   
     var_ptr%xfang  =>     var%xfang  
     var_ptr%refl   =>     var%refl   
     var_ptr%taul   =>     var%taul   
     var_ptr%canst1 =>     var%canst1 
     var_ptr%dleaf  =>     var%dleaf  
     var_ptr%ejmax  =>     var%ejmax  
     var_ptr%meth   =>     var%meth   
     var_ptr%frac4  =>     var%frac4  
     var_ptr%xalbnir=>     var%xalbnir
     var_ptr%rp20   =>     var%rp20   
     var_ptr%rpcoef =>     var%rpcoef 
     var_ptr%rs20   =>     var%rs20   
     var_ptr%shelrb =>     var%shelrb 
     var_ptr%vegcf  =>     var%vegcf  
     var_ptr%tminvj =>     var%tminvj 
     var_ptr%tmaxvj =>     var%tmaxvj 
     var_ptr%vbeta  =>     var%vbeta  
     var_ptr%vcmax  =>     var%vcmax  
     var_ptr%extkn  =>     var%extkn  
     var_ptr%wai    =>     var%wai    
     var_ptr%a1gs   =>     var%a1gs   
     var_ptr%d0gs   =>     var%d0gs   
     var_ptr%alpha  =>     var%alpha  
     var_ptr%convex =>     var%convex 
     var_ptr%cfrd   =>     var%cfrd   
     var_ptr%gswmin =>     var%gswmin 
     var_ptr%conkc0 =>     var%conkc0 
     var_ptr%conko0 =>     var%conko0 
     var_ptr%ekc    =>     var%ekc    
     var_ptr%eko    =>     var%eko    
     var_ptr%g0     =>     var%g0     
     var_ptr%g1     =>     var%g1     
     var_ptr%froot  =>     var%froot  
     var_ptr%rootbeta =>     var%rootbeta
     var_ptr% GAMMA =>     var% GAMMA 
     var_ptr% f10   =>     var% f10   
     var_ptr% zr    =>     var% zr    
     var_ptr% clitt =>     var% clitt 
     var_ptr%csoil  =>     var%csoil  
     var_ptr%cplant =>     var%cplant 
     var_ptr%ratecp =>     var%ratecp 
     var_ptr%ratecs =>     var%ratecs 
 
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE assoc_veg_params_type_cbl

SUBROUTINE nullify_veg_params_type_cbl(var_ptr)

  !No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook
  
IMPLICIT NONE
  
!Arguments
TYPE(veg_params_type_cbl), INTENT(INOUT) :: var_ptr
  
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
  
CHARACTER(LEN=*), PARAMETER :: RoutineName='nullify_veg_params_type_cbl'
  
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY(    var_ptr%iveg   )
NULLIFY(    var_ptr%hc     )
NULLIFY(    var_ptr%vlai   )
NULLIFY(    var_ptr%xfang  ) 
NULLIFY(    var_ptr%refl   ) 
NULLIFY(    var_ptr%taul   )
NULLIFY(    var_ptr%canst1 )
NULLIFY(    var_ptr%dleaf  )
NULLIFY(    var_ptr%ejmax  )
NULLIFY(    var_ptr%meth   )
NULLIFY(    var_ptr%frac4  )
NULLIFY(    var_ptr%xalbnir)
NULLIFY(    var_ptr%rp20   )
NULLIFY(    var_ptr%rpcoef )
NULLIFY(    var_ptr%rs20   )
NULLIFY(    var_ptr%shelrb )
NULLIFY(    var_ptr%vegcf  )
NULLIFY(    var_ptr%tminvj )
NULLIFY(    var_ptr%tmaxvj )
NULLIFY(    var_ptr%vbeta  )
NULLIFY(    var_ptr%vcmax  )
NULLIFY(    var_ptr%extkn  )
NULLIFY(    var_ptr%wai    )
NULLIFY(    var_ptr%a1gs   )
NULLIFY(    var_ptr%d0gs   )
NULLIFY(    var_ptr%alpha  )
NULLIFY(    var_ptr%convex )
NULLIFY(    var_ptr%cfrd   )
NULLIFY(    var_ptr%gswmin )
NULLIFY(    var_ptr%conkc0 )
NULLIFY(    var_ptr%conko0 )
NULLIFY(    var_ptr%ekc    )
NULLIFY(    var_ptr%eko    )
NULLIFY(    var_ptr%g0     )
NULLIFY(    var_ptr%g1     )
NULLIFY(    var_ptr%froot  )
NULLIFY(    var_ptr%rootbeta )
NULLIFY(    var_ptr% GAMMA )
NULLIFY(    var_ptr% f10   )
NULLIFY(    var_ptr% zr    )
NULLIFY(    var_ptr% clitt )
NULLIFY(    var_ptr%csoil  )
NULLIFY(    var_ptr%cplant )
NULLIFY(    var_ptr%ratecp )
NULLIFY(    var_ptr%ratecs )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nullify_veg_params_type_cbl
 
SUBROUTINE init_veg_cbl( mp, frac_surft, veg, vegin, L_tile_pts )

USE vegin_pars_mod_cbl,  ONLY: vegin_type
USE ancil_info,           ONLY: nsurft, land_pts

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   init_cables veg parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER ::mp
TYPE(veg_params_data_type_cbl), TARGET :: veg
TYPE(vegin_type),               TARGET :: vegin
LOGICAL :: L_tile_pts(land_pts,nsurft)
INTEGER :: JSurfaceTypeID(land_pts,nsurft)  
REAL:: frac_surft(land_pts,nsurft)  
INTEGER :: i
INTEGER :: h

!local var to pack surface type:
JSurfaceTypeID = 0
DO i = 1,nsurft
  IF ( frac_surft(1,i) > 0 ) JSurfaceTypeID(:,i) = i
END DO

veg%iveg = PACK( JSurfaceTypeID, L_tile_pts)

! Prescribe parameters for current gridcell based on veg/soil type
! (which may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid
  veg%taul(h,1)   = vegin%taul(1,veg%iveg(h))
  veg%taul(h,2)   = vegin%taul(2,veg%iveg(h))
  veg%refl(h,1)   = vegin%refl(1,veg%iveg(h))
  veg%refl(h,2)   = vegin%refl(2,veg%iveg(h))
  veg%cplant(h,1)   = vegin%cplant(1,veg%iveg(h))
  veg%cplant(h,2)   = vegin%cplant(2,veg%iveg(h))
  veg%cplant(h,3)   = vegin%cplant(3,veg%iveg(h))
  veg%csoil(h,1)   = vegin%csoil(1,veg%iveg(h))
  veg%csoil(h,2)   = vegin%csoil(2,veg%iveg(h))
  veg%ratecp(h,1)   = vegin%ratecp(1,veg%iveg(h))
  veg%ratecp(h,2)   = vegin%ratecp(2,veg%iveg(h))
  veg%ratecp(h,3)   = vegin%ratecp(3,veg%iveg(h))
  veg%ratecs(h,1)   = vegin%ratecs(1,veg%iveg(h))
  veg%ratecs(h,2)   = vegin%ratecs(2,veg%iveg(h))
  veg%hc(h)       = vegin%hc(veg%iveg(h))
  veg%xfang(h)    = vegin%xfang(veg%iveg(h))
  veg%frac4(h)    = vegin%frac4(veg%iveg(h))
  veg%canst1(h)   = vegin%canst1(veg%iveg(h))
  veg%dleaf(h)    = vegin%dleaf(veg%iveg(h))
  veg%vcmax(h)    = vegin%vcmax(veg%iveg(h))
  veg%ejmax(h)    = vegin%ejmax(veg%iveg(h))
  veg%vbeta(h)    = vegin%vbeta(veg%iveg(h))
  veg%xalbnir(h)  = vegin%xalbnir(veg%iveg(h))
  veg%rp20(h)     = vegin%rp20(veg%iveg(h))
  veg%rpcoef(h)   = vegin%rpcoef(veg%iveg(h))
  veg%rs20(h)     = vegin%rs20(veg%iveg(h))
  veg%shelrb(h)   = vegin%shelrb(veg%iveg(h))
  veg%wai(h)      = vegin%wai(veg%iveg(h))
  veg%a1gs(h)     = vegin%a1gs(veg%iveg(h))
  veg%d0gs(h)     = vegin%d0gs(veg%iveg(h))
  veg%vegcf(h)    = vegin%vegcf(veg%iveg(h))
  veg%extkn(h)    = vegin%extkn(veg%iveg(h))
  veg%tminvj(h)   = vegin%tminvj(veg%iveg(h))
  veg%tmaxvj(h)   = vegin%tmaxvj(veg%iveg(h))
  veg%g0(h)       = vegin%g0(veg%iveg(h)) ! Ticket #56
  veg%g1(h)       = vegin%g1(veg%iveg(h)) ! Ticket #56
  veg%a1gs(h)   = vegin%a1gs(veg%iveg(h))
  veg%d0gs(h)   = vegin%d0gs(veg%iveg(h))
  veg%alpha(h)  = vegin%alpha(veg%iveg(h))
  veg%convex(h) = vegin%convex(veg%iveg(h))
  veg%cfrd(h)   = vegin%cfrd(veg%iveg(h))
  veg%gswmin(h) = vegin%gswmin(veg%iveg(h))
  veg%conkc0(h) = vegin%conkc0(veg%iveg(h))
  veg%conko0(h) = vegin%conko0(veg%iveg(h))
  veg%ekc(h)    = vegin%ekc(veg%iveg(h))
  veg%eko(h)    = vegin%eko(veg%iveg(h))
  veg%rootbeta(h)  = vegin%rootbeta(veg%iveg(h))
  veg%zr(h)       = vegin%zr(veg%iveg(h))
  veg%clitt(h)    = vegin%clitt(veg%iveg(h))
END DO ! over each veg patch in land point

END SUBROUTINE init_veg_cbl

END MODULE veg_params_type_mod_cbl
