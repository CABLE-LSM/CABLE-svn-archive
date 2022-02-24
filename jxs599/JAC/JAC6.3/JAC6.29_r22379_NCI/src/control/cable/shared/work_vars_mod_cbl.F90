MODULE work_vars_mod_cbl

!------------------------------------------------------------------------------
! Description:
!
!   Declares/(de)allocates/assigns CABLE "working" variables
!   These are vars req'd to be kept across the surf_couple* pathways
!   and/or timesteps. Some will be elevated to rognostics, others will be
!   removed via rewriting of the algorithm where req'd
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC :: alloc_work_vars_cbl
PUBLIC :: assoc_work_vars_cbl
PUBLIC :: work_vars_data_type
PUBLIC :: work_vars_type
PRIVATE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WORK_VARS_MOD_CBL'

TYPE :: work_vars_data_type
  
  !fields returned @ hydrology (surf_couple_extra) level of JULES
  !computed in CABLE @ implicit (surf_couple_implicit) level 
  REAL, ALLOCATABLE, PUBLIC :: SNOW_TILE(:,:) 
  REAL, ALLOCATABLE, PUBLIC :: LYING_SNOW(:) 
  REAL, ALLOCATABLE, PUBLIC :: SURF_ROFF(:) 
  REAL, ALLOCATABLE, PUBLIC :: SUB_SURF_ROFF(:) 
  REAL, ALLOCATABLE, PUBLIC :: TOT_TFALL(:) 

END TYPE work_vars_data_type

TYPE :: work_vars_type

  !fields returned @ hydrology (surf_couple_extra) level of JULES
  !computed in CABLE @ implicit (surf_couple_implicit) level 
  REAL, POINTER, PUBLIC :: SNOW_TILE(:,:)
  REAL, POINTER, PUBLIC :: LYING_SNOW(:)
  REAL, POINTER, PUBLIC :: SURF_ROFF(:) 
  REAL, POINTER, PUBLIC :: SUB_SURF_ROFF(:) 
  REAL, POINTER, PUBLIC :: TOT_TFALL(:) 

END TYPE work_vars_type

CONTAINS

!===============================================================================
SUBROUTINE alloc_work_vars_cbl(land_pts, nsurft, sm_levels, lsm_id, cable,     &
                                work_data_cbl )

!Replacements for the argument list
USE grid_constants_mod_cbl,       ONLY: nsnl

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook
USE jules_print_mgr,          ONLY: jules_message, jules_print, PrNorm
USE ereport_mod,              ONLY: ereport

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: land_pts, nsurft,sm_levels
INTEGER, INTENT(IN) :: lsm_id, cable

TYPE(work_vars_data_type), INTENT(IN OUT) :: work_data_cbl

!-----------------------------------------------------------------------
! Local variables for error trapping
!-----------------------------------------------------------------------
INTEGER ::                                                                     &
  ERROR        = 0,                                                            &
                       ! Variable for trapping the error from each
                       ! individual call to allocate
  error_sum    = 0,                                                            &

                       ! Variable to track the sum of all errors
                       ! resulting from calls to allocate. Hence we
                       ! know that everything was successful if and
                       ! only if this is zero at the end
  errcode
                       ! Variable to use in error report

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! CABLE vars to be initialized
IF ( lsm_id == cable ) THEN

  ALLOCATE( work_data_cbl%SNOW_TILE(land_pts, nsurft), STAT = ERROR )
  ALLOCATE( work_data_cbl%LYING_SNOW(land_pts),        STAT = ERROR )
  ALLOCATE( work_data_cbl%SURF_ROFF(land_pts),         STAT = ERROR )
  ALLOCATE( work_data_cbl%TOT_TFALL(land_pts),         STAT = ERROR )
  ALLOCATE( work_data_cbl%SUB_SURF_ROFF(land_pts),     STAT = ERROR )
  error_sum = error_sum + ERROR

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  IF ( error_sum /= 0 )                                                        &
  CALL ereport(RoutineName, errcode,                                           &
              "Error allocating CABLE work type")

ELSE

  ALLOCATE( work_data_cbl%SNOW_TILE(1,1), STAT = ERROR )
  ALLOCATE( work_data_cbl%LYING_SNOW(1),  STAT = ERROR )
  ALLOCATE( work_data_cbl%SURF_ROFF(1),   STAT = ERROR )
  ALLOCATE( work_data_cbl%TOT_TFALL(1),   STAT = ERROR )
  ALLOCATE( work_data_cbl%SUB_SURF_ROFF(1), STAT = ERROR )
  error_sum = error_sum + ERROR

  !-----------------------------------------------------------------------
  ! Write out an error if there was one
  !-----------------------------------------------------------------------
  IF ( error_sum /= 0 )                                                        &
  CALL ereport(RoutineName, errcode,                                           &
              "Error allocating CABLE work type")

END IF

work_data_cbl% SNOW_TILE      (:,:)       = 0.0
work_data_cbl% LYING_SNOW     (:)         = 0.0
work_data_cbl% SURF_ROFF      (:)         = 0.0
work_data_cbl% SUB_SURF_ROFF  (:)         = 0.0
work_data_cbl% TOT_TFALL      (:)         = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE alloc_work_vars_cbl

!===============================================================================
SUBROUTINE dealloc_work_vars_cbl(work_data_cbl )

!Common Non-science modules
USE parkind1,                 ONLY: jprb, jpim
USE yomhook,                  ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_data_type), INTENT(IN OUT) :: work_data_cbl

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEALLOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! CABLE vars to be initialized via JULES i/o
DEALLOCATE( work_data_cbl %SNOW_TILE     )
DEALLOCATE( work_data_cbl %LYING_SNOW    )
DEALLOCATE( work_data_cbl %SURF_ROFF     )
DEALLOCATE( work_data_cbl %SUB_SURF_ROFF )
DEALLOCATE( work_data_cbl %TOT_TFALL     )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dealloc_work_vars_cbl

!===============================================================================
SUBROUTINE assoc_work_vars_cbl(work_cbl, work_data_cbl )

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_type),      INTENT(IN OUT)         :: work_cbl
TYPE(work_vars_data_type), INTENT(IN OUT), TARGET :: work_data_cbl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASSOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL nullify_assoc_work_vars_cbl(work_cbl)

work_cbl% SNOW_TILE     => work_data_cbl% SNOW_TILE    
work_cbl% LYING_SNOW    => work_data_cbl% LYING_SNOW   
work_cbl% SURF_ROFF     => work_data_cbl% SURF_ROFF    
work_cbl% SUB_SURF_ROFF => work_data_cbl% SUB_SURF_ROFF
work_cbl% TOT_TFALL     => work_data_cbl% TOT_TFALL    

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE assoc_work_vars_cbl

!===============================================================================
SUBROUTINE nullify_assoc_work_vars_cbl(work_cbl)

  !No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
TYPE(work_vars_type), INTENT(IN OUT) :: work_cbl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NULLIFY_ASSOC_WORK_VARS_CBL'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

NULLIFY( work_cbl% SNOW_TILE      )
NULLIFY( work_cbl% LYING_SNOW     )
NULLIFY( work_cbl% SURF_ROFF      )
NULLIFY( work_cbl% SUB_SURF_ROFF  )
NULLIFY( work_cbl% TOT_TFALL      )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE nullify_assoc_work_vars_cbl

END MODULE work_vars_mod_cbl

