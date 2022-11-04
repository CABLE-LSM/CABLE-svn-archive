MODULE popdriver_casa_mod

IMPLICIT NONE

#define  UM_BUILD YES

#ifndef UM_BUILD

CONTAINS

SUBROUTINE POPdriver(casaflux,casabal,veg, POP)

  USE cable_def_types_mod
  USE cable_common_module, only: cable_runtime
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module,  ONLY: CurYear, CABLE_USER
  USE TypeDef,              ONLY: i4b, dp
  USE POPMODULE,            ONLY: POPStep
  USE POP_TYPES,            ONLY: POP_TYPE


  IMPLICIT NONE


  TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters
  TYPE (casa_flux),           INTENT(IN) :: casaflux
  TYPE (casa_balance),        INTENT(IN) :: casabal
  TYPE(POP_TYPE),             INTENT(INOUT) :: POP

  INTEGER                                   :: it, nit
  REAL(dp)                               :: StemNPP(mp,2)
  REAL(dp), allocatable :: NPPtoGPP(:)
  REAL(dp), allocatable ::  LAImax(:)  , Cleafmean(:),  Crootmean(:)
  CHARACTER                                 :: cyear*4
  CHARACTER                                 :: ncfile*99
  !! vh_js !!
  INTEGER, allocatable :: Iw(:) ! array of indices corresponding to woody (shrub or forest) tiles

  ! INTEGER, INTENT(IN) :: wlogn
  INTEGER , parameter :: wlogn=6

  if (.NOT.Allocated(LAIMax)) allocate(LAIMax(mp))
  if (.NOT.Allocated(Cleafmean))  allocate(Cleafmean(mp))
  if (.NOT.Allocated(Crootmean)) allocate(Crootmean(mp))
  if (.NOT.Allocated(NPPtoGPP)) allocate(NPPtoGPP(mp))
  if (.NOT.Allocated(Iw)) allocate(Iw(POP%np))

  IF (cable_user%CALL_POP .and. POP%np.gt.0) THEN ! CALL_POP
     Iw = POP%Iwood

     StemNPP(:,1) = casaflux%stemnpp
     StemNPP(:,2) = 0.0
     WHERE (casabal%FCgppyear > 1.e-5 .and. casabal%FCnppyear > 1.e-5  )
        NPPtoGPP = casabal%FCnppyear/casabal%FCgppyear
     ELSEWHERE
        NPPtoGPP = 0.5
     ENDWHERE
     LAImax = casabal%LAImax
     Cleafmean = casabal%cleafmean
     Crootmean = casabal%Crootmean

     CALL POPStep(pop, max(StemNPP(Iw,:)/1000.,0.0001), int(veg%disturbance_interval(Iw,:), i4b),&
          real(veg%disturbance_intensity(Iw,:),dp)      ,&
          max(LAImax(Iw),0.001), Cleafmean(Iw), Crootmean(Iw), NPPtoGPP(Iw))


  ENDIF ! CALL_POP

END SUBROUTINE POPdriver

#endif

END MODULE popdriver_casa_mod
