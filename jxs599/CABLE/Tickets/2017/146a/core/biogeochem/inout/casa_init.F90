SUBROUTINE casa_init(casabiome,casamet,casaflux,casapool,casabal,veg,phen)
! mst not used (BP sep2010)
!! for first time reading file *_1220.csv  (BP may2010)
!SUBROUTINE casa_init(mst,casapool,casabal,veg)
!!SUBROUTINE casa_init(mst,casapool,casabal)
!! end addition (BP may2010)
!  initialize some values in phenology parameters and leaf growth phase
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
! for first time reading file *_1220.csv  (BP may2010)
  USE cable_def_types_mod
  USE cable_io_vars_module, ONLY: landpt, patch
  USE cable_common_module, only: cable_user

! end addition (BP may2010)
  IMPLICIT NONE
!  INTEGER,        INTENT(IN)    :: mst
  TYPE (casa_biome),   INTENT(IN)    :: casabiome
  TYPE (casa_met),     INTENT(INOUT) :: casamet
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
! for first time reading file *_1220.csv  (BP may2010)
  TYPE (veg_parameter_type), INTENT(IN) :: veg
  TYPE (phen_variable),   INTENT(INOUT) :: phen
  REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
  REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
  REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc
! end addition (BP may2010)

  ! local variables
  INTEGER   :: np,npt,npz
  INTEGER   :: nyearz,ivtz,istz,isoz
  REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz
  LOGICAL   :: EXRST

  
if (.NOT.cable_user%casa_fromzero) THEN
   PRINT *, 'initial pool from restart file'
ENDIF
  PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp
  !phen%phase = 2

  !CLN initialise all !!!!! THIS NEEDS FIXING because of e.g. ICE-WATER
  casaflux%Cgpp         = 0.
  casaflux%Cnpp         = 0.
  casaflux%Crp          = 0.
  casaflux%Crgplant     = 0.
 ! casaflux%Nminfix      = 0.
  casaflux%Nminuptake   = 0.
  casaflux%Plabuptake   = 0.
  casaflux%Clabloss     = 0.
  casaflux%fracClabile  = 0.
  casaflux%stemnpp      = 0.
  casaflux%frac_sapwood = 0.
  casaflux%sapwood_area = 0.
  casaflux%FluxCtohwp = 0.
  casaflux%FluxCtoClear = 0.
  casaflux%fracCalloc   = 0.
  casaflux%fracNalloc   = 0.
  casaflux%fracPalloc   = 0.
  casaflux%Crmplant     = 0.
  casaflux%kplant       = 0.

  casaflux%fromPtoL     = 0.

  casaflux%Cnep         = 0.
  casaflux%Crsoil       = 0.
  casapool%dClabiledt = 0.0
  !casaflux%Nmindep      =  casaflux%Nmindep /2.0
 !casaflux%Nmindep      = 0.
  casaflux%Nminloss     = 0.
  casaflux%Nminleach    = 0.
  casaflux%Nupland      = 0.
  casaflux%Nlittermin   = 0.
  casaflux%Nsmin        = 0.
  casaflux%Nsimm        = 0.
  casaflux%Nsnet        = 0.
  !casaflux%fNminloss    = 0.
  !casaflux%fNminleach   = 0.
  !casaflux%Pdep         = 0.
  !casaflux%Pwea         = 0.
  casaflux%Pleach       = 0.
  casaflux%Ploss        = 0.
  casaflux%Pupland      = 0.
  casaflux%Plittermin   = 0.
  casaflux%Psmin        = 0.
  casaflux%Psimm        = 0.
  casaflux%Psnet        = 0.
!  casaflux%fPleach      = 0. !vh ! this should be a parameter, not a flux variable
  casaflux%kplab        = 0.
  casaflux%kpsorb       = 0.
  casaflux%kpocc        = 0.
  casaflux%kmlabp       = 0.
!  casaflux%Psorbmax     = 0. !vh ! this should be a paramter, not a flux variable

  casaflux%klitter      = 0.
  casaflux%ksoil        = 0.
  casaflux%fromLtoS     = 0.
  casaflux%fromStoS     = 0.
  casaflux%fromLtoCO2   = 0.
  casaflux%fromStoCO2   = 0.
  casaflux%FluxCtolitter= 0.
  casaflux%FluxNtolitter= 0.
  casaflux%FluxPtolitter= 0.
  casaflux%FluxCtosoil  = 0.
  casaflux%FluxNtosoil  = 0.
  casaflux%FluxPtosoil  = 0.
  casaflux%FluxCtoCO2   = 0.

  casaflux%Cplant_turnover = 0.

  phen%doyphase(:,1) = -50
  phen%doyphase(:,2) = phen%doyphase(:,1) +14
  phen%doyphase(:,3) = 367
  phen%doyphase(:,4) = phen%doyphase(:,3) + 14
  phen%phase(:) = 2
  phen%phen(:) = 1
  phen%aphen(:) = 0
  !CLN add more if necessary

  IF (initcasa==1) THEN
     if (.NOT.cable_user%casa_fromzero) THEN
#ifndef UM_BUILD
        CALL READ_CASA_RESTART_NC (  casamet, casapool, casaflux, phen )
#endif
     ELSE
        WRITE(*,*)'casa_init: not using restart file!'
        WRITE(*,*)'Using input from readbiome.!!!'
        WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
        casaflux%frac_sapwood(:) = 1.0
        casaflux%sapwood_area(:) = 0.0
     ENDIF
  ENDIF
 WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
 WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
 WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
!!$IF (initcasa==1) THEN
!!$     INQUIRE( FILE=TRIM(casafile%cnpipool), EXIST=EXRST )
!!$!! vh_js!!
!!$     IF ( EXRST ) THEN
!!$
!!$           PRINT*, ' Reading cnppoolOutfile as input: ,',casafile%cnpipool
!!$
!!$    OPEN(99,file=casafile%cnpipool)
!!$
!!$    DO npt =1, mp
!!$       SELECT CASE(icycle)
!!$       CASE(1)
!!$          !! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt) , &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt) ,casapool%cplant(npt,:) ,  &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:), &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt)
!!$
!!$
!!$             ELSE
!!$              READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt) , &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt) ,casapool%cplant(npt,:) ,  &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$          ENDIF
!!$
!!$
!!$       CASE(2)
!!$!! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!!$
!!$          ELSE
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$
!!$          ENDIF
!!$       CASE(3)
!!$!! vh_js !!
!!$          IF (cable_user%CALL_POP) THEN
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
!!$                  casapool%pplant(npt,:),casapool%plitter(npt,:),      &
!!$                  casapool%psoil(npt,:),casapool%psoillab(npt),        &
!!$                  casapool%psoilsorb(npt),casapool%psoilocc(npt)
!!$          ELSE
!!$             READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
!!$                  casamet%glai(npt),slaz,phen%phase(npt), &
!!$                  phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
!!$                  casapool%clabile(npt),casapool%cplant(npt,:),   &
!!$                  casapool%clitter(npt,:),casapool%csoil(npt,:),       &
!!$                  casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
!!$                  casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
!!$                  casapool%pplant(npt,:),casapool%plitter(npt,:),      &
!!$                  casapool%psoil(npt,:),casapool%psoillab(npt),        &
!!$                  casapool%psoilsorb(npt),casapool%psoilocc(npt)
!!$             casaflux%frac_sapwood(:) = 1.0
!!$             casaflux%sapwood_area(:) = 0.0
!!$
!!$
!!$          ENDIF
!!$       END SELECT
!!$       IF (ABS(patch(npt)%longitude - lonz) > 0.9 .OR. &
!!$            ABS(patch(npt)%latitude  - latz) > 0.9) THEN
!!$          PRINT *, 'patch(npt)%longitude, lonz:', patch(npt)%longitude, lonz
!!$          PRINT *, 'patch(npt)%latitude,  latz:', patch(npt)%latitude,  latz
!!$          PRINT *, 'npt = ', npt
!!$          STOP
!!$       ENDIF
!!$    ENDDO
!!$    CLOSE(99)
!!$
!!$
!!$ ELSE
!!$ !! vh_js !!
!!$    WRITE(*,*)'No valid restart file for casa_init found.'
!!$    WRITE(*,*)'Using input from readbiome.!!!'
!!$    WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
!!$    casaflux%frac_sapwood(:) = 1.0
!!$    casaflux%sapwood_area(:) = 0.0
!!$
!!$
!!$ ENDIF  ! IF (EXRST)

!!$ENDIF
!92 format(5(i6,2x),5(f18.6,3x),2(i6,',',2x),',',2x,100(f18.6,3x))
92    format(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),',',2x,100(f18.6,',',2x))


if(initcasa==0) then
   nyearz = 1
   do npt=1,mp
      casamet%lon(npt) = patch(npt)%longitude
      casamet%lat(npt) = patch(npt)%latitude
   enddo
endif

  ! reset labile C pool,comment out by Q.Zhang 10/09/2011
  !  casapool%clabile    = 0.0
  ! check pool sizes
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0
  casabal%FCgppyear=0.0;casabal%FCrpyear=0.0
  casabal%FCnppyear=0;casabal%FCrsyear=0.0;casabal%FCneeyear=0.0
  !vh !
  WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
  IF (icycle==1) THEN
    casapool%Nplant(:,:) = casapool%cplant(:,:) * casapool%ratioNCplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%Nsoil(:,:)/casapool%ratioNPsoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF

  IF (icycle >=1) THEN
     casapool%nplant     = MAX(1.e-6,casapool%nplant)
     casapool%nlitter    = MAX(1.e-6,casapool%nlitter)
     casapool%nsoil      = MAX(1.e-6,casapool%nsoil)
     casapool%nsoilmin   = MAX(1.e-6,casapool%nsoilmin)
     casabal%nplantlast  = casapool%nplant
     casabal%nlitterlast = casapool%nlitter
     casabal%nsoillast   = casapool%nsoil
     casabal%nsoilminlast= casapool%nsoilmin
     casabal%sumnbal     = 0.0
     casabal%FNdepyear=0.0;casabal%FNfixyear=0.0;casabal%FNsnetyear=0.0
     casabal%FNupyear=0.0;casabal%FNleachyear=0.0;casabal%FNlossyear=0.0
     !vh !
     WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
  ENDIF

  IF (icycle >=1) THEN
     casapool%pplant       = MAX(1.0e-7,casapool%pplant)
     casapool%plitter      = MAX(1.0e-7,casapool%plitter)
     casapool%psoil        = MAX(1.0e-7,casapool%psoil)
     casapool%Psoillab     = MAX(1.0e-7,casapool%psoillab)  ! was 2.0, changed according to  YP
     casapool%psoilsorb    = MAX(1.0e-7,casapool%psoilsorb) ! was 10.0, -
     casapool%psoilocc     = MAX(1.0e-7,casapool%psoilocc)  ! was 50.0, -
     casabal%pplantlast    = casapool%pplant
     casabal%plitterlast   = casapool%plitter
     casabal%psoillast     = casapool%psoil
     casabal%psoillablast  = casapool%psoillab
     casabal%psoilsorblast = casapool%psoilsorb
     casabal%psoilocclast  = casapool%psoilocc
     casabal%sumpbal       = 0.0
     casabal%FPweayear=0.0;casabal%FPdustyear=0.0; casabal%FPsnetyear=0.0
     casabal%FPupyear=0.0;casabal%FPleachyear=0.0;casabal%FPlossyear=0.0
     !vh !
     WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
  EndIF


END SUBROUTINE casa_init


