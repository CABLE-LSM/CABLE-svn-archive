!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: Input and output code for CASA-CNP when run offline
!          ACCESS version may use some of this code but split into different files?
!
! Contact: Yingping.Wang@csiro.au and Bernard.Pak@csiro.au
!
! History: Developed for offline code.  Expect to re-write for MPI and ACCESS 
!          versions
!
!
! ==============================================================================
! casa_inout.f90
!
! the following routines are used when "casacnp" is coupled to "cable"
!   casa_readbiome
!   casa_readphen
!   casa_readpoint   (removed, now done in parameter_module)
!   casa_init
!   casa_poolout
!   casa_cnpflux  (not used?)
MODULE casa_inout_module

CONTAINS

SUBROUTINE casa_readphen(veg,casamet,phen)
  ! read in the tabulated modis-derived leaf phenology data
  ! for latitude bands of 79.75 to -55.25
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module, ONLY : knode_gl
  IMPLICIT NONE
  TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
  TYPE (casa_met),           INTENT(IN)    :: casamet
  TYPE (phen_variable),      INTENT(INOUT) :: phen

  ! local variables
  INTEGER, PARAMETER            :: nphen=8! was 10(IGBP). changed by Q.Zhang @01/12/2011
  INTEGER np,nx,ilat
  INTEGER, DIMENSION(271,mvtype) :: greenup, fall,  phendoy1
  INTEGER, DIMENSION(nphen)     :: greenupx,fallx,xphendoy1
  INTEGER, DIMENSION(nphen)     :: ivtx
  REAL(r_2), DIMENSION(271)     :: xlat

  ! initilize for evergreen PFTs
  greenup(:,:) = -50
  fall(:,:)    = 367
  phendoy1(:,:)= 2

  OPEN(101,file=casafile%phen)

  if (knode_gl==0) then
    print *, '  '; print *, 'CASA_log:'
    print *, '  Opened file - '
    print *, '  ', trim(casafile%phen)
    print *, '  for reading phen vars.'
    print *, 'End CASA_log:'; print *, '  '
  endif

  READ(101,*)
  READ(101,*) (ivtx(nx),nx=1,nphen) ! fixed at 10, as only 10 of 17 IGBP PFT
                                    ! have seasonal leaf phenology
  DO ilat=271,1,-1
    READ(101,*) xlat(ilat),(greenupx(nx),nx=1,nphen), &
                (fallx(nx),nx=1,nphen),(xphendoy1(nx),nx=1,nphen)
    DO nx=1,nphen
      greenup(ilat,ivtx(nx)) = greenupx(nx)
      fall(ilat,ivtx(nx))    = fallx(nx)
      phendoy1(ilat,ivtx(nx))= xphendoy1(nx)
    ENDDO
  ENDDO

  DO np=1,mp
    ilat=(casamet%lat(np)+55.25)/0.5+1
    ilat= MIN(271,MAX(1,ilat))
    phen%phase(np) = phendoy1(ilat,veg%iveg(np))
    phen%doyphase(np,1) = greenup(ilat,veg%iveg(np)) ! DOY for greenup
    phen%doyphase(np,2) = phen%doyphase(np,1) +14    ! DOY for steady LAI
    phen%doyphase(np,3) = fall(ilat,veg%iveg(np))    ! DOY for leaf senescence
    phen%doyphase(np,4) = phen%doyphase(np,3) +14    ! DOY for minimal LAI season
    IF (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
    IF (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365

  ENDDO

END SUBROUTINE casa_readphen

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
  
! end addition (BP may2010)
  IMPLICIT NONE
!  INTEGER,        INTENT(IN)    :: mst
  TYPE (casa_biome),   INTENT(IN)    :: casabiome
  TYPE (casa_met),     INTENT(INOUT) :: casamet
TYPE (casa_flux)   , INTENT(INOUT)        :: casaflux
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
! for first time reading file *_1220.csv  (BP may2010)
  TYPE (veg_parameter_type), INTENT(IN) :: veg
  TYPE (phen_variable),   INTENT(INOUT) :: phen
  REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
  REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
  REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc

  ! local variables
  INTEGER   :: np,npt,npz
  INTEGER   :: nyearz,ivtz,istz,isoz
  REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz

  PRINT *, 'initial pool from ',TRIM(casafile%cnpipool)
  PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp
  !phen%phase = 2

  IF (initcasa==1) THEN
    OPEN(99,file=casafile%cnpipool)

    DO npt =1, mp
      SELECT CASE(icycle)
      CASE(1)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),  &
                   casapool%clitter(npt,:),casapool%csoil(npt,:)
      CASE(2)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt)
      CASE(3)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz, &
                   casamet%glai(npt),slaz,phen%phase(npt), &
                   casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
                   casapool%pplant(npt,:),casapool%plitter(npt,:),      &
                   casapool%psoil(npt,:),casapool%psoillab(npt),        &
                   casapool%psoilsorb(npt),casapool%psoilocc(npt)
      END SELECT 
      IF (ABS(patch(npt)%longitude - lonz) > 0.9 .OR. &
          ABS(patch(npt)%latitude  - latz) > 0.9) THEN
        PRINT *, 'patch(npt)%longitude, lonz:', patch(npt)%longitude, lonz
        PRINT *, 'patch(npt)%latitude,  latz:', patch(npt)%latitude,  latz
        PRINT *, 'npt = ', npt
        STOP
      ENDIF
    ENDDO
    CLOSE(99)
  ENDIF 
92  format(5(i6,3x),5(f15.6,3x),i6,3x,100(f15.6,3x))

  ! reset labile C pool,comment out by Q.Zhang 10/09/2011
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

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Nsoilmin(:) = 2.5
  ENDIF 
    
  IF (icycle >1) THEN
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
  ENDIF 

  IF (icycle >2) THEN
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
  EndIF 

END SUBROUTINE casa_init


SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  INTEGER,               INTENT(IN)    :: ktau
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

  ! Soiltype     soilnumber soil P(g P/m2)
  ! Alfisol	1	61.3
  ! Andisol	2	103.9
  ! Aridisol	3	92.8
  ! Entisol	4	136.9
  ! Gellisol	5	98.2
  ! Histosol	6	107.6
  ! Inceptisol	7	84.1
  ! Mollisol	8	110.1
  ! Oxisol	9	35.4	
  ! Spodosol	10	41.0	
  ! Ultisol	11	51.5	
  ! Vertisol	12	190.6
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

  PRINT *, 'Within casa_poolout, mp = ', mp
  nout=103
  OPEN(nout,file=casafile%cnpepool)
  PRINT *, 'Opened file ', casafile%cnpepool

  casabal%sumcbal=MIN(9999.0,MAX(-9999.0,casabal%sumcbal))
  casabal%sumnbal=MIN(9999.0,MAX(-9999.0,casabal%sumnbal))
  casabal%sumpbal=MIN(9999.0,MAX(-9999.0,casabal%sumpbal))

  DO npt =1, mp
    nso = casamet%isorder(npt)
    totpsoil(npt) = psorder(nso) *xpsoil50(nso)

    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0 
    ENDIF 

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
      casapool%pplant(npt,:) = casapool%ratiopcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
    ENDIF 

    WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
        casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
        casabiome%sla(veg%iveg(npt)), phen%phase(npt), casapool%clabile(npt), &
        casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
        casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
        casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
        casapool%plitter(npt,:), casapool%psoil(npt,:),         &
        casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
        casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
  ENDDO

  CLOSE(nout)

92    format(5(i6,',',2x),5(f15.6,',',2x),i6,',',2x,100(f15.6,',',2x))
END SUBROUTINE casa_poolout

! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters 
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  INTEGER,               INTENT(IN)    :: myear

  ! local variables
  INTEGER  npt,nout
  REAL(r_2) xyear, totGPP, totNPP

  totGPP =0.0
  totNPP =0.0
  nout=104
  xyear=1.0/FLOAT(myear)
  casabal%FCgppyear=casabal%FCgppyear * xyear
  casabal%FCnppyear=casabal%FCnppyear * xyear
  casabal%FCrsyear=casabal%FCrsyear * xyear
  casabal%FCneeyear=casabal%FCneeyear * xyear
  casabal%FNdepyear=casabal%FNdepyear * xyear
  casabal%FNfixyear=casabal%FNfixyear * xyear
  casabal%FNsnetyear=casabal%FNsnetyear * xyear
  casabal%FNupyear=casabal%FNupyear * xyear
  casabal%FNleachyear=casabal%FNleachyear * xyear
  casabal%FNlossyear=casabal%FNlossyear * xyear
  casabal%FPweayear=casabal%FPweayear * xyear
  casabal%FPdustyear=casabal%FPdustyear * xyear
  casabal%FPsnetyear=casabal%FPsnetyear * xyear
  casabal%FPupyear=casabal%FPupyear * xyear
  casabal%FPleachyear=casabal%FPleachyear * xyear
  casabal%FPlossyear=casabal%FPlossyear * xyear

  print *, 'writing CNP fluxes out to file ', casafile%cnpflux
  OPEN(nout,file=casafile%cnpflux)
    DO npt =1,mp
      SELECT CASE(icycle)
      CASE(1)

        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%Fcnppyear(npt),  &
            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &

      CASE(2)
        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
        casabal%FCnppyear(npt),casabal%FCrsyear(npt), casabal%FCneeyear(npt), &
        casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

      CASE(3)
        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
        casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
       casabal%FCnppyear(npt),casabal%FCrsyear(npt),   casabal%FCneeyear(npt),&
       casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt),&
       casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt),&
       casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt),&
       casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)

      END SELECT 
      totGPP = totGPP+casabal%Fcgppyear(npt)* casamet%areacell(npt)
      totNPP = totNPP+casabal%Fcnppyear(npt)* casamet%areacell(npt)
    ENDDO

    print *, 'totGPP global = ', totGPP*(1.0e-15)
    print *, 'totNPP global = ', totNPP*(1.0e-15)
  CLOSE(nout)
92    format(5(i6,',',2x),100(f15.6,',',2x))
END SUBROUTINE casa_fluxout

! clitterinput and csoilinput are for Julie Tang; comment out (BP apr2010)
SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,zeroflux)
USE cable_def_types_mod
USE casadimension
USE casaparm
USE casavariable
IMPLICIT NONE
TYPE (casa_flux),    INTENT(INOUT) :: casaflux
TYPE (casa_pool),    INTENT(INOUT) :: casapool
TYPE (casa_balance), INTENT(INOUT) :: casabal
LOGICAL :: zeroflux
INTEGER n

  casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
  casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
  casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp   * deltpool
  casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
  casabal%FCneeyear = casabal%FCneeyear &
                    + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
 
!  DO n=1,3
!    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
!    csoilinput(:,n) = csoilinput(:,n) + casaflux%fluxCtosoil(:,n) * deltpool
!    !csoilinput(:,n) = csoilinput(:,n)+casaflux%fluxCtolitter(:,n)*deltpool
!  ENDDO

  IF (icycle >1) THEN
    casabal%FNdepyear   = casabal%FNdepyear   + casaflux%Nmindep    * deltpool
    casabal%FNfixyear   = casabal%FNfixyear   + casaflux%Nminfix    * deltpool
    casabal%FNsnetyear  = casabal%FNsnetyear  + casaflux%Nsnet      * deltpool
    casabal%FNupyear    = casabal%FNupyear    + casaflux%Nminuptake * deltpool
    casabal%FNleachyear = casabal%FNleachyear + casaflux%Nminleach  * deltpool
    casabal%FNlossyear  = casabal%FNlossyear  + casaflux%Nminloss   * deltpool
  ENDIF 

  IF (icycle >2) THEN
    casabal%FPweayear   = casabal%FPweayear   + casaflux%Pwea       * deltpool
    casabal%FPdustyear  = casabal%FPdustyear  + casaflux%Pdep       * deltpool
    casabal%FPsnetyear  = casabal%FPsnetyear  + casaflux%Psnet      * deltpool
    casabal%FPupyear    = casabal%FPupyear    + casaflux%Plabuptake * deltpool
    casabal%FPleachyear = casabal%FPleachyear + casaflux%Pleach     * deltpool  
    casabal%FPlossyear  = casabal%FPlossyear  + casaflux%Ploss      * deltpool 
  ENDIF 

END SUBROUTINE casa_cnpflux

END MODULE casa_inout_module
