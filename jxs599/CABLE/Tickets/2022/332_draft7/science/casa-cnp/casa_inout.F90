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

SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
! mst actually not used in this routine (BP sep2010)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module, ONLY : knode_gl
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mvtype)       :: slawright
  REAL(r_2), DIMENSION(mvtype)       :: leafage,frootage,woodage
  REAL(r_2), DIMENSION(mvtype)       :: totroot
  REAL(r_2), DIMENSION(mvtype)       :: cwdage,metage,strage
  REAL(r_2), DIMENSION(mvtype)       :: micage,slowage,passage,clabileage
  REAL(r_2), DIMENSION(mvtype,mplant):: ratioCNplant
  REAL(r_2), DIMENSION(mvtype,msoil) :: ratioCNsoil,ratioCNsoilmin,ratioCNsoilmax
  REAL(r_2), DIMENSION(ms)           :: depthsoila,depthsoilb
  REAL(r_2), DIMENSION(mvtype)       :: xfNminloss, xfNminleach, xnfixrate
  REAL(r_2), DIMENSION(mvtype)       :: cleaf,cwood,cfroot,      &
                                     cmet,cstr,ccwd,          &
                                     cmic,cslow,cpass
  REAL(r_2), DIMENSION(mvtype)       :: nleaf,nwood,nfroot,      &
                                     nmet,nstr,ncwd,          &
                                     nmic,nslow,npass,xnsoilmin
  REAL(r_2), DIMENSION(mvtype)       :: xpleaf, xpwood, xpfroot, &
                                     xpmet, xpstr, xpcwd,     &
                                     xpmic,xpslow,xppass,xplab,xpsorb,xpocc
  REAL(r_2), DIMENSION(mso)       :: xkmlabp,xpsorbmax,xfPleach
  REAL(r_2), DIMENSION(mso,msoil) :: ratioNPsoil
  REAL(r_2), DIMENSION(mvtype)       :: xfherbivore,xxkleafcoldmax, xxkleafdrymax
  REAL(r_2), DIMENSION(mvtype)       :: xkuplabp
  REAL(r_2), DIMENSION(mvtype,ms)    :: fracroot 
  REAL(r_2) ::  xratioNPleafmin,xratioNPleafmax,         &
                xratioNPwoodmin,xratioNPwoodmax,         &
                xratioNPfrootmin,xratioNPfrootmax
  INTEGER :: i,iv1,nv,ns,npt,iv,is,iso
  INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10

  OPEN(101,file=casafile%cnpbiome)

  if (knode_gl==0) then    
    print *, '  '; print *, 'CASA_log:'
    print *, '  Opened file - '
    print *, '  ', trim(casafile%cnpbiome)
    print *, '  for reading cnpbiome vars.'
    print *, 'End CASA_log:'; print *, '  '
  endif

  DO i=1,3
    READ(101,*) 
  ENDDO
  
  DO nv=1,mvtype
    READ(101,*) nv0,casabiome%ivt2(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv1,casabiome%kroot(nv),casabiome%rootdepth(nv),      &
                casabiome%kuptake(nv),casabiome%krootlen(nv),         &
                casabiome%kminN(nv), casabiome%kuplabP(nv),           &
                xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                metage(nv),strage(nv),cwdage(nv),  &
                micage(nv),slowage(nv),passage(nv),clabileage(nv),slawright(nv) 
  ENDDO  

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv2, &
                casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
                casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
                casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
  ENDDO 

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv2, ratioCNplant(nv,leaf),ratioCNplant(nv,wood),   &
         ratioCNplant(nv,froot),                                         &
         casabiome%ftransNPtoL(nv,leaf), casabiome%ftransNPtoL(nv,wood), &
         casabiome%ftransNPtoL(nv,froot),                                & 
         casabiome%fracligninplant(nv,leaf),                             &
         casabiome%fracligninplant(nv,wood),                             &
         casabiome%fracligninplant(nv,froot),                            &
         ratioCNsoil(nv,mic),ratioCNsoil(nv,slow),ratioCNsoil(nv,pass),  &
         ratioCNsoilmin(nv,mic),ratioCNsoilmin(nv,slow),ratioCNsoilmin(nv,pass),  &
         ratioCNsoilmax(nv,mic),ratioCNsoilmax(nv,slow),ratioCNsoilmax(nv,pass),  &
         casabiome%glaimax(nv),casabiome%glaimin(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv3, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
                cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv4, &
         phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv), &
         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv6, &
      casabiome%ratioNCplantmin(nv,leaf),casabiome%ratioNCplantmax(nv,leaf), &
      casabiome%ratioNCplantmin(nv,wood),casabiome%ratioNCplantmax(nv,wood), &
      casabiome%ratioNCplantmin(nv,froot),casabiome%ratioNCplantmax(nv,froot), &
      xfNminloss(nv), xfNminleach(nv),xnfixrate(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                nmet(nv),nstr(nv), ncwd(nv), &
                nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv8,xratioNPleafmin,xratioNPleafmax,      &
         xratioNPwoodmin,xratioNPwoodmax,                      &
         xratioNPfrootmin,xratioNPfrootmax,                    &
         casabiome%ftransPPtoL(nv,leaf), casabiome%ftransPPtoL(nv,wood), &
         casabiome%ftransPPtoL(nv,froot)
    casabiome%ratioPcplantmin(nv,leaf)  = 1.0/(xratioNPleafmax*ratioCNplant(nv,leaf))
    casabiome%ratioPcplantmax(nv,leaf)  = 1.0/(xratioNPleafmin*ratioCNplant(nv,leaf))
    casabiome%ratioPcplantmin(nv,wood)  = 1.0/(xratioNPwoodmax*ratioCNplant(nv,wood))
    casabiome%ratioPcplantmax(nv,wood)  = 1.0/(xratioNPwoodmin*ratioCNplant(nv,wood))
    casabiome%ratioPcplantmin(nv,froot) = 1.0/(xratioNPfrootmax*ratioCNplant(nv,froot))
    casabiome%ratioPcplantmax(nv,froot) = 1.0/(xratioNPfrootmin*ratioCNplant(nv,froot))
  ENDDO

  READ(101,*)
  READ(101,*)
  DO iso=1,mso
    READ(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
                ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv10, &
         xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
         xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
  ENDDO
  CLOSE(101)

  fracroot   = 0.0
  depthsoila = 0.0
  depthsoilb = 0.0
  DO ns=1,ms
    depthsoilb(ns) = depthsoilb(ns) + soil%zse(ns)
    IF (ns==1) THEN
      depthsoila(ns) = 0.0
    ELSE
      depthsoila(ns) = depthsoilb(ns-1)
    ENDIF        
  ENDDO

  DO nv=1,mvtype
! use the value from Wright et al. (2004) (read in) instead of equation 
    casabiome%sla(nv)             = slawright(nv)
    casabiome%fraclabile(nv,leaf) = deltcasa*0.6    !1/day
    casabiome%fraclabile(nv,froot)= deltcasa*0.4    !1/day
    casabiome%fraclabile(nv,wood) = deltcasa*0.0
    casabiome%plantrate(nv,leaf)  = deltcasa/(leafage(nv)*(1.0-xfherbivore(nv)))
    casabiome%plantrate(nv,froot) = deltcasa/frootage(nv)
    casabiome%plantrate(nv,wood)  = deltcasa/woodage(nv)
    casabiome%litterrate(nv,metb) = deltcasa/metage(nv)
    casabiome%litterrate(nv,str)  = deltcasa/strage(nv)
    casabiome%litterrate(nv,cwd)  = deltcasa/cwdage(nv)
    casabiome%soilrate(nv,mic)    = deltcasa/micage(nv)
    casabiome%soilrate(nv,slow)   = deltcasa/slowage(nv)
    casabiome%soilrate(nv,pass)   = deltcasa/passage(nv)
    casabiome%xkleafcoldmax(nv)   = deltcasa * xxkleafcoldmax(nv)
    casabiome%xkleafdrymax(nv)    = deltcasa * xxkleafdrymax(nv)
    casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa 
    casabiome%kclabrate(nv)       = deltcasa/clabileage(nv)
  ENDDO


  DO npt = 1, mp
    iv1=veg%iveg(npt)
    iso=casamet%isorder(npt)
    ! The following to be commented out when coupled to CABLE
    casamet%iveg2(npt) =casabiome%ivt2(iv1)
    casamet%lnonwood(npt) = 1
    casapool%cplant(npt,wood)  = 0.0
    casapool%clitter(npt,cwd)  = 0.0
    casapool%nplant(npt,wood)  = 0.0
    casapool%nlitter(npt,cwd)  = 0.0
    casapool%pplant(npt,wood)  = 0.0
    casapool%plitter(npt,cwd)  = 0.0
    IF (casamet%iveg2(npt)==forest.or.casamet%iveg2(npt)==shrub) THEN 
      casamet%lnonwood(npt) = 0
      casapool%cplant(npt,wood)  = cwood(iv1) 
      casapool%clitter(npt,cwd)  = ccwd(iv1)
      casapool%nplant(npt,wood)  = nwood(iv1) 
      casapool%nlitter(npt,cwd)  = ncwd(iv1)
      casapool%pplant(npt,wood)  = xpwood(iv1)
      casapool%plitter(npt,cwd)  = xpcwd(iv1)
    ENDIF 
    casapool%cplant(npt,leaf)     = cleaf(iv1)
    casapool%cplant(npt,froot)    = cfroot(iv1)
    casapool%clabile(npt)         = 0.0
    casapool%clitter(npt,metb)     = cmet(iv1)
    casapool%clitter(npt,str)     = cstr(iv1)
    casapool%csoil(npt,mic)       = cmic(iv1)
    casapool%csoil(npt,slow)      = cslow(iv1)
    casapool%csoil(npt,pass)      = cpass(iv1)
    IF (icycle==1) THEN
      casapool%rationcplant(npt,:)  = 1.0/ratiocnplant(iv1,:)
    ENDIF 

    ! initializing glai in case not reading pool file (eg. during spin)
    casamet%glai(npt) = MAX(casabiome%glaimin(iv1), &
                        casabiome%sla(iv1) * casapool%cplant(npt,leaf))

    casaflux%fNminloss(npt)   = xfNminloss(iv1) 
    ! comment out by ypw 12/07/2009
    casaflux%fNminleach(npt)  = 10.0*xfNminleach(iv1) * deltcasa
    casapool%nplant(npt,leaf) = nleaf(iv1)
    casapool%nplant(npt,froot)= nfroot(iv1)
    casapool%nlitter(npt,metb) = nmet(iv1)
    casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
    casapool%nsoil(npt,mic)   = nmic(iv1)
    casapool%nsoil(npt,slow)  = nslow(iv1)
    casapool%nsoil(npt,pass)  = npass(iv1) 
    casapool%nsoilmin(npt)    = xnsoilmin(iv1) 
    casapool%pplant(npt,leaf) = xpleaf(iv1)
    casapool%pplant(npt,froot)= xpfroot(iv1) 
    casapool%plitter(npt,metb) = xpmet(iv1)
    casapool%plitter(npt,str) = cstr(iv1)*ratioPCstrfix
    casapool%psoil(npt,mic)   = xpmic(iv1)
    casapool%psoil(npt,slow)  = xpslow(iv1)
    casapool%psoil(npt,pass)  = xppass(iv1)
    casapool%psoillab(npt)    = xplab(iv1)
    casapool%psoilsorb(npt)   = xpsorb(iv1)
    casapool%psoilocc(npt)    = xpocc(iv1)
    casaflux%kmlabp(npt)      = xkmlabp(iso)
    casaflux%psorbmax(npt)    = xpsorbmax(iso)
    casaflux%fpleach(npt)     = xfPleach(iso)/(365.0)

    casapool%rationcplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
    casapool%ratiopcplant(npt,:)  = casabiome%ratioPcplantmax(iv1,:)
    casapool%rationclitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratiopclitter(npt,:) = casapool%plitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
    casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))
    casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
    casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
    casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)

  ENDDO

  IF (icycle==1) THEN
    casapool%nplant(:,:)  = casapool%cplant(:,:) * casapool%rationcplant(:,:)
  ELSE
    casapool%Nsoil(:,:)   = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)   = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
                            /(casaflux%kmlabp(:)+casapool%psoillab(:))
  ENDIF 
      

END SUBROUTINE casa_readbiome

SUBROUTINE casa_init(casabiome,casamet,casapool,casabal,veg,phen)
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
SUBROUTINE casa_cnpflux(casaflux,casabal)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_balance), INTENT(INOUT) :: casabal
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
