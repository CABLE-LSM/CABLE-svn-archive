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
!   casa_cnpflux  (zeros casabal quantites on doy 1 and updates casabal at end of biogeochem)
!CABLE_LSM:This has to be commented for offline
#define UM_BUILD YES
MODULE casa_inout_module

USE casavariable, ONLY : casafile

contains

SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
USE cable_def_types_mod
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
USE cable_common_module, ONLY : cable_runtime
USE cable_common_module, ONLY : knode_gl,cable_user 
  IMPLICIT NONE
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2), DIMENSION(mvtype)       :: leafage,frootage,woodage
  REAL(r_2), DIMENSION(mvtype)       :: totroot
  REAL(r_2), DIMENSION(mvtype)       :: cwdage,metage,strage
  REAL(r_2), DIMENSION(mvtype)       :: micage,slowage,passage,clabileage,slax
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
    INTEGER :: nv0,nv1,nv2,nv3,nv4,nv5,nv6,nv7,nv8,nv9,nv10,nv11,nv12
    REAL(r_2), DIMENSION(mvtype)       :: xxnpmax,xq10soil,xxkoptlitter,xxkoptsoil,xprodptase, &
         xcostnpup,xmaxfinelitter,xmaxcwd,xnintercept,xnslope
    REAL(r_2), DIMENSION(mso)          :: xxkplab,xxkpsorb,xxkpocc
  
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
            micage(nv),slowage(nv),passage(nv),clabileage(nv),slax(nv)
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


       casabiome%ratioNPplantmin(nv,leaf)  = xratioNPleafmin
       casabiome%ratioNPplantmax(nv,leaf)  = xratioNPleafmax
       casabiome%ratioNPplantmin(nv,wood)  = xratioNPwoodmin
       casabiome%ratioNPplantmax(nv,wood)  = xratioNPwoodmax
       casabiome%ratioNPplantmin(nv,froot) = xratioNPfrootmin
       casabiome%ratioNPplantmax(nv,froot) = xratioNPfrootmax

  ENDDO

  READ(101,*)
  READ(101,*)
  DO iso=1,mso
    READ(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
            ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass), &
            xxkplab(iso),xxkpsorb(iso),xxkpocc(iso)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv10, &
         xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
         xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
  ENDDO

    !@@@@@@@@@@@@@@@@@@@@@@@@@
    READ(101,*)
    READ(101,*)
    DO nv=1,mvtype
       READ(101,*) nv11, &
            xxnpmax(nv),xq10soil(nv),xxkoptlitter(nv),xxkoptsoil(nv),xprodptase(nv), &
            xcostnpup(nv),xmaxfinelitter(nv),xmaxcwd(nv),xnintercept(nv),xnslope(nv)
    ENDDO
    !@@@@@@@@@@@@@@@@@@@@@

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
       casabiome%sla(nv)             = slax(nv)
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

       !@@@@@@@@@@@@@@@@@
       casabiome%xnpmax(nv)          = xxnpmax(nv)
       casabiome%q10soil(nv)         = xq10soil(nv)
       casabiome%xkoptlitter(nv)     = xxkoptlitter(nv)
       casabiome%xkoptsoil(nv)       = xxkoptsoil(nv)
       casabiome%prodptase(nv)       = xprodptase(nv)/365.0   ! convert from yearly to daily
       casabiome%costnpup(nv)        = xcostnpup(nv)
       casabiome%maxfinelitter(nv)   = xmaxfinelitter(nv)
       casabiome%maxcwd(nv)          = xmaxcwd(nv)
       casabiome%nintercept(nv)      = xnintercept(nv)
       casabiome%nslope(nv)          = xnslope(nv)
       !@@@@@@@@@@@@@@
    ENDDO

    !@@@@@@@@@@@@@@
    DO ns=1,mso
       casabiome%xkplab(ns)          =  xxkplab(ns)
       casabiome%xkpsorb(ns)         =  xxkpsorb(ns)
       casabiome%xkpocc(ns)          =  xxkpocc(ns)
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
          !! vh_js !!
          IF (cable_user%CALL_POP) THEN  ! initialise very small wood pool, so POP can start from zero.
             casapool%cplant(npt,wood) = 0.01
             casapool%nplant(npt,wood)= casabiome%ratioNCplantmin(iv1,wood)* casapool%cplant(npt,wood)
             casapool%pplant(npt,wood)= casabiome%ratioPCplantmin(iv1,wood)* casapool%cplant(npt,wood)
          ENDIF
          !! vh_js

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
    casapool%plitter(npt,str) = casapool%nlitter(npt,str)/ratioNPstrfix
    casapool%psoil(npt,mic)   = xpmic(iv1)
    casapool%psoil(npt,slow)  = xpslow(iv1)
    casapool%psoil(npt,pass)  = xppass(iv1)
    casapool%psoillab(npt)    = xplab(iv1)
    casapool%psoilsorb(npt)   = xpsorb(iv1)
    casapool%psoilocc(npt)    = xpocc(iv1)
    casaflux%kmlabp(npt)      = xkmlabp(iso)
    casaflux%psorbmax(npt)    = xpsorbmax(iso)
      IF(cable_runtime%esm15_casa_inout ) THEN 
    casaflux%fpleach(npt)     = xfPleach(iso)/(365.0)
      ELSE
        casaflux%fpleach(npt)     = xfPleach(iso) 
      ENDIF  

       ! use the PFT-specific C:N:P stoichiometry
    casapool%rationcplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
       casapool%ratioNPplant(npt,:)  = casabiome%ratioNPplantmin(iv1,:)
       casapool%ratioPCplant(npt,:)  = 1.0/(ratioCNplant(iv1,:) *casabiome%ratioNPplantmin(iv1,:) )

       casapool%ratioNClitter(npt,metb) = casapool%ratioNCplant(npt,leaf)  * casabiome%ftransNPtoL(iv1,leaf)
       casapool%ratioNClitter(npt,str)  = casapool%ratioNCplant(npt,froot) * casabiome%ftransNPtoL(iv1,froot)
       casapool%ratioNClitter(npt,cwd)  = casapool%ratioNCplant(npt,wood)  * casabiome%ftransNPtoL(iv1,wood)

       casapool%ratioPClitter(npt,metb) = casapool%ratioPCplant(npt,leaf)  * casabiome%ftransPPtoL(iv1,leaf)
       casapool%ratioPClitter(npt,str)  = casapool%ratioPCplant(npt,froot) * casabiome%ftransPPtoL(iv1,froot)
       casapool%ratioPClitter(npt,cwd)  = casapool%ratioPCplant(npt,wood)  * casabiome%ftransPPtoL(iv1,wood)

       casapool%ratioNPlitter(npt,metb) = casapool%ratioNClitter(npt,metb)/(casapool%ratioPClitter(npt,metb) +1.0e-10)
       casapool%ratioNPlitter(npt,str)  = casapool%ratioNClitter(npt,str)/(casapool%ratioPClitter(npt,str) +1.0e-10)
       casapool%ratioNPlitter(npt,cwd)  = casapool%ratioNClitter(npt,cwd)/(casapool%ratioPClitter(npt,cwd) +1.0e-10)

    casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
       casapool%ratioNPsoil(npt,:)   = ratioNPsoil(iso,:)
    casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))

    casapool%ratioNCsoilmin(npt,:)   = 1.0/ratioCNsoilmax(iv1,:)
    casapool%ratioNCsoilmax(npt,:)   = 1.0/ratioCNsoilmin(iv1,:)
    casapool%ratioNCsoilnew(npt,:)   = casapool%ratioNCsoilmax(npt,:)
  ENDDO

    IF(icycle<2) THEN
    casapool%nplant(:,:)  = casapool%cplant(:,:) * casapool%rationcplant(:,:)
       casapool%Nlitter(:,:) = casapool%Clitter(:,:) * casapool%ratioNClitter(:,:)
       casapool%Nsoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioNCsoil(:,:) 
    ENDIF
    IF(icycle<3) THEN
       casapool%Pplant(:,:)  = casapool%Cplant(:,:)  * casapool%ratioPCplant(:,:)
       casapool%Plitter(:,:) = casapool%Clitter(:,:) * casapool%ratioPClitter(:,:)
       casapool%Psoil(:,:)   = casapool%Csoil(:,:)   * casapool%ratioPCsoil(:,:) 
    casapool%psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
                            /(casaflux%kmlabp(:)+casapool%psoillab(:))
  ENDIF 
      

END SUBROUTINE casa_readbiome

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

#ifndef UM_BUILD
SUBROUTINE casa_init(casabiome,casamet,casaflux,casapool,casabal,veg,phen)
    USE casadimension
    USE casaparm
    USE casavariable
    USE phenvariable
    USE cable_def_types_mod
    USE cable_io_vars_module, ONLY: landpt, patch
    USE cable_common_module, ONLY: cable_user

    IMPLICIT NONE

    TYPE (casa_biome),   INTENT(IN)    :: casabiome
    TYPE (casa_met),     INTENT(INOUT) :: casamet
    TYPE (casa_flux),    INTENT(INOUT) :: casaflux
    TYPE (casa_pool),    INTENT(INOUT) :: casapool
    TYPE (casa_balance), INTENT(INOUT) :: casabal
    TYPE (veg_parameter_type), INTENT(IN) :: veg
    TYPE (phen_variable),   INTENT(INOUT) :: phen
    REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
    REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
    REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc

    ! local variables
    INTEGER   :: np,npt,npz
    INTEGER   :: nyearz,ivtz,istz,isoz
    REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz
    LOGICAL   :: EXRST


    IF (.NOT.cable_user%casa_fromzero) THEN
       PRINT *, 'initial pool from restart file'
    ENDIF
    PRINT *, 'icycle,initcasa,mp ', icycle,initcasa,mp

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
!Ticket #204 - rm phen% clobbing here AND incorrectly so anyway

    IF (initcasa==1) THEN
       IF (.NOT.cable_user%casa_fromzero) THEN
          CALL READ_CASA_RESTART_NC (  casamet, casapool, casaflux, phen )
          WRITE(*,*)'Using input from readbiome.!!!'
          WRITE(*,*) 'initialising frac_sapwood=1 and sapwood_area = 0)'
          casaflux%frac_sapwood(:) = 1.0
          casaflux%sapwood_area(:) = 0.0
       ENDIF
    ENDIF
92  FORMAT(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),',',2x,100(f18.6,',',2x))
    WHERE(casamet%lnonwood==1) casapool%cplant(:,WOOD) = 0.0
    WHERE(casamet%lnonwood==1) casapool%nplant(:,WOOD) = 0.0
    WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0


    IF(initcasa==0) THEN
       nyearz = 1
       DO npt=1,mp
          casamet%lon(npt) = patch(npt)%longitude
          casamet%lat(npt) = patch(npt)%latitude
       ENDDO
    ENDIF

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
       WHERE(casamet%lnonwood==1) casapool%pplant(:,WOOD) = 0.0
    ENDIF
    
    casapool%cwoodprod=0.0; casapool%nwoodprod=0.0;casapool%pwoodprod=0.0

END SUBROUTINE casa_init

SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
    USE cable_common_module, ONLY: cable_user
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
    !
    ! estimated based on Yang, Post and Jain (2013)
    !   Soiltype     soilnumber soil P(g P/m2  top 50 cm)
    !   Alfisol     1       400
    !   Andisol     2       426
    !   Aridisol    3       352
    !   Entisol     4       490
    !   Gellisol    5       403
    !   Histosol    6       441
    !   Inceptisol  7       501
    !   Mollisol    8       358
    !   Oxisol      9       96
    !   Spodosol    10      364
    !   Ultisol     11      272
    !   Vertisol    12      430
    !  DATA psorder/400.0,426.0,352.0,490.0,403.0,441.0,501.0,358.0,96.0,364.0,272.0,430.0/
    !  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
    !  DATA fracpLab/0.07,0.04,0.08,0.10,0.08,0.10,0.12,0.05,0.05,0.06,0.06,0.05/
    !  DATA fracPsorb/0.30,0.44,0.69,0.53,0.37,0.14,0.24,0.32,0.15,0.21,0.17,0.35/
    !  DATA fracPocc/0.38,0.22,0.18,0.22,0.38,0.42,0.23,0.44,0.60,0.30,0.51,0.48/
    !  DATA fracPorg/0.25,0.30,0.05,0.15,0.17,0.34,0.41,0.19,0.20,0.43,0.26,0.12/
    !  DATA xpsoil50/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/

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
       IF(casamet%iveg2(npt)>0 ) THEN
    IF (icycle<2) THEN
      casapool%nplant(npt,:) = casapool%rationcplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0 
             IF(casamet%iveg2(npt)==grass) THEN
                casapool%nplant(npt,wood) = 0.0
                casapool%nlitter(npt,cwd) = 0.0
             ENDIF
    ENDIF 

    IF (icycle<3) THEN
      casabal%sumpbal(npt)   = 0.0
             casapool%pplant(npt,:)  = casapool%Nplant(npt,:)/casapool%ratioNPplant(npt,:)
             casapool%plitter(npt,:) = casapool%Nlitter(npt,:)/(casapool%ratioNPlitter(npt,:)+1.0e-10)
             casapool%psoil(npt,:)   = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
      casapool%psoillab(npt) = totpsoil(npt) *fracpLab(nso)
      casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
      casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(nso)
             IF(casamet%iveg2(npt)==grass) THEN
                casapool%pplant(npt,wood) = 0.0
                casapool%plitter(npt,cwd) = 0.0
             ENDIF
    ENDIF 
       ELSE
          casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0; casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
          casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0; casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
          casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0; casapool%psoil(npt,:) = 0.0
          casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0; casapool%psoilocc(npt) = 0.0
          casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
       ENDIF

       IF (cable_user%CALL_POP) THEN

          WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt) ,     &
               casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
               casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
               casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
               casaflux%frac_sapwood(npt), casaflux%sapwood_area(npt), &
               casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
               casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
               casapool%plitter(npt,:), casapool%psoil(npt,:),         &
               casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
               casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)


       ELSE
    WRITE(nout,92) ktau,npt,veg%iveg(npt),soil%isoilm(npt),     &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
        casamet%areacell(npt)*(1.0e-9),casamet%glai(npt),       &
               casabiome%sla(veg%iveg(npt)), phen%phase(npt), &
               phen%doyphase(npt,3), phen%phen(npt), phen%aphen(npt), &
               casapool%clabile(npt), &
        casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
        casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
        casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
        casapool%plitter(npt,:), casapool%psoil(npt,:),         &
        casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
        casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
       ENDIF

  ENDDO

  CLOSE(nout)

92  FORMAT(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),100(f18.6,',',2x))
END SUBROUTINE casa_poolout

! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
SUBROUTINE casa_fluxout(myear,veg,soil,casabal,casamet)

  USE cable_def_types_mod
USE cable_common_module, ONLY : knode_gl,cable_user, cable_runtime 
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
    casabal%FCrmleafyear=casabal%FCrmleafyear * xyear
    casabal%FCrmwoodyear=casabal%FCrmwoodyear * xyear
    casabal%FCrmrootyear=casabal%FCrmrootyear * xyear
    casabal%FCrgrowyear=casabal%FCrgrowyear * xyear
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
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
            casabal%Fcnppyear(npt),  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt)  ! ,           &

      CASE(2)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
            casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt),  &
               casabal%FCnppyear(npt),                                 &
               casabal%FCrmleafyear(npt),casabal%FCrmwoodyear(npt),     &
               casabal%FCrmrootyear(npt),casabal%FCrgrowyear(npt),     &
               casabal%FCrsyear(npt), casabal%FCneeyear(npt),          &
                                !        clitterinput(npt,:),csoilinput(npt,:), &
        casabal%FNdepyear(npt),casabal%FNfixyear(npt),casabal%FNsnetyear(npt), &
        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)

      CASE(3)
          WRITE(nout,*) myear,npt,veg%iveg(npt),soil%isoilm(npt), &
        casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt),  &
        casamet%areacell(npt)*(1.0e-9),casabal%Fcgppyear(npt), &
               casabal%FCnppyear(npt),                                  &
               casabal%Fcrmleafyear(npt),casabal%Fcrmwoodyear(npt),     &
               casabal%Fcrmrootyear(npt),casabal%Fcrgrowyear(npt),     &
               casabal%FCrsyear(npt),   casabal%FCneeyear(npt),         &
                                !        clitterinput(npt,:),csoilinput(npt,:), &
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
#endif

  SUBROUTINE casa_cnpflux(casaflux,casapool,casabal,zeroflux)
  USE cable_def_types_mod
USE cable_common_module, ONLY : knode_gl,cable_user 
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
    TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
    LOGICAL :: zeroflux
  INTEGER n

IF( .NOT. cable_runtime%esm15_casa_cnpflux) THEN 
    IF(zeroflux) THEN
       casabal%FCgppyear    = 0.0
       casabal%FCrpyear     = 0.0
       casabal%FCrmleafyear = 0.0
       casabal%FCrmwoodyear = 0.0
       casabal%FCrmrootyear = 0.0
       casabal%FCrgrowyear  = 0.0
       casabal%FCnppyear    = 0.0
       casabal%FCrsyear     = 0.0
       casabal%FCneeyear    = 0.0
       casabal%dCdtyear    = 0.0


       casabal%FNdepyear    = 0.0
       casabal%FNfixyear    = 0.0
       casabal%FNsnetyear   = 0.0
       casabal%FNupyear     = 0.0
       casabal%FNleachyear  = 0.0
       casabal%FNlossyear   = 0.0

       casabal%FPweayear   = 0.0
       casabal%FPdustyear  = 0.0
       casabal%FPsnetyear  = 0.0
       casabal%FPupyear    = 0.0
       casabal%FPleachyear = 0.0
       casabal%FPlossyear  = 0.0

       casaflux%FluxCtohwp = 0.0
       casaflux%FluxNtohwp = 0.0
       casaflux%FluxPtohwp = 0.0
       casaflux%FluxCtoclear = 0.0
       casaflux%FluxNtoclear = 0.0
       casaflux%FluxPtoclear = 0.0
       casaflux%CtransferLUC = 0.02

       casaflux%Crp(:)   = casaflux%Crmplant(:,leaf) + casaflux%Crmplant(:,wood) + casaflux%Crmplant(:,froot) + casaflux%Crgplant(:)
       casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
       casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
       casabal%FCrmleafyear(:)  = casabal%FCrmleafyear(:)  + casaflux%Crmplant(:,leaf)    * deltpool
       casabal%FCrmwoodyear(:)  = casabal%FCrmwoodyear(:)  + casaflux%Crmplant(:,wood)    * deltpool
       casabal%FCrmrootyear(:)  = casabal%FCrmrootyear(:)  + casaflux%Crmplant(:,froot)   * deltpool
       casabal%FCrgrowyear      = casabal%FCrgrowyear      + casaflux%Crgplant            * deltpool
       ! change made ypwang 17-nov-2013 to accoutn for change in labile carbon pool  size
       casabal%FCnppyear        = casabal%FCnppyear + (casaflux%Cnpp+casapool%dClabiledt)   * deltpool
       casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
       casabal%FCneeyear = casabal%FCneeyear &
            + (casaflux%Cnpp+casapool%dClabiledt-casaflux%Crsoil) * deltpool
       casabal%dCdtyear =  casabal%dCdtyear + (casapool%Ctot-casapool%Ctot_0)*deltpool

       !  DO n=1,3
       !    clitterinput(:,n)= clitterinput(:,n) + casaflux%kplant(:,n) * casapool%cplant(:,n) * deltpool
       !    csoilinput(:,n) = csoilinput(:,n) + casaflux%fluxCtosoil(:,n) * deltpool
       !    !csoilinput(:,n) = csoilinput(:,n)+casaflux%fluxCtolitter(:,n)*deltpool
       !  ENDDO

    ELSE
  casaflux%Crp = 0.

  casabal%FCgppyear = casabal%FCgppyear + casaflux%Cgpp   * deltpool
  casabal%FCrpyear  = casabal%FCrpyear  + casaflux%Crp    * deltpool
  casabal%FCnppyear = casabal%FCnppyear + casaflux%Cnpp   * deltpool
  casabal%FCrsyear  = casabal%FCrsyear  + casaflux%Crsoil * deltpool
  casabal%FCneeyear = casabal%FCneeyear &
                    + (casaflux%Cnpp-casaflux%Crsoil) * deltpool
 
ENDIF

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

#ifndef UM_BUILD
  SUBROUTINE WRITE_CASA_RESTART_NC ( casamet, casapool, casaflux, phen, CASAONLY )

    USE casavariable, ONLY : casa_met, casa_pool, casa_flux, icycle, mplant, mlitter, msoil
    USE cable_common_module
  USE casa_ncdf_module, ONLY: HANDLE_ERR
   
    USE cable_def_types_mod, ONLY: met_type, mp
    USE phenvariable
    USE netcdf

    IMPLICIT NONE


    TYPE (casa_met),  INTENT(IN) :: casamet
    TYPE (casa_pool),  INTENT(IN) :: casapool
    TYPE (casa_flux),           INTENT(IN) :: casaflux
    TYPE (phen_variable),       INTENT(IN) :: phen

    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS
    INTEGER*4   :: FILE_ID, land_ID, plnt_ID, litt_ID, soil_ID, i
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50

    ! ! 1 dim arrays (npt )
    ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
    !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
    ! ! 2 dim arrays (npt,mplant)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
    ! ! 2 dim arrays (npt,mlitter)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
    ! ! 2 dim arrays (npt,msoil)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)

    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(12) :: A1
    CHARACTER(len=20),DIMENSION(2) :: AI1
    ! 2 dim arrays (npt,mplant)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,mlitter)
    CHARACTER(len=20),DIMENSION(3) :: A3
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3) :: A4
    INTEGER*4 :: VID1(SIZE(A1)), VIDI1(SIZE(AI1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))

    mp4=INT(mp,fmp4)
    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'glai'
    A1(4) = 'clabile'
    A1(5) = 'psoillab'
    A1(6) = 'psoilsorb'
    A1(7) = 'psoilocc'
    A1(8) = 'frac_sapwood'
    A1(9) = 'sapwood_area'
    A1(10) = 'phen'
    A1(11) = 'aphen'
    A1(12) = 'nsoilmin'

    AI1(1) = 'phase'
    AI1(2) = 'doyphase3'


    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'

    ! Get File-Name
    WRITE(CYEAR, FMT='(I4)') CurYear + 1

    IF (LEN( TRIM(casafile%cnpepool) ) .GT. 0) THEN
       fname=TRIM(casafile%cnpepool)
    ELSE
       fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
    ENDIF
    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) 'writing casa restart', fname
    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", "01/01/"//CYEAR  )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle", icycle  )
    IF ( CASAONLY ) THEN
       dum = 'CASA-ONLY run'
    ELSE
       dum = 'CABLE-CASA coupled run'
    ENDIF
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp4     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(AI1)
       STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID/),VIDI1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A3)
       STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID/),VID3(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A4)
       STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID/),VID4(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(1), casamet%lat )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(2), casamet%lon )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT VARS
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(3), casamet%glai )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(4), casapool%clabile )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    STATUS = NF90_PUT_VAR(FILE_ID, VID1(8), casaflux%frac_sapwood )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(9), casaflux%sapwood_area )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), phen%phen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), phen%aphen )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casapool%Nsoilmin )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), phen%phase )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(2), phen%doyphase(:,3) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter  )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(5), casapool%psoillab )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(6), casapool%psoilsorb )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(7), casapool%psoilocc )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant  )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  END SUBROUTINE WRITE_CASA_RESTART_NC

#ifndef UM_BUILD
  SUBROUTINE READ_CASA_RESTART_NC (  casamet, casapool, casaflux,phen )

    USE CASAVARIABLE
    USE phenvariable
    USE CABLE_COMMON_MODULE
    USE casa_ncdf_module, ONLY: HANDLE_ERR
    USE CABLE_DEF_TYPES_MOD, ONLY: MET_TYPE, r_2, mp
    USE netcdf

    IMPLICIT NONE

    !INTEGER, INTENT(in)    :: YEAR
    TYPE (casa_met) , INTENT(inout) :: casamet
    TYPE (casa_pool), INTENT(inout) :: casapool
    TYPE (casa_flux), INTENT(inout) :: casaflux
    TYPE (phen_variable),       INTENT(INOUT) :: phen

    INTEGER*4 :: mp4
    INTEGER*4, PARAMETER   :: pmp4 =0
    INTEGER, PARAMETER   :: fmp4 = KIND(pmp4)
    INTEGER*4   :: STATUS, i
    INTEGER*4   :: FILE_ID, dID, land_dim, mp_dim, ml_dim, ms_dim, mw_dim
    CHARACTER :: FRST_IN*99, CYEAR*4, CDATE*12, RSTDATE*12, FNAME*99

    ! ! 1 dim arrays (npt )
    ! CHARACTER(len=20),DIMENSION(7), PARAMETER :: A1 = (/ 'latitude', 'longitude', 'glai', &
    !      'clabile', 'psoillab','psoilsorb','psoilocc' /)
    ! ! 2 dim arrays (npt,mplant)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A2 = (/ 'cplant' , 'nplant' , 'pplantc' /)
    ! ! 2 dim arrays (npt,mlitter)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A3 = (/ 'clitter', 'nlitter', 'plitter' /)
    ! ! 2 dim arrays (npt,msoil)
    ! CHARACTER(len=20),DIMENSION(3), PARAMETER :: A4 = (/ 'csoil', 'nsoil', 'psoil' /)
    REAL(r_2), DIMENSION(mp)          :: LAT, LON, TMP
    REAL(r_2)                         :: TMP2(mp,mplant),TMP3(mp,mlitter),TMP4(mp,msoil)

    ! 1 dim arrays (npt )
    CHARACTER(len=20),DIMENSION(12) :: A1
    CHARACTER(len=20),DIMENSION(2) :: AI1
    ! 2 dim arrays (npt,mplant)
    CHARACTER(len=20),DIMENSION(3) :: A2
    ! 2 dim arrays (npt,mlitter)
    CHARACTER(len=20),DIMENSION(3) :: A3
    ! 2 dim arrays (npt,msoil)
    CHARACTER(len=20),DIMENSION(3) :: A4
    ! 2-d array (npt,mwood)
    CHARACTER(len=20),DIMENSION(3) :: A5

    INTEGER :: VID1(SIZE(A1)), VID2(SIZE(A2)), VID3(SIZE(A3)), VID4(SIZE(A4))
    LOGICAL            ::  EXISTFILE, EXISTFILE1
    mp4=INT(mp,fmp4)
    A1(1) = 'latitude'
    A1(2) = 'longitude'
    A1(3) = 'glai'
    A1(4) = 'clabile'
    A1(5) = 'psoillab'
    A1(6) = 'psoilsorb'
    A1(7) = 'psoilocc'
    A1(8) = 'frac_sapwood'
    A1(9) = 'sapwood_area'
    A1(10) = 'phen'
    A1(11) = 'aphen'
    A1(12) = 'nsoilmin'

    AI1(1) = 'phase'
    AI1(2) = 'doyphase3'

    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'

    A5(1) = 'cwoodprod'
    A5(2) = 'nwoodprod'
    A5(3) = 'pwoodprod'


    fname =  TRIM(casafile%cnpipool)
    INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE )
    IF (EXISTFILE) THEN
       STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       PRINT *, 'initial pool from restart file: ', fname
    ELSE
       WRITE(*,*) 'CASA restart file:', TRIM(fname), ' does not exist'
       fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
            '_casa_rst.nc'
       INQUIRE( FILE=TRIM(fname), EXIST=EXISTFILE1 )
       IF (EXISTFILE1) THEN
          STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          PRINT *, 'initial pool from restart file: ', fname
       ELSE
          WRITE(*,*) 'CASA restart file:', TRIM(fname), ' does not exist either'
          WRITE(*,*) 'Set cable_user%CASA_fromZero to true to initialise without restart file.'
          WRITE(*,*) 'Otherwise set casafile%cnpipool to netcdf restart file name in cable.nml'
          STOP
       ENDIF
    ENDIF

    ! TIME
    STATUS = NF90_GET_ATT( FILE_ID, NF90_GLOBAL, "Valid restart date", RSTDATE )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$
    WRITE(CYEAR, FMT="(I4)") CurYear
    CDATE = '01/01/'//CYEAR
    ! compare current year with restart year (only for non-site type met data)
    IF ( CDATE .NE. RSTDATE .AND. &
         TRIM(cable_user%MetType).NE.'' .AND. TRIM(cable_user%MetType).NE.'site' ) THEN
       WRITE(*,*)"Restart Date in rst file doesn't match start date of Run!"
       WRITE(*,*)"File: "//RSTDATE//' Run: '//CDATE
       ! STOP
    ENDIF

    ! DIMS
    STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'mplant', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mp_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'mlitter', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ml_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_DIMID( FILE_ID, 'msoil', dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=ms_dim )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    if(l_landuse) then
       STATUS = NF90_INQ_DIMID( FILE_ID, 'mwood', dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mw_dim )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    endif

    IF ( land_dim .NE. SIZE(casamet%lon) .OR. mp_dim .NE. mplant .OR. &
         ml_dim   .NE. mlitter             .OR. ms_dim .NE. msoil ) THEN
       WRITE(*,*)"Dimension misfit!"
       WRITE(*,*)"Restart file      Run"
       WRITE(*,*)"# points  ",land_dim,"     ",SIZE(casamet%lon)
       WRITE(*,*)"# mplant  ",mp_dim,"     ",mplant
       WRITE(*,*)"# mlitter ",ml_dim,"     ",mlitter
       WRITE(*,*)"# msoil   ",ms_dim,"     ",msoil
       STOP
    ENDIF

    ! LAT & LON
    STATUS = NF90_INQ_VARID( FILE_ID, A1(1), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LAT )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    STATUS = NF90_INQ_VARID( FILE_ID, A1(2), dID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_GET_VAR( FILE_ID, dID, LON )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! CHECK FOR VALID LONS

    ! READ 1-dimensional fields
    DO i = 3, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A1(i)))
       CASE ('glai'      ) ; casamet%glai       = TMP
       CASE ('clabile'   ) ; casapool%clabile   = TMP
       CASE ('frac_sapwood' ) ; casaflux%frac_sapwood  = TMP
       CASE ( 'sapwood_area' ) ; casaflux%sapwood_area  = TMP
       CASE ( 'phen' ) ; phen%phen  = TMP
       CASE ( 'aphen' ) ; phen%aphen  = TMP
       CASE ( 'nsoilmin' ) ; casapool%Nsoilmin  = TMP
       END SELECT
    END DO
    IF (icycle==3) THEN
       DO i = 3, SIZE(A1)
          STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A1(i)))
          CASE ('psoillab'  ) ; casapool%psoillab  = TMP
          CASE ('psoilsorb' ) ; casapool%psoilsorb = TMP
          CASE ('psoilocc'  ) ; casapool%psoilocc  = TMP
          END SELECT
       END DO
    ENDIF

    DO i = 1, SIZE(AI1)
       STATUS = NF90_INQ_VARID( FILE_ID, AI1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(AI1(i)))
       CASE ( 'phase' ) ; phen%phase  = TMP
       CASE ( 'doyphase3' ) ; phen%doyphase(:,3)  = TMP
       END SELECT
    END DO

    ! READ 2-dimensional fields (mplant)
    DO i = 1, SIZE(A2)
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('cplant' ) ; casapool%cplant = TMP2
       CASE ('nplant' ) ; casapool%nplant = TMP2
       END SELECT
    END DO


    IF (icycle==3) THEN
       DO i = 1, SIZE(A2)
          STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A2(i)))
          CASE ('pplant' ) ; casapool%pplant = TMP2
          END SELECT
       END DO
    ENDIF

    ! READ 2-dimensional fields (mlitter)
    DO i = 1, SIZE(A3)
       STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A3(i)))
       CASE ('clitter' ) ; casapool%clitter = TMP3
       CASE ('nlitter' ) ; casapool%nlitter = TMP3
       END SELECT
    END DO

    IF (icycle==3) THEN

       DO i = 1, SIZE(A3)
          STATUS = NF90_INQ_VARID( FILE_ID, A3(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP3 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A3(i)))
          CASE ('plitter' ) ; casapool%plitter = TMP3
          END SELECT
       END DO


    ENDIF

    ! READ 2-dimensional fields (msoil)
    DO i = 1, SIZE(A4)
       STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       SELECT CASE ( TRIM(A4(i)))
       CASE ('csoil' ) ; casapool%csoil = TMP4
       CASE ('nsoil' ) ; casapool%nsoil = TMP4
       END SELECT
    END DO

    if(l_landuse) then

       DO i = 1, SIZE(A5)
          STATUS = NF90_INQ_VARID( FILE_ID, A5(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
           STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          SELECT CASE ( TRIM(A4(i)))
             CASE ('cwoodprod' ) ; casapool%cwoodprod = TMP4
             CASE ('nwoodprod' ) ; casapool%nwoodprod = TMP4
             CASE ('pwoodprod' ) ; casapool%pwoodprod = TMP4
          END SELECT
        END DO

    endif

    IF (icycle==3) THEN
       DO i = 1, SIZE(A4)
          STATUS = NF90_INQ_VARID( FILE_ID, A4(i), dID )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_GET_VAR( FILE_ID, dID, TMP4 )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          SELECT CASE ( TRIM(A4(i)))
          CASE ('psoil' ) ; casapool%psoil = TMP4
          END SELECT
       END DO
    ENDIF

    STATUS = NF90_CLOSE( FILE_ID )

  END SUBROUTINE READ_CASA_RESTART_NC
#endif
  SUBROUTINE WRITE_CASA_OUTPUT_NC ( veg, casamet, casapool, casabal, casaflux, &
       CASAONLY, ctime, FINAL )

    USE CASAVARIABLE
    USE CABLE_COMMON_MODULE
  USE casa_ncdf_module, ONLY: HANDLE_ERR


    USE cable_def_types_mod, ONLY: veg_parameter_type

    USE netcdf

    IMPLICIT NONE

    TYPE (casa_met) ,   INTENT(in) :: casamet
    TYPE (casa_pool),   INTENT(in) :: casapool
    TYPE (casa_balance),INTENT(in) :: casabal
    TYPE (casa_flux),   INTENT(in) :: casaflux
    TYPE (veg_parameter_type),  INTENT(IN) :: veg  ! vegetation parameters

    INTEGER   :: STATUS, ctime
    INTEGER   :: land_ID, plnt_ID, litt_ID, soil_ID, t_ID, i
    LOGICAL   :: CASAONLY, FINAL
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! ! 1 dim arrays (mp )
    ! CHARACTER(len=20),DIMENSION(2), PARAMETER :: A0 = (/ 'latitude', 'longitude' /)
    ! ! 2 dim arrays (mp,t)
    ! CHARACTER(len=20),DIMENSION(44),PARAMETER :: A1 = (/ 'glai', 'clabile',      &
    !      'psoillab','psoilsorb','psoilocc', 'sumcbal','sumnbal','sumpbal','Cgpp',&
    !      'Cnpp','stemnpp','Crp','Crgplant','Nminfix','Plabuptake','Clabloss',    &
    !      'fraclabile','Cnep','Crsoil','Nmindep','Nminloss','Nminleach',          &
    !      'Nupland','Nlittermin','Nsmin','Nsnet','fNMinloss','fNMinleach','Pdep', &
    !      'pwea','Pleach','Ploss','Pupland','Plittermin','Psmin','Psimm','Psnet', &
    !      'fPleach','kPlab','kPsorb','kpocc','kmlabP','Psorbmax','FluxCtoco2'/)
    ! ! 3 dim arrays (mp,mplant,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A2 = (/ 'cplant' , 'nplant' ,   &
    !      'pplantc','fracCalloc','fracNalloc','fracPalloc','kplant','Crmplant'/)
    ! ! 3 dim arrays (mp,mlitter,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A3 = (/ 'clitter', 'nlitter',   &
    !      'plitter','klitter','fromL2CO2','FluxCtolitter','FluxNtolitter',        &
    !      'FluxPtolitter' /)
    ! ! 3 dim arrays (mp,msoil,t)
    ! CHARACTER(len=20),DIMENSION(8), PARAMETER :: A4 = (/ 'csoil','nsoil','psoil',&
    !      'ksoil','fromStoCO2','FluxCtosoil','FluxNtosoil','FluxPxtosoil'/)

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim arrays (mp,t)
    CHARACTER(len=20),DIMENSION(51):: A1
    ! 3 dim arrays (mp,mplant,t)
    CHARACTER(len=20),DIMENSION(8) :: A2
    ! 3 dim arrays (mp,mlitter,t)
    CHARACTER(len=20),DIMENSION(8) :: A3
    ! 3 dim arrays (mp,msoil,t)
    CHARACTER(len=20),DIMENSION(8) :: A4

    ! 4 dim arrays (mp,mlitter,mplant,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A5 = (/ 'fromPtoL'/)
    ! 4 dim arrays (mp,msoil,mlitter,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A6 = (/ 'fromLtoS'/)
    ! 4 dim arrays (mp,msoil,msoil,t)
    CHARACTER(len=20),DIMENSION(1), PARAMETER :: A7 = (/ 'fromStoS'/)

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VID2(SIZE(A2)),VID3(SIZE(A3))
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5)),VID6(SIZE(A6)),VID7(SIZE(A7))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    LOGICAL   :: EXRST
    CHARACTER(len=50) :: RecordDimName
    ! temporary fixer: ypwang 10-6-2021
    REAL(r_2), DIMENSION(mp,mplant)   :: tempFCrmyear
    REAL(R_2), DIMENSION(mp)          :: tempFCrgyear

    tempFCrmyear(:,1)  = casabal%FCrmleafyear(:)   ! leaf
    tempFCrmyear(:,2)  = casabal%FCrmwoodyear(:)   ! wood
    tempFCrmyear(:,3)  = casabal%FCrmrootyear(:)   ! froot
    tempFCrgyear(:)    = casabal%FCrgrowyear(:)



    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'glai'
    A1(2) = 'clabile'
    A1(3) = 'psoillab'
    A1(4) = 'psoilsorb'
    A1(5) = 'psoilocc'
    A1(6) = 'sumcbal'
    A1(7) = 'sumnbal'
    A1(8) = 'sumpbal'
    A1(9) = 'Cgpp'
    A1(10) = 'Cnpp'
    A1(11) = 'stemnpp'
    A1(12) = 'Crp'
    A1(13) = 'Crgplant'
    A1(14) = 'Nminfix'
    A1(15) = 'Plabuptake'
    A1(16) = 'Clabloss'
    A1(17) = 'fraclabile'
    A1(18) = 'Cnep'
    A1(19) = 'Crsoil'
    A1(20) = 'Nmindep'
    A1(21) = 'Nminloss'
    A1(22) = 'Nminleach'
    A1(23) = 'Nupland'
    A1(24) = 'Nlittermin'
    A1(25) = 'Nsmin'
    A1(26) = 'Nsimm'
    A1(27) = 'Nsnet'
    A1(28) = 'fNMinloss'
    A1(29) = 'Pdep'
    A1(30) = 'pwea'
    A1(31) = 'Pleach'
    A1(32) = 'Ploss'
    A1(33) = 'Pupland'
    A1(34) = 'Plittermin'
    A1(35) = 'Psmin'
    A1(36) = 'Psimm'
    A1(37) = 'Psnet'
    A1(38) = 'fPleach'
    A1(39) = 'kPlab'
    A1(40) = 'kPsorb'
    A1(41) = 'kpocc'
    A1(42) = 'kmlabP'
    A1(43) = 'Psorbmax'
    A1(44) = 'FluxCtoco2'
    A1(45) = 'FCgppyear'
    A1(46) = 'FCrpyear'
    A1(47) = 'FCnppyear'
    A1(48) = 'FCrsyear'
    A1(49) = 'FCNeeyear'
    A1(50) = 'vcmax'
    A1(51) = 'Nsoilmin'

    A2(1) = 'cplant'
    A2(2) = 'nplant'
    A2(3) = 'pplant'
    A2(4) = 'fracCalloc'
    A2(5) = 'fracNalloc'
    A2(6) = 'fracPalloc'
    A2(7) = 'kplant'
    A2(8) = 'Crmplant'

    A3(1) = 'clitter'
    A3(2) = 'nlitter'
    A3(3) = 'plitter'
    A3(4) = 'klitter'
    A3(5) = 'fromL2CO2'
    A3(6) = 'FluxCtolitter'
    A3(7) = 'FluxNtolitter'
    A3(8) = 'FluxPtolitter'

    A4(1) = 'csoil'
    A4(2) = 'nsoil'
    A4(3) = 'psoil'
    A4(4) = 'ksoil'
    A4(5) = 'fromStoCO2'
    A4(6) = 'FluxCtosoil'
    A4(7) = 'FluxNtosoil'
    A4(8) = 'FluxPxtosoil'


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       IF (TRIM(cable_user%MetType).NE.'' ) THEN

          WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          IF (CABLE_USER%YEARSTART.LT.1000.AND.CABLE_USER%YEAREND.LT.1000) THEN
             WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          ELSEIF (CABLE_USER%YEARSTART.LT.1000) THEN
             WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
          ENDIF
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
               TRIM(dum)//'_casa_out.nc'
       ELSE
          ! site data
          fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_casa_out.nc'
       ENDIF
       INQUIRE( FILE=TRIM( fname ), EXIST=EXRST )
       EXRST = .FALSE.
       IF ( EXRST ) THEN
          STATUS = NF90_open(fname, mode=nf90_write, ncid=FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          CALL1 = .FALSE.

          STATUS = nf90_inq_dimid(FILE_ID, 'time', t_id)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          !CRM        status = nf90_inquire_dimension(FILE_ID, t_id,name = RecordDimName, len = CNT)
          !CRM        if (status /= nf90_noerr) call handle_err(status)
          !CRM        CNT = CNT+1

          STATUS = nf90_inq_varid(FILE_ID, 'time', VIDTime)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          DO i = 1, SIZE(A0)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A0(i)),VID0(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A1)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A1(i)), VID1(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A2)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A2(i)) , VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A3)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A3(i)) ,VID3(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A4)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A4(i)) ,VID4(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A5)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A5(i)), VID5(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A6)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A6(i)), VID6(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A7)
             STATUS = nf90_inq_varid(FILE_ID,TRIM(A7(i)),VID7(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

       ELSE
          ! Create NetCDF file:
          STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Put the file in define mode:
          STATUS = NF90_redef(FILE_ID)

          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Icycle"   , icycle  )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )
          IF ( CASAONLY ) THEN
             dum = 'CASA-ONLY run'
          ELSE
             dum = 'CABLE-CASA coupled run'
          ENDIF
          STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "Run-Type", TRIM(dum) )

          ! Define dimensions:
          ! Land (number of points)
          STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'mplant' , mplant , plnt_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'mlitter', mlitter, litt_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'msoil'  , msoil  , soil_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          STATUS = NF90_def_dim(FILE_ID, 'time'   , NF90_UNLIMITED, t_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          ! Define variables
          STATUS = NF90_def_var(FILE_ID,'time' ,NF90_INT,(/t_ID/),VIDtime )
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

          DO i = 1, SIZE(A0)
             STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A1)
             STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID,t_ID/),VID1(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,plnt_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A3)
             STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,litt_ID,t_ID/),VID3(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A4)
             STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,soil_ID,t_ID/),VID4(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A5)
             STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT, &
                  (/land_ID,litt_ID,plnt_ID,t_ID/),VID5(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A6)
             STATUS = NF90_def_var(FILE_ID,TRIM(A6(i)) ,NF90_FLOAT, &
                  (/land_ID,soil_ID,litt_ID,t_ID/),VID6(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          DO i = 1, SIZE(A7)
             STATUS = NF90_def_var(FILE_ID,TRIM(A7(i)) ,NF90_FLOAT, &
                  (/land_ID,soil_ID,soil_ID,t_ID/),VID7(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO

          ! End define mode:
          STATUS = NF90_enddef(FILE_ID)
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


          ! PUT LAT / LON ( mp )
          STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), casamet%lat )
          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

          STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), casamet%lon )
          IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

          CALL1 = .FALSE.
       ENDIF !( EXRST )
    ENDIF

    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    SELECT CASE(icycle)
    CASE(1)
       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), REAL(casaflux%Cnpp), start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casaflux%Crsoil,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    CASE(2)


       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casabal%FNfixyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       ! STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       ! IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casabal%FPupyear, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casabal%FCrsyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casabal%FNdepyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casabal%FNlossyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casabal%FNleachyear,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casabal%FNleachyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casabal%FNsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casabal%FPdustyear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casabal%FPweayear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casabal%FPleachyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casabal%FPlossyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casabal%FPsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoco2,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    CASE(3)
       ! PUT 2D VARS ( mp, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), casamet%glai,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), casapool%clabile,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), casapool%psoillab,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), casapool%psoilsorb,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), casapool%psoilocc,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), casabal%sumcbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), casabal%sumnbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), casabal%sumpbal,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), casaflux%Cgpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), casaflux%Cnpp,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(11), casaflux%stemnpp,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(12), casaflux%Crp,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(13), tempFCrgyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(14), casabal%FNfixyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       !  IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       !  STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casaflux%Nminuptake, start=(/ 1, CNT /), count=(/ mp, 1 !/) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(15), casabal%FPupyear, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(16), casaflux%Clabloss,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(17), casaflux%fracClabile,start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(18), casaflux%Cnep,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(19), casabal%FCrsyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(20), casabal%FNdepyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(21), casabal%FNlossyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(22), casabal%FNleachyear,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(23), casabal%FNupyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(24), casaflux%Nlittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(25), casaflux%Nsmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(26), casaflux%Nsimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(27), casabal%FNsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(28), casaflux%fNminloss,  start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(29), casabal%FPdustyear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(30), casabal%FPweayear,       start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(31), casabal%FPleachyear,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(32), casabal%FPlossyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(33), casaflux%Pupland,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(34), casaflux%Plittermin, start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(35), casaflux%Psmin,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(36), casaflux%Psimm,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(37), casabal%FPsnetyear,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(38), casaflux%fPleach,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(39), casaflux%kplab,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(40), casaflux%kpsorb,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(41), casaflux%kpocc,      start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(42), casaflux%kmlabP,     start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(43), casaflux%Psorbmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(44), casaflux%FluxCtoCo2,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(45), casabal%FCgppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(46), casabal%FCrpyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(47), casabal%FCnppyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(48), casabal%FCrsyear,    start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(49), casabal%FCneeyear,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(50), veg%vcmax,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID1(51), casapool%Nsoilmin,   start=(/ 1, CNT /), count=(/ mp, 1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


       ! PUT 3D VARS ( mp, mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), casapool%cplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2), casapool%nplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), casapool%pplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), casaflux%fracCalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(5), casaflux%fracNalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(6), casaflux%fracPalloc,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(7), casaflux%kplant,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(8), tempFCrmyear,   &
            start=(/ 1,1,CNT /), count=(/ mp,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), casapool%clitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(2), casapool%nlitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(3), casapool%plitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(4), casaflux%klitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(5), casaflux%fromLtoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(6), casaflux%FluxCtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(7), casaflux%FluxNtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID3(8), casaflux%FluxPtolitter,   &
            start=(/ 1,1,CNT /), count=(/ mp,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 3D VARS ( mp, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), casapool%csoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), casapool%nsoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), casapool%psoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), casaflux%ksoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), casaflux%fromStoCO2,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), casaflux%FluxCtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), casaflux%FluxNtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID4(8), casaflux%FluxPtosoil,   &
            start=(/ 1,1,CNT /), count=(/ mp,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, mlitter,mplant, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), casaflux%fromPtoL,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,mlitter,mplant,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, mlitter, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID6(1), casaflux%fromLtoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,mlitter,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       ! PUT 4D VARS ( mp, msoil, msoil, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID7(1), casaflux%fromStoS,   &
            start=(/ 1,1,1,CNT /), count=(/ mp,msoil,msoil,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    END SELECT

    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " Casa Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_CASA_OUTPUT_NC
#endif

END MODULE casa_inout_module
