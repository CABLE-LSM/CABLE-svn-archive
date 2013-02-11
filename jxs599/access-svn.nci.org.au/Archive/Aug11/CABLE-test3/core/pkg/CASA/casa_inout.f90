! casa_inout.f90
!
! Model development by YingPing Wang, CSIRO Marine and Atmospheric Research.
! Coupling to Mk3L by Bernard Pak,    CSIRO Marine and Atmospheric Research.
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! the following routines are used when "casacnp" is coupled to "cable"
!   casa_readbiome
!   casa_readphen
!   casa_readpoint
!   casa_init
!   casa_poolout
!   casa_cnpflux  (not used?)
!   biogeochem

SUBROUTINE casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
! added mvtype (=mvt) and mstype (=mst) to define_dimensions (BP sep2010)
! mst actually not used in this routine (BP sep2010)
!SUBROUTINE casa_readbiome(mvt,mst,veg,soil, &
!                          casabiome,casapool,casaflux,casamet,phen)
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
!  INTEGER(i_d),               INTENT(IN)    :: mvt,mst
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
  REAL(r_2), DIMENSION(mvtype)       :: micage,slowage,passage
  REAL(r_2), DIMENSION(mvtype,mplant):: ratioCNplant
  REAL(r_2), DIMENSION(mvtype,msoil) :: ratioCNsoil
  REAL(r_2), DIMENSION(ms)        :: depthsoila,depthsoilb
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
  DO i=1,3
    READ(101,*) 
  ENDDO
  
  DO nv=1,mvtype
    READ(101,*) nv0,casabiome%ivt2(nv)
!     PRINT *, nv,nv0,casabiome%ivt2(nv)
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv1,casabiome%kroot(nv),casabiome%rootdepth(nv),      &
                casabiome%kuptake(nv),casabiome%krootlen(nv),         &
                casabiome%kminN(nv), casabiome%kuplabP(nv),           &
                xfherbivore(nv),leafage(nv),woodage(nv),frootage(nv), &
                metage(nv),strage(nv),cwdage(nv),  &
                micage(nv),slowage(nv),passage(nv) 
!     PRINT *, 'nv1',nv,nv1
  ENDDO  

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv2, &
                casabiome%fracnpptoP(nv,leaf),casabiome%fracnpptoP(nv,wood), &
                casabiome%fracnpptoP(nv,froot),casabiome%rmplant(nv,leaf),   &
                casabiome%rmplant(nv,wood),casabiome%rmplant(nv,froot)
!     PRINT *, 'nv2', nv2
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
         ratioCNsoil(nv,mic),ratioCNsoil(nv,slow),ratioCNsoil(nv,pass)
!,                            &
!         xfherbivore(nv),casabiome%ratiofrootleaf(nv),                   &
!         casabiome%glaimax(nv)
!     PRINT *, 'nv22',nv2
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv3, cleaf(nv),cwood(nv),cfroot(nv),cmet(nv),   &
                cstr(nv),ccwd(nv), cmic(nv), cslow(nv),cpass(nv)
!     PRINT *, 'nv3',nv3
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv4, &
         phen%TKshed(nv),xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv), &
         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
!     PRINT *, 'nv4',nv4
  ENDDO
!  READ(101,*)
!  READ(101,*)
!  DO nv=1,mvtype
!    READ(101,*) nv5, &
!         xxkleafcoldmax(nv),casabiome%xkleafcoldexp(nv),   &
!         xxkleafdrymax(nv),casabiome%xkleafdryexp(nv)
!  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv6, &
      casabiome%ratioNCplantmin(nv,leaf),casabiome%ratioNCplantmax(nv,leaf), &
      casabiome%ratioNCplantmin(nv,wood),casabiome%ratioNCplantmax(nv,wood), &
      casabiome%ratioNCplantmin(nv,froot),casabiome%ratioNCplantmax(nv,froot), &
      xfNminloss(nv), xfNminleach(nv),xnfixrate(nv)
!     PRINT *, 'nv6',nv6
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv7,nleaf(nv),nwood(nv),nfroot(nv), &
                nmet(nv),nstr(nv), ncwd(nv), &
                nmic(nv),nslow(nv),npass(nv),xnsoilmin(nv)
!     PRINT *, 'nv7',nv7
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv8,xratioNPleafmin,xratioNPleafmax,      &
         xratioNPwoodmin,xratioNPwoodmax,                      &
         xratioNPfrootmin,xratioNPfrootmax,                    &
         casabiome%ftransPPtoL(nv,leaf), casabiome%ftransPPtoL(nv,wood), &
         casabiome%ftransPPtoL(nv,froot)
    casabiome%ratioPcplantmin(nv,leaf)  = 1.0/(xratioNPleafmin*ratioCNplant(nv,leaf))
    casabiome%ratioPcplantmax(nv,leaf)  = 1.0/(xratioNPleafmax*ratioCNplant(nv,leaf))
    casabiome%ratioPcplantmin(nv,wood)  = 1.0/(xratioNPwoodmin*ratioCNplant(nv,wood))
    casabiome%ratioPcplantmax(nv,wood)  = 1.0/(xratioNPwoodmax*ratioCNplant(nv,wood))
    casabiome%ratioPcplantmin(nv,froot) = 1.0/(xratioNPfrootmin*ratioCNplant(nv,froot))
    casabiome%ratioPcplantmax(nv,froot) = 1.0/(xratioNPfrootmax*ratioCNplant(nv,froot))
!     PRINT *, 'nv8',nv8
  ENDDO

  READ(101,*)
  READ(101,*)
  DO iso=1,mso
    READ(101,*) nv9,xkmlabp(iso),xpsorbmax(iso),xfPleach(iso), &
                ratioNPsoil(iso,mic),ratioNPsoil(iso,slow),ratioNPsoil(iso,pass)
!     PRINT *, 'nv9',nv9
  ENDDO

  READ(101,*)
  READ(101,*)
  DO nv=1,mvtype
    READ(101,*) nv10, &
         xpleaf(nv),xpwood(nv),xpfroot(nv),xpmet(nv),xpstr(nv),xpcwd(nv), &
         xpmic(nv),xpslow(nv),xppass(nv),xplab(nv),xpsorb(nv),xpocc(nv)
!     PRINT *, 'nv10',nv10
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
    casabiome%sla(nv)             = 0.025 * (leafage(nv)**(-0.5)) ! see eqn A1 of Arora and Boer, GCB, 2005
!    casabiome%fherbivore(nv)      = deltcasa*xfherbivore(nv)
    casabiome%fraclabile(nv,leaf) = deltcasa*0.6    !1/day   !what is that for????
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
!    casabiome%kuplabp(nv)         = xkuplabp(nv)
    casabiome%rmplant(nv,:)       = casabiome%rmplant(nv,:)*deltcasa 
  ENDDO

  PRINT *, 'casabiome%ivt2 = ', casabiome%ivt2

  DO npt = 1, mp
    iv1=veg%iveg(npt)
    iso=casamet%isorder(npt)
    ! The following to be commented out when coupled to CABLE
!    veg%froot(npt,:) =fracroot(iv1,:)
!    PRINT *, 'npt,iv1,iso ', npt,iv1, iso
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

    casaflux%fNminloss(npt)   = xfNminloss(iv1) 
    ! comment out by ypw 12/07/2009
    casaflux%fNminleach(npt)  = 10.0*xfNminleach(iv1) * deltcasa
!    casaflux%fNminleach(npt)  = xfNminleach(iv1) 
    casapool%nplant(npt,leaf) = nleaf(iv1)
    casapool%nplant(npt,froot)= nfroot(iv1)
    casapool%nlitter(npt,metb) = nmet(iv1)
!    casapool%nlitter(npt,str) = nstr(iv1)
    casapool%nlitter(npt,str) = cstr(iv1)*ratioNCstrfix
    casapool%nsoil(npt,mic)   = nmic(iv1)
    casapool%nsoil(npt,slow)  = nslow(iv1)
    casapool%nsoil(npt,pass)  = npass(iv1) 
    casapool%nsoilmin(npt)    = xnsoilmin(iv1) 
    casapool%pplant(npt,leaf) = xpleaf(iv1)
    casapool%pplant(npt,froot)= xpfroot(iv1) 
    casapool%plitter(npt,metb) = xpmet(iv1)
!    casapool%plitter(npt,str) = xpstr(iv1)
    casapool%plitter(npt,str) = cstr(iv1)*ratioPCstrfix
    casapool%psoil(npt,mic)   = xpmic(iv1)
    casapool%psoil(npt,slow)  = xpslow(iv1)
    casapool%psoil(npt,pass)  = xppass(iv1)
    casapool%psoillab(npt)    = xplab(iv1)
    casapool%psoilsorb(npt)   = xpsorb(iv1)
    casapool%psoilocc(npt)    = xpocc(iv1)
    casaflux%kmlabp(npt)      = xkmlabp(iso)
    casaflux%psorbmax(npt)    = xpsorbmax(iso)
    casaflux%fpleach(npt)     = xfPleach(iso)
    casaflux%Nminfix(npt)     = xnfixrate(iv1)/365.0  

    casapool%rationcplant(npt,:)  = 1.0/ratioCNplant(iv1,:)
    casapool%ratiopcplant(npt,:)  = casabiome%ratioPcplantmax(iv1,:)
    casapool%rationclitter(npt,:) = casapool%nlitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratiopclitter(npt,:) = casapool%plitter(npt,:)/(casapool%clitter(npt,:)+1.0e-10)
    casapool%ratioNCsoil(npt,:)   = 1.0/ratioCNsoil(iv1,:)
    casapool%ratioPCsoil(npt,:)   = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iso,:))


!    casapool%rationcplant(npt,:) = casapool%nplant(npt,:)/(casapool%cplant(npt,:)+1.0e-10)
!    casapool%ratiopcplant(npt,:) = casapool%pplant(npt,:)/(casapool%cplant(npt,:)+1.0e-10)
!    casapool%ratioNCsoil(npt,:)  = 1.0/ratioCNsoil(iv1,:)
!    casapool%ratioPCsoil(npt,:)  = 1.0/(ratioCNsoil(iv1,:)*ratioNPsoil(iv1,:))
  ENDDO

  IF (icycle==1) THEN
    casapool%nplant(:,:)  = casapool%cplant(:,:) * casapool%rationcplant(:,:)
  ELSE
    casapool%Nsoil(:,:)   = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)   = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
    casapool%psoilsorb(:) = casaflux%psorbmax(:) * casapool%psoillab(:) &
                            /(casaflux%kmlabp(:)+casapool%psoillab(:))
  ENDIF 
      
!  DO npt=1,mp
!    IF (veg%iveg(npt)==12) PRINT *, npt, veg%iveg(npt), &
!         casapool%Psoil(npt,:),casapool%psoilsorb(npt), &
!         casaflux%psorbmax(npt),casapool%psoillab(npt),casaflux%kmlabp(npt)
!  ENDDO

END SUBROUTINE casa_readbiome

SUBROUTINE casa_readphen(veg,casamet,phen)
! added mvtype (=mvt) and mstype to define_dimensions (BP sep2010)
!SUBROUTINE casa_readphen(mvt,veg,casamet,phen)
  ! read in the tabulated modis-derived leaf phenology data
  ! for latitude bands of 79.75 to -55.25
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
!  INTEGER(i_d),              INTENT(IN)    :: mvt
  TYPE (veg_parameter_type), INTENT(IN)    :: veg  ! vegetation parameters
  TYPE (casa_met),           INTENT(IN)    :: casamet
  TYPE (phen_variable),      INTENT(INOUT) :: phen

  ! local variables
  INTEGER np,nx,ilat
  INTEGER, DIMENSION(271,mvtype) :: greenup, fall,  phendoy1
  INTEGER, DIMENSION(10)        :: greenupx,fallx,xphendoy1
  INTEGER, DIMENSION(10)        :: ivtx(10)
  REAL(r_2), DIMENSION(271)     :: xlat

  ! initilize for evergreen PFTs
  greenup(:,:) = -50
  fall(:,:)    = 367
  phendoy1(:,:)= 2

  OPEN(101,file=casafile%phen)
  READ(101,*)
  READ(101,*) (ivtx(nx),nx=1,10)   ! fixed at 10, as only 10 of 17 IGBP PFT
                                   ! have seasonal leaf phenology
  DO ilat=271,1,-1
    READ(101,*) xlat(ilat),(greenupx(nx),nx=1,10), &
                (fallx(nx),nx=1,10),(xphendoy1(nx),nx=1,10)
    DO nx=1,10
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
    phen%doyphase(np,4) = phen%doyphase(np,3) +14   ! DOY for minimal LAI season
    IF (phen%doyphase(np,2) > 365) phen%doyphase(np,2)=phen%doyphase(np,2)-365
    IF (phen%doyphase(np,4) > 365) phen%doyphase(np,4)=phen%doyphase(np,4)-365
  ENDDO

END SUBROUTINE casa_readphen

SUBROUTINE casa_readpoint(veg,soil,casaflux,casamet,rad)
! added mvtype (=mvt) and mstype to define_dimensions (BP sep2010)
!SUBROUTINE casa_readpoint(mvt,veg,soil,casaflux,casamet,patch,rad)
  USE define_dimensions
  USE define_types
  USE io_variables, ONLY: patch
  USE casaparm
  USE casadimension
  USE casavariable
  IMPLICIT NONE
!  INTEGER(i_d),               INTENT(IN)    :: mvt
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (radiation_type),      INTENT(IN)    :: rad

  ! local variables
  INTEGER :: np,nland
  REAL(r_2) :: annNdep,annNfix,annPwea,annPdust
  REAL(r_2) :: annNfert,annPfert   ! not really used yet
  INTEGER, DIMENSION(mp) :: vtypex,stypex
  INTEGER :: nlandx,ivtigbp, inPatch
     
  OPEN(101,file=casafile%cnppoint,FORM='FORMATTED')
  READ(101,*) inPatch
  PRINT * ,'Within casa_readpoint, mp = ', mp
  PRINT * ,'Input file has ', inPatch, ' patches.'

  np = 0
  DO nland=1,inPatch
    np = np + 1
    READ(101,*) &
        nlandx,ivtigbp,stypex(np),casamet%isorder(np), &
               casamet%lat(np),casamet%lon(np),casamet%areacell(np), &
               annNfix,annNdep,annNfert,annPwea,annPdust,annPfert
!    PRINT * , nlandx,ivtigbp,stypex(np),veg%iveg(np),soil%isoilm(np), &
!              patch(np)%frac,patch(np)%latitude,patch(np)%longitude, &
!              casamet%lat(np),casamet%lon(np)
    IF (ivtigbp == 0) ivtigbp = iceland

    IF (ABS(casamet%lat(np) - patch(np)%latitude) < 0.1 .AND. &
        ABS(casamet%lon(np) - patch(np)%longitude) < 0.1) THEN
      IF (ivtigbp /= veg%iveg(np) .OR. stypex(np) /= soil%isoilm(np)) THEN
        PRINT * ,'Check why iveg, isoil do not match'
        STOP
      ELSE
        casaflux%Nmindep(np) = annNdep/365.0
        casaflux%Nminfix(np) = annNfix/365.0
        casaflux%Pdep(np)    = annPdust/365.0     ! gP/m2/day
        casaflux%Pwea(np)    = annPwea/365.0      ! gP/m2/day
        IF (mvtype==17) THEN
          vtypex(np)  = ivtigbp  ! for running IGBP veg type only
        END IF 
      END IF
    ELSE
      PRINT * ,'Check why lat, lon do not match'
      STOP
    END IF

    IF (veg%iveg(np)==12 .OR. veg%iveg(np)==14) casaflux%Pdep(np)= &
       casaflux%Pdep(np)+0.7/365.0    ! P fertilizer =13 Mt P globally in 1994

  ENDDO 
  CLOSE(101)

END SUBROUTINE casa_readpoint

SUBROUTINE casa_init(casapool,casabal,veg)
! mst not used (BP sep2010)
!! for first time reading file *_1220.csv  (BP may2010)
!SUBROUTINE casa_init(mst,casapool,casabal,veg)
!!SUBROUTINE casa_init(mst,casapool,casabal)
!! end addition (BP may2010)
!  initialize some values in phenology parameters and leaf growth phase
  USE define_dimensions 
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
! for first time reading file *_1220.csv  (BP may2010)
  USE define_types
  USE io_variables, ONLY: landpt, patch
! end addition (BP may2010)
  IMPLICIT NONE
!  INTEGER(i_d),        INTENT(IN)    :: mst
  TYPE (casa_pool),    INTENT(INOUT) :: casapool
  TYPE (casa_balance), INTENT(INOUT) :: casabal
! for first time reading file *_1220.csv  (BP may2010)
  TYPE (veg_parameter_type), INTENT(IN) :: veg
  REAL(r_2) :: clabile,cplant(3),clitter(3),csoil(3)
  REAL(r_2) :: nplant(3),nlitter(3),nsoil(3),nsoilmin,pplant(3)
  REAL(r_2) :: plitter(3),psoil(3),psoillab,psoilsorb,psoilocc
! end addition (BP may2010)

  ! local variables
  INTEGER   :: np,npt,npz
  INTEGER   :: nyearz,ivtz,istz,isoz
  REAL(r_2) :: latz,lonz,areacellz,glaiz,slaz

  PRINT *, 'initial pool file: ',TRIM(casafile%cnpipool)
  PRINT *, 'initcasa,mp ', initcasa,mp
  !phen%phase = 2
  IF (initcasa==1) THEN
    OPEN(99,file=casafile%cnpipool)
!! for first time reading file *_1220.csv  (BP may2010)
!    PRINT *, 'Checking mland within casa_init = ', mland
!    DO npt =1, mland
!      READ(99,*) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz, &
!                   slaz,clabile,cplant(:),clitter(:),csoil(:),        &
!                   nplant(:),nlitter(:),nsoil(:),nsoilmin,pplant(:),  &
!                   plitter(:),psoil(:),psoillab,psoilsorb,psoilocc
!      ! check veg type and lat, lon
!      IF (ivtz /= veg%iveg(landpt(npt)%cstart)) THEN
!        PRINT *, 'npt,ivtz,veg%iveg: ', npt,ivtz,veg%iveg(landpt(npt)%cstart)
!        PRINT *, '         lat,lon : ', latz,lonz
!      ENDIF
!      IF ( ABS(latz - patch(landpt(npt)%cstart)%latitude) > 0.1  .OR. &
!           ABS(lonz - patch(landpt(npt)%cstart)%longitude)> 0.1 ) THEN
!        PRINT *, 'npt,latz,lonz = ', npt,latz,lonz
!        PRINT *, 'patch,lat,lon = ', landpt(npt)%cstart, &
!                                     patch(landpt(npt)%cstart)%latitude, &
!                                     patch(landpt(npt)%cstart)%longitude
!      ENDIF
!      ! put landpoint data into patch
!      DO np = landpt(npt)%cstart, landpt(npt)%cend
!!       IF (veg%iveg(landpt(npt)%cstart) == 11 .OR. &
!!           veg%iveg(landpt(npt)%cstart) == 13 .OR. &
!!           veg%iveg(landpt(npt)%cstart) == 15 .OR. &
!!           veg%iveg(landpt(npt)%cstart) == 17) THEN
!        casapool%clabile  (np)   = clabile
!        casapool%cplant   (np,:) = cplant(:)
!        casapool%clitter  (np,:) = clitter(:)
!        casapool%csoil    (np,:) = csoil(:)
!        casapool%nplant   (np,:) = nplant(:)
!        casapool%nlitter  (np,:) = nlitter(:)
!        casapool%nsoil    (np,:) = nsoil(:)
!        casapool%nsoilmin (np)   = nsoilmin
!        casapool%pplant   (np,:) = pplant(:)
!        casapool%plitter  (np,:) = plitter(:)
!        casapool%psoil    (np,:) = psoil(:)
!        casapool%psoillab (np)   = psoillab
!        casapool%psoilsorb(np)   = psoilsorb
!        casapool%psoilocc (np)   = psoilocc
!!       ELSE
!
!!       ENDIF
!      ENDDO
!    ENDDO

    DO npt =1, mp
      SELECT CASE(icycle)
      CASE(1)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz, &
                   slaz, casapool%clabile(npt),casapool%cplant(npt,:),  &
                   casapool%clitter(npt,:),casapool%csoil(npt,:)
      CASE(2)
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz, &
                   slaz,casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt)
      CASE(3)
!        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz
        READ(99,92) nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz, &
                   slaz,casapool%clabile(npt),casapool%cplant(npt,:),   &
                   casapool%clitter(npt,:),casapool%csoil(npt,:),       &
                   casapool%nplant(npt,:),casapool%nlitter(npt,:),      &
                   casapool%nsoil(npt,:),casapool%nsoilmin(npt),        &
                   casapool%pplant(npt,:),casapool%plitter(npt,:),      &
                   casapool%psoil(npt,:),casapool%psoillab(npt),        &
                   casapool%psoilsorb(npt),casapool%psoilocc(npt)
      END SELECT 
!      PRINT *, 'npt,nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz:'
!      PRINT *, npt,nyearz,npz,ivtz,istz,isoz,latz,lonz,areacellz,glaiz
    ENDDO
    CLOSE(99)
  ENDIF 
92  format(5(i6,3x),100(f15.6,3x))

  ! reset labile C pool
  casapool%clabile    = 0.0    
  ! check pool sizes
  casapool%cplant     = MAX(0.0,casapool%cplant)
  casapool%clitter    = MAX(0.0,casapool%clitter)
  casapool%csoil      = MAX(0.0,casapool%csoil)
  casabal%cplantlast  = casapool%cplant
  casabal%clitterlast = casapool%clitter
  casabal%csoillast   = casapool%csoil
  casabal%clabilelast = casapool%clabile
  casabal%sumcbal     = 0.0

  IF (icycle==1) THEN
    casapool%nplant(:,:) = casapool%cplant(:,:) * casapool%rationcplant(:,:)
    casapool%Nsoil(:,:)  = casapool%ratioNCsoil(:,:) * casapool%Csoil(:,:)
    casapool%Psoil(:,:)  = casapool%ratioPCsoil(:,:) * casapool%Csoil(:,:)
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
  ENDIF 

  IF (icycle >2) THEN
    casapool%pplant       = MAX(1.0e-7,casapool%pplant)
    casapool%plitter      = MAX(1.0e-7,casapool%plitter)
    casapool%psoil        = MAX(1.0e-7,casapool%psoil)
    casapool%Psoillab     = MAX(2.0,casapool%psoillab)
    casapool%psoilsorb    = MAX(10.0,casapool%psoilsorb)
    casapool%psoilocc     = MAX(50.0,casapool%psoilocc)
    casabal%pplantlast    = casapool%pplant
    casabal%plitterlast   = casapool%plitter
    casabal%psoillast     = casapool%psoil       
    casabal%psoillablast  = casapool%psoillab
    casabal%psoilsorblast = casapool%psoilsorb
    casabal%psoilocclast  = casapool%psoilocc
    casabal%sumpbal       = 0.0
  EndIF 

END SUBROUTINE casa_init


SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  INTEGER(i_d),               INTENT(IN)    :: ktau
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

!  WRITE(*,91) nyear,cplantsum,clittersum,csoilsum 
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
        casabiome%sla(veg%iveg(npt)), casapool%clabile(npt),    &
        casapool%cplant(npt,:),casapool%clitter(npt,:),casapool%csoil(npt,:), &
        casapool%nplant(npt,:),casapool%nlitter(npt,:),casapool%nsoil(npt,:), &
        casapool%nsoilmin(npt),casapool%pplant(npt,:),          &
        casapool%plitter(npt,:), casapool%psoil(npt,:),         &
        casapool%psoillab(npt),casapool%psoilsorb(npt),casapool%psoilocc(npt), &
        casabal%sumcbal(npt),casabal%sumnbal(npt),casabal%sumpbal(npt)
  ENDDO

  CLOSE(nout)

91    format(i6,100(g9.3,2x))
92    format(5(i6,',',2x),100(f15.6,',',2x))
END SUBROUTINE casa_poolout

! casa_fluxout output data for Julie Tang; comment out (BP apr2010)
!SUBROUTINE casa_fluxout(myear,clitterinput,csoilinput)
!  USE define_dimensions
!  USE define_types
!  USE cableDeclare, ONLY: veg, soil
!  USE casadimension
!  USE casaparm
!  USE casavariable
!  USE phenvariable
!  USE casaDeclare
!  IMPLICIT NONE
!  INTEGER(i_d), INTENT(IN)      :: myear
!  REAL(r_2),    INTENT(IN) :: clitterinput(mp,3),csoilinput(mp,3)
!
!  ! local variables
!  INTEGER  npt,nout
!  REAL(r_2) xyear
!
!  nout=104
!  xyear=1.0/FLOAT(myear)
!  casabal%FCnppyear=casabal%FCnppyear * xyear
!  casabal%FCrsyear=casabal%FCrsyear * xyear
!  casabal%FCneeyear=casabal%FCneeyear * xyear
!  casabal%FNdepyear=casabal%FNdepyear * xyear
!  casabal%FNfixyear=casabal%FNfixyear * xyear
!  casabal%FNsnetyear=casabal%FNsnetyear * xyear
!  casabal%FNupyear=casabal%FNupyear * xyear
!  casabal%FNleachyear=casabal%FNleachyear * xyear
!  casabal%FNlossyear=casabal%FNlossyear * xyear
!  casabal%FPweayear=casabal%FPweayear * xyear
!  casabal%FPdustyear=casabal%FPdustyear * xyear
!  casabal%FPsnetyear=casabal%FPsnetyear * xyear
!  casabal%FPupyear=casabal%FPupyear * xyear
!  casabal%FPleachyear=casabal%FPleachyear * xyear
!  casabal%FPlossyear=casabal%FPlossyear * xyear
!  clitterinput = clitterinput * xyear
!  csoilinput   = csoilinput   * xyear
!
!  OPEN(nout,file=casafile%cnpflux)
!    DO npt =1,mp
!      SELECT CASE(icycle)
!      CASE(1)
!        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
!            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
!            casamet%areacell(npt)*(1.0e-9),casabal%Fcnppyear(npt),  &
!            casabal%Fcrsyear(npt),casabal%Fcneeyear(npt),           &
!            clitterinput(npt,:),csoilinput(npt,:)
!
!      CASE(2)
!        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),    &
!            casamet%isorder(npt),casamet%lat(npt),casamet%lon(npt), &
!            casamet%areacell(npt)*(1.0e-9), &
!        casabal%FCnppyear(npt),casabal%FCrsyear(npt),   casabal%FCneeyear(npt), &
!        clitterinput(npt,:),csoilinput(npt,:), &
!        casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt), &
!        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt)
!
!      CASE(3)
!        WRITE(nout,92) myear,npt,veg%iveg(npt),soil%isoilm(npt),casamet%isorder(npt), &
!        casamet%lat(npt),casamet%lon(npt),casamet%areacell(npt)*(1.0e-9), &
!        casabal%FCnppyear(npt),casabal%FCrsyear(npt),   casabal%FCneeyear(npt), &
!        clitterinput(npt,:),csoilinput(npt,:), &
!        casabal%FNdepyear(npt),casabal%FNfixyear(npt),  casabal%FNsnetyear(npt), &
!        casabal%FNupyear(npt), casabal%FNleachyear(npt),casabal%FNlossyear(npt), &
!        casabal%FPweayear(npt),casabal%FPdustyear(npt), casabal%FPsnetyear(npt), &
!        casabal%FPupyear(npt), casabal%FPleachyear(npt),casabal%FPlossyear(npt)
!
!      END SELECT 
!    ENDDO
!
!  CLOSE(nout)
!92    format(5(i6,',',2x),100(f15.6,',',2x))
!END SUBROUTINE casa_fluxout

! clitterinput and csoilinput are for Julie Tang; comment out (BP apr2010)
!SUBROUTINE casa_cnpflux(clitterinput,csoilinput)
SUBROUTINE casa_cnpflux(casaflux,casabal)
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  TYPE (casa_flux),    INTENT(INOUT) :: casaflux
  TYPE (casa_balance), INTENT(INOUT) :: casabal
!  REAL(r_2), INTENT(INOUT) :: clitterinput(mp,3),csoilinput(mp,3)
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

SUBROUTINE biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
  USE define_dimensions
  USE define_types
  USE casadimension
  USE casa_cnp
  IMPLICIT NONE
  INTEGER(i_d), INTENT(IN)    :: ktau
  REAL(r_1),    INTENT(IN)    :: dels
  INTEGER(i_d), INTENT(IN)    :: idoy
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL(r_2),    DIMENSION(mp) :: xnplimit,xNPuptake
  REAL(r_2),    DIMENSION(mp) :: xklitter,xksoil,xkNlimiting
  REAL(r_2),    DIMENSION(mp) :: xkleafcold,xkleafdry,xkleaf
  INTEGER  npt

!  PRINT *, 'Within biogeochem, mp = ', mp
  xKNlimiting = 1.0
  call phenology(idoy,phen)
!  PRINT *, 'calling avgsoil', mp
  call avgsoil(veg,soil,casamet)
!  PRINT *, 'calling casa_rplant', mp
  call casa_rplant(veg,casabiome,casapool,casaflux,casamet)

!  PRINT *, 'calling allocation', mp
  call casa_allocation(veg,soil,casabiome,casaflux,casamet,phen)

!  PRINT *, 'calling xrateplant',mp
  call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                       casamet,phen)

!  PRINT *, 'calling coeffplant',mp
  call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                       casaflux,casamet)

!  PRINT *, 'calling xnp',mp,idoy
  call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)

  IF (icycle>1) THEN
!    PRINT *, 'calling nuptake'
    call casa_nuptake(veg,casabiome,casapool,casaflux,casamet)
    IF (icycle >2) call casa_puptake(veg,casabiome,casapool,casaflux,casamet)
  ENDIF 

!  PRINT *, 'calling delplant',mp
  call casa_delplant(veg,casabiome,casapool,casaflux,casamet)

  IF (icycle>1) call casa_xkN(xkNlimiting,casapool,casaflux,casamet)

!  PRINT *, 'calling xratesoil',mp
  call casa_xratesoil(xklitter,xksoil,veg,soil,casamet)

!  PRINT *, 'calling coeffsoil'
  call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

!  PRINT *, 'calling delsoil'
  call casa_delsoil(veg,casapool,casaflux,casamet)
!  PRINT *, 'nsoilmin,xNPuptake xnp,xkn', casapool%Nsoilmin(1),xKNlimiting(1), &
!           xNPuptake(1),casaflux%Nsmin(1),casaflux%Nsimm(1), &
!           casaflux%Nsnet(1),casapool%Nsoilmin(1)*deltpool

!  PRINT *, 'calling cnpcycle'
  call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)

  IF (icycle==1) call casa_ndummy(casapool)

!  PRINT *, 'calling cnpbal'
  call casa_cnpbal(casapool,casaflux,casabal)

  ! Trick added (BP 28jun2010)
  casapool%Nsoilmin(:) = MAX(0.5, casapool%Nsoilmin(:))
  casapool%Psoillab(:) = MAX(0.1, casapool%Psoillab(:))

END SUBROUTINE biogeochem

