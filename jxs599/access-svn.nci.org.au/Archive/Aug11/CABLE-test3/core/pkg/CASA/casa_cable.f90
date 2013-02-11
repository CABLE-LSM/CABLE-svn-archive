! casa_cable.f90
!
! Model development by YingPing Wang, CSIRO Marine and Atmospheric Research.
! Coupling to Mk3L by Bernard Pak,    CSIRO Marine and Atmospheric Research.
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! the following routines are used when "casacnp" is coupled to "cable"
!   bgcdriver - interface between casacnp and cable
!   sumcflux  - accumulating carbon fluxes
!   spincasacnp - for offline spinup only, NOT USED with Mk3L

SUBROUTINE bgcdriver(ktau,kstart,kend,dels,met,ssoil,canopy,veg,soil, &
                     casabiome,casapool,casaflux,casamet,casabal,phen)

  USE define_dimensions
  USE define_types
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE

  INTEGER(i_d),      INTENT(IN) :: ktau ! integration step number
  INTEGER(i_d),      INTENT(IN) :: kstart ! starting value of ktau
  INTEGER(i_d),      INTENT(IN) :: kend ! total # timesteps in run
  REAL(r_1),         INTENT(IN) :: dels ! time setp size (s)
  TYPE (met_type), INTENT(INOUT)       :: met  ! met input variables
  TYPE (soil_snow_type), INTENT(INOUT) :: ssoil ! soil and snow variables
  TYPE (canopy_type), INTENT(INOUT) :: canopy ! vegetation variables
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),          INTENT(INOUT) :: casabiome
  TYPE (casa_pool),           INTENT(INOUT) :: casapool
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux
  TYPE (casa_met),            INTENT(INOUT) :: casamet
  TYPE (casa_balance),        INTENT(INOUT) :: casabal
  TYPE (phen_variable),       INTENT(INOUT) :: phen

  ! local variables
  INTEGER(i_d)                  :: idoy ! day of year (1-365)
  INTEGER(i_d)                  :: ktauday

    ktauday=int(24.0*3600.0/dels)
    idoy = mod(ktau/ktauday,365)
    if(idoy==0) idoy=365
    phen%phase = 2

    if(ktau == kstart) then
       casamet%tairk  = 0.0
       casamet%tsoil  = 0.0
       casamet%moist  = 0.0
       casaflux%cgpp  = 0.0
! add initializations (BP jul2010)
       casaflux%Crsoil   = 0.0
       casaflux%crgplant = 0.0
       casaflux%crmplant = 0.0
!       casaflux%crmplant(:,leaf) = 0.0
! end changes (BP jul2010)
    endif
    if(mod(ktau,ktauday)==1) then
       casamet%tairk  = met%tk
       casamet%tsoil = ssoil%tgg
       casamet%moist = ssoil%wb
       casaflux%cgpp = (-canopy%fpn+canopy%frday)*dels
       casaflux%crmplant(:,leaf) = canopy%frday*dels
    else
       casamet%tairk  =casamet%tairk + met%tk
       casamet%tsoil = casamet%tsoil + ssoil%tgg
       casamet%moist = casamet%moist + ssoil%wb
       casaflux%cgpp = casaflux%cgpp + (-canopy%fpn+canopy%frday)*dels
       casaflux%crmplant(:,leaf) =casaflux%crmplant(:,leaf) + canopy%frday*dels
    endif

    if(mod((ktau-kstart+1),ktauday)==0) then
       casamet%tairk  =casamet%tairk/FLOAT(ktauday)
       casamet%tsoil=casamet%tsoil/FLOAT(ktauday)
       casamet%moist=casamet%moist/FLOAT(ktauday)
       WRITE(57,*) 'calling biogeochem',ktau,idoy,mp
       CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
    endif
! BP moved this to interface so that it will not just output one latitude band
! for Mk3L run.
!    if(ktau==kend) then
!!       PRINT *, 'Before casa_poolout'
!       call casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
!                         casabal,phen)
!!       PRINT *, 'After casa_poolout'
!    endif
111 format(i6,100(f12.5,2x))
END SUBROUTINE bgcdriver

SUBROUTINE sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
                    soil, ssoil, sum_flux, veg, met, casaflux)
! added mvtype and mstype to define_dimensions (BP sep2010)
!SUBROUTINE sumcflux(ktau, kstart, kend, dels, mvtype, mstype, bgc, canopy,  &
!                    soil, ssoil, sum_flux, veg, met, casaflux)

  USE define_dimensions
  USE define_types
  USE carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  IMPLICIT NONE
  INTEGER(i_d), INTENT(IN)    :: ktau ! integration step number
  INTEGER(i_d), INTENT(IN)    :: kstart ! starting value of ktau
  INTEGER(i_d), INTENT(IN)    :: kend ! total # timesteps in run
!  INTEGER(i_d), INTENT(IN)    :: mvtype  ! Number of veg types
!  INTEGER(i_d), INTENT(IN)    :: mstype ! Number of soil types
  REAL(r_1),    INTENT(IN)    :: dels ! time setp size (s)
  TYPE (bgc_pool_type),       INTENT(INOUT) :: bgc
  TYPE (canopy_type),         INTENT(INOUT) :: canopy
  TYPE (soil_parameter_type), INTENT(INOUT) :: soil
  TYPE (soil_snow_type),      INTENT(INOUT) :: ssoil
  TYPE (sum_flux_type),       INTENT(INOUT) :: sum_flux
  TYPE (met_type),            INTENT(IN)    :: met    
  TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
  TYPE (casa_flux),           INTENT(INOUT) :: casaflux

    if(icycle<=0) then
       CALL plantcarb(veg,bgc,met,canopy)
       CALL soilcarb(soil, ssoil, veg, bgc, met, canopy)
       CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc)
! mvtype and mstype are now declared in define_dimensions (BP sep2010)
!       CALL soilcarb(soil, ssoil, veg, bgc, met, canopy, mstype)
!       CALL carbon_pl(dels, soil, ssoil, veg, canopy, bgc, mvtype)
    else
       canopy%frp(:) = (casaflux%crmplant(:,wood)+casaflux%crmplant(:,froot) &
                        +casaflux%crgplant(:))/86400.0
       canopy%frs(:) = casaflux%Crsoil(:)/86400.0
       canopy%frpw(:)= casaflux%crmplant(:,wood)/86400.0
       canopy%frpr(:)= casaflux%crmplant(:,froot)/86400.0
    endif
    if(ktau == kstart) then
       sum_flux%sumpn  = canopy%fpn*dels
       sum_flux%sumrd  = canopy%frday*dels
       sum_flux%dsumpn = canopy%fpn*dels
       sum_flux%dsumrd = canopy%frday*dels
       sum_flux%sumrpw = canopy%frpw*dels
       sum_flux%sumrpr = canopy%frpr*dels
       sum_flux%sumrp  = canopy%frp*dels
       sum_flux%dsumrp = canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = canopy%frs*dels
    else
       sum_flux%sumpn  = sum_flux%sumpn  + canopy%fpn*dels
       sum_flux%sumrd  = sum_flux%sumrd  + canopy%frday*dels
       sum_flux%dsumpn = sum_flux%dsumpn + canopy%fpn*dels
       sum_flux%dsumrd = sum_flux%dsumrd + canopy%frday*dels
       sum_flux%sumrpw = sum_flux%sumrpw + canopy%frpw*dels
       sum_flux%sumrpr = sum_flux%sumrpr + canopy%frpr*dels
       sum_flux%sumrp  = sum_flux%sumrp  + canopy%frp*dels
       sum_flux%dsumrp = sum_flux%dsumrp + canopy%frp*dels
    ! canopy%frs set in soilcarb
       sum_flux%sumrs = sum_flux%sumrs+canopy%frs*dels
    endif
    ! Set net ecosystem exchange after adjustments to frs:
    canopy%fnee = canopy%fpn + canopy%frs + canopy%frp
    canopy%fnpp = -1.0* canopy%fpn - canopy%frp
!    write(*,101) ktau,casaflux%Crsoil(:)
!    write(*,101) ktau,dels,veg%vlai,veg%vcmax*1.0e6,casaflux%cgpp,canopy%fpn*1.0e6/12.0,canopy%frp*1.0e6/12.0,canopy%frs*1.0e6/12.0,canopy%fnee*1.0e6/12.0
101  format(i6,2x,100(f12.5,2x))
    if(ktau==kend) then
       PRINT *, 'carbon fluxes'
       PRINT *, 'sumpn', sum_flux%sumpn
       PRINT *, 'sumrd', sum_flux%sumrd
       PRINT *, 'sumrp', sum_flux%sumrp
       PRINT *, 'sumrs', sum_flux%sumrs
       PRINT *, 'npp', sum_flux%sumpn+sum_flux%sumrp
       PRINT *, 'nee', sum_flux%sumpn+sum_flux%sumrp+sum_flux%sumrs
     !  PRINT *, 'carbon pools', leaf,wood,froot
     !  PRINT *,  casapool%cplant(1,2),casaflux%crmplant(1,wood),casaflux%Crsoil(1)
     !  PRINT *, 'respiration rate'
     !  PRINT *,  casabiome%rmplant(1,2)*365.0
    endif

END SUBROUTINE sumcflux


SUBROUTINE spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                       casaflux,casamet,casabal,phen)
  USE define_dimensions 
  USE define_types
  USE carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  REAL(r_1),    INTENT(IN)    :: dels
  INTEGER(i_d), INTENT(IN)    :: kstart
  INTEGER(i_d), INTENT(IN)    :: kend
  INTEGER(i_d), INTENT(IN)    :: mloop
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  REAL, DIMENSION(:,:,:), ALLOCATABLE   :: xtsoil,xmoist
  REAL, DIMENSION(:,:),   ALLOCATABLE   :: xcgpp,xcrmleaf
  REAL, DIMENSION(:,:),   ALLOCATABLE   :: xtairk
  INTEGER :: ktauday,nloop,idoy,nday,ktaux,ktauy
  REAL    :: xlai
    
  ktauday=int(24.0*3600.0/dels)
  nday=(kend-kstart+1)/ktauday
  PRINT *, 'nday mp ms ',nday,kend,kstart,mp,ms
  ALLOCATE(xtsoil(nday,mp,ms),xmoist(nday,mp,ms))
  ALLOCATE(xcgpp(nday,mp),xcrmleaf(nday,mp))
  ALLOCATE(xtairk(nday,mp))

  OPEN(111,file=casafile%cnpmet)
  DO idoy=1,nday
    READ(111,*) ktaux,xlai,xtairk(idoy,:),xtsoil(idoy,:,:),xmoist(idoy,:,:), &
                xcgpp(idoy,:),xcrmleaf(idoy,:)
  ENDDO
  close(111)
  nloop = 0 
  DO nloop=1,mloop
    DO idoy=1,nday
      ktauy=idoy*ktauday
      casamet%tairk(:)=xtairk(idoy,:)
      casamet%tsoil(:,:)=xtsoil(idoy,:,:)
      casamet%moist(:,:)=xmoist(idoy,:,:)
      casaflux%cgpp(:) = xcgpp(idoy,:)
      casaflux%crmplant(:,leaf) =xcrmleaf(idoy,:)
      call biogeochem(ktauy,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                      casamet,casabal,phen)
    ENDDO
    IF((nloop+10)>mloop) THEN
      WRITE(*,151), nloop, casapool%cplant,casapool%clitter,casapool%csoil
    ENDIF
  ENDDO
151 FORMAT(i6,100(f12.5,2x))
END SUBROUTINE spincasacnp

