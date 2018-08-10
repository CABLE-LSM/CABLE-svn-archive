MODULE landuse_constant
  USE cable_def_types_mod,  ONLY: mvtype,mstype,mland,ms,msn,nrb
  USE casadimension,        ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  IMPLICIT NONE
  integer, parameter        :: dp =selected_real_kind(8)
  !  
  INTEGER,   PARAMETER      :: mstate   = 12  ! number of land use states
  INTEGER,   PARAMETER      :: mvmax    = 21  ! UM setup require 27, worry later
  INTEGER,   PARAMETER      :: mharvw   = 5

  real,          PARAMETER                   :: fracmin  = 0.0
  REAL(kind=dp), PARAMETER, DIMENSION(mwood) :: fwoodprod    =(/0.3,0.4,0.4/)
  REAL(kind=dp), PARAMETER, DIMENSION(mwood) :: ratewoodprod =(/1.0,0.1,0.01/)
  REAL(kind=dp), PARAMETER, DIMENSION(mwood) :: fracwoodseed =(/0.4,0.3,0.3/)
  REAL(kind=dp), PARAMETER, DIMENSION(mwood) :: fracgrassseed=(/0.6,0.0,0.4/)
  REAL,          PARAMETER                   :: fseedling = 0.001
END MODULE landuse_constant

MODULE landuse_variable
  use landuse_constant
  IMPLICIT NONE
  character*120      :: fpft,fxpft,fxluh2cable
  TYPE landuse_type
    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: phase_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: patchfrac_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: clabile_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: cplant_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: clitter_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: csoil_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: cwoodprod_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nplant_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nlitter_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nsoil_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: nsoilmin_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nwoodprod_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: pplant_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: plitter_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: psoil_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoillab_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoilsorb_x
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoilocc_x
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: pwoodprod_x

    INTEGER,   DIMENSION(:,:),       ALLOCATABLE :: phase_y
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: patchfrac_y
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: clabile_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: cplant_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: clitter_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: csoil_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: cwoodprod_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nplant_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nlitter_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nsoil_y 
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: nsoilmin_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: nwoodprod_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: pplant_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: plitter_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: psoil_y
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoillab_y
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoilsorb_y
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: psoilocc_y
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: pwoodprod_y

  ! landuse data
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: primary
    REAL(kind=dp), DIMENSION(:,:),   ALLOCATABLE :: fharvw 
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: xluh2cable 
    REAL(kind=dp), DIMENSION(:,:,:), ALLOCATABLE :: atransit
  END TYPE landuse_type

  CONTAINS

  SUBROUTINE landuse_setup(mland,luc)
   use landuse_constant
   IMPLICIT NONE
   TYPE(landuse_type), INTENT(INOUT)  :: luc
   integer  mland
! Allocate temprary pools.
  ALLOCATE(luc%phase_x(mland,mvmax),            &
           luc%patchfrac_x(mland,mvmax),        &
           luc%cplant_x(mland,mvmax,mplant),    &
           luc%nplant_x(mland,mvmax,mplant),    &
           luc%pplant_x(mland,mvmax,mplant),    &
           luc%clitter_x(mland,mvmax,mlitter),  &
           luc%nlitter_x(mland,mvmax,mlitter),  &
           luc%plitter_x(mland,mvmax,mlitter),  &
           luc%csoil_x(mland,mvmax,msoil),      &
           luc%nsoil_x(mland,mvmax,msoil),      &
           luc%psoil_x(mland,mvmax,msoil),      &
           luc%clabile_x(mland,mvmax),          &
           luc%nsoilmin_x(mland,mvmax),         &
           luc%psoillab_x(mland,mvmax),         &
           luc%psoilsorb_x(mland,mvmax),        &
           luc%psoilocc_x(mland,mvmax),         &
           luc%cwoodprod_x(mland,mvmax,mwood),  &
           luc%nwoodprod_x(mland,mvmax,mwood),  &
           luc%pwoodprod_x(mland,mvmax,mwood),  &
           luc%phase_y(mland,mvmax),            &
           luc%patchfrac_y(mland,mvmax),        &
           luc%cplant_y(mland,mvmax,mplant),    &
           luc%nplant_y(mland,mvmax,mplant),    &
           luc%pplant_y(mland,mvmax,mplant),    &
           luc%clitter_y(mland,mvmax,mlitter),  &
           luc%nlitter_y(mland,mvmax,mlitter),  &
           luc%plitter_y(mland,mvmax,mlitter),  &
           luc%csoil_y(mland,mvmax,msoil),      &
           luc%nsoil_y(mland,mvmax,msoil),      &
           luc%psoil_y(mland,mvmax,msoil),      &
           luc%clabile_y(mland,mvmax),          &
           luc%nsoilmin_y(mland,mvmax),         &
           luc%psoillab_y(mland,mvmax),         &
           luc%psoilsorb_y(mland,mvmax),        &
           luc%psoilocc_y(mland,mvmax),         &
           luc%cwoodprod_y(mland,mvmax,mwood),  &
           luc%nwoodprod_y(mland,mvmax,mwood),  &
           luc%pwoodprod_y(mland,mvmax,mwood),  &
           luc%primary(mland,mvtype),           &
           luc%fharvw(mland,mharvw),            &
           luc%xluh2cable(mland,mvmax,mstate),  &
           luc%atransit(mland,mvmax,mvmax))

! Initialize temporary variables
           luc%phase_x   = -1;    luc%patchfrac_x=0.0
           luc%cplant_x  = 0.0;   luc%nplant_x   = 0.0;    luc%pplant_x   = 0.0
           luc%clitter_x = 0.0;   luc%nlitter_x  = 0.0;    luc%plitter_x  = 0.0
           luc%csoil_x   = 0.0;   luc%nsoil_x    = 0.0;    luc%psoil_x    = 0.0
           luc%clabile_x = 0.0;   luc%nsoilmin_x = 0.0;    luc%psoillab_x = 0.0
           luc%psoilsorb_x = 0.0; luc%psoilocc_x = 0.
           luc%cwoodprod_x=0.0;   luc%nwoodprod_x=0.0;     luc%pwoodprod_x=0.0

           luc%phase_y   = -1;    luc%patchfrac_y=0.0
           luc%cplant_y  = 0.0;   luc%nplant_y   = 0.0;    luc%pplant_y   = 0.0
           luc%clitter_y = 0.0;   luc%nlitter_y  = 0.0;    luc%plitter_y  = 0.0
           luc%csoil_y   = 0.0;   luc%nsoil_y    = 0.0;    luc%psoil_y    = 0.0
           luc%clabile_y = 0.0;   luc%nsoilmin_y = 0.0;    luc%psoillab_y = 0.0
           luc%psoilsorb_y = 0.0; luc%psoilocc_y = 0.0
           luc%cwoodprod_y=0.0;   luc%nwoodprod_y=0.0;     luc%pwoodprod_y=0.0

           luc%primary = 0.0;     luc%fharvw=0.0;          luc%xluh2cable=0.0;    luc%atransit=0.0
   END SUBROUTINE landuse_setup

   SUBROUTINE landuse_deallocate(luc)
   IMPLICIT NONE
   TYPE(landuse_type), INTENT(INOUT)  :: luc
! Allocate temprary pools.
   DEALLOCATE(luc%phase_x,           &
              luc%patchfrac_x,        &
              luc%cplant_x,          &
              luc%nplant_x,          &
              luc%pplant_x,          &
              luc%clitter_x,         &
              luc%nlitter_x,         &
              luc%plitter_x,         &
              luc%csoil_x,           &
              luc%nsoil_x,           &
              luc%psoil_x,           &
              luc%clabile_x,         &
              luc%nsoilmin_x,        &
              luc%psoillab_x,        &
              luc%psoilsorb_x,       &
              luc%psoilocc_x,        &
              luc%cwoodprod_x,       &
              luc%nwoodprod_x,       &
              luc%pwoodprod_x,       &
              luc%phase_y,           &
              luc%patchfrac_y,        &
              luc%cplant_y,          &
              luc%nplant_y,          &
              luc%pplant_y,          &
              luc%clitter_y,         &
              luc%nlitter_y,         &
              luc%plitter_y,         &
              luc%csoil_y,           &
              luc%nsoil_y,           &
              luc%psoil_y,           &
              luc%clabile_y,         &
              luc%nsoilmin_y,        &
              luc%psoillab_y,        &
              luc%psoilsorb_y,       &
              luc%psoilocc_y,        &
              luc%cwoodprod_y,       &
              luc%nwoodprod_y,       &
              luc%pwoodprod_y,       &
              luc%primary,           &
              luc%fharvw,            &
              luc%xluh2cable,        &
              luc%atransit)

   END SUBROUTINE landuse_deallocate

END MODULE landuse_variable

MODULE landuse_patch
 use landuse_constant
 implicit none

 TYPE landuse_mp
   ! generic patch properties
   integer,        dimension(:),        allocatable   :: iveg,isoil,soilorder,phase,isflag,meth      !(mp)
   real(kind=dp),  dimension(:),        allocatable   :: zse                                         !(ms)    
   real(kind=dp),  dimension(:),        allocatable   :: patchfrac,areacell,lai,slai                 !(mp)
   ! the following assigned in biophysical parameter lookup table (vegetation)
   real(kind=dp),  dimension(:),        allocatable   :: hc,xfang,dleaf,frac4,ejmax,vcmax,rp20,rpcoef, &
                                                         wai,shelrb,vegcf,extkn,tminvj,tmaxvj,vbeta, &
                                                         canst1,xalbnir                              !(mp)
   real(kind=dp),  dimension(:,:),      allocatable   :: albsoil                                     !(mp,nrb) 
   ! the following assigned using biophysical parameter lookup table (soil)
   real(kind=dp),  dimension(:),        allocatable   :: clay,sand,silt,ssat,sfc,swilt,bch,hyds,sucs,css,rhosoil !(mp)

   real(kind=dp),  dimension(:),        allocatable   :: rtsoil,wetfac,tss,runoff,rnof1,rnof2,cansto,       &
                                                         ghflux,sghflux,ga,dgdtg,fev,fes,fhs,wbtot0,osnowd0,trad,  &
                                                         za_uv,za_tq,ssdnn,snowd,snage,osnowd     !(mp)
   real(kind=dp),  dimension(:,:),      allocatable   :: albedo,albsoilsn                             !(mp,nrb)
   real(kind=dp),  dimension(:,:),      allocatable   :: tggsn,ssdn,smass,sdepth                      !(snow,mp)
   real(kind=dp),  dimension(:,:),      allocatable   :: tgg,wb,wbice,gammzz,froot                    !(ms,mp)
   ! biophysical variables
   ! biogeochemical variables
   real(kind=dp),  dimension(:),        allocatable   :: sumcbal,sumnbal,sumpbal                      !(mp)                   
   real(kind=dp),  dimension(:),        allocatable   :: clabile,nsoilmin,psoillab,psoilsorb,psoilocc
   real(kind=dp),  dimension(:,:),      allocatable   :: cplant,nplant,pplant
   real(kind=dp),  dimension(:,:),      allocatable   :: clitter,nlitter,plitter
   real(kind=dp),  dimension(:,:),      allocatable   :: csoil,nsoil,psoil
   real(kind=dp),  dimension(:,:),      allocatable   :: cwoodprod,nwoodprod,pwoodprod

 END TYPE landuse_mp
 
 CONTAINS

   SUBROUTINE landuse_allocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)
   integer    mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood
   TYPE(landuse_mp), INTENT(INOUT)  :: lucmp

     allocate(lucmp%iveg(mpx),lucmp%isoil(mpx),lucmp%soilorder(mpx),lucmp%phase(mpx),lucmp%meth(mpx),lucmp%isflag(mpx))
     allocate(lucmp%zse(ms))
     allocate(lucmp%patchfrac(mpx),lucmp%areacell(mpx),lucmp%lai(mpx),lucmp%slai(mpx))
     allocate(lucmp%hc(mpx),lucmp%xfang(mpx),lucmp%dleaf(mpx),lucmp%frac4(mpx),lucmp%ejmax(mpx),lucmp%vcmax(mpx),lucmp%rp20(mpx),lucmp%rpcoef(mpx), &
              lucmp%wai(mpx),lucmp%shelrb(mpx),lucmp%vegcf(mpx),lucmp%extkn(mpx),lucmp%tminvj(mpx),lucmp%tmaxvj(mpx),lucmp%vbeta(mpx), &
              lucmp%canst1(mpx),lucmp%xalbnir(mpx))                              
     allocate(lucmp%albsoil(mpx,nrb))
     allocate(lucmp%clay(mpx),lucmp%sand(mpx),lucmp%silt(mpx),lucmp%ssat(mpx),lucmp%sfc(mpx),lucmp%swilt(mpx), &
              lucmp%bch(mpx),lucmp%hyds(mpx),lucmp%sucs(mpx),lucmp%css(mpx),lucmp%rhosoil(mpx))

     allocate(lucmp%rtsoil(mpx),lucmp%wetfac(mpx),lucmp%tss(mpx),lucmp%runoff(mpx),lucmp%rnof1(mpx),lucmp%rnof2(mpx), &
              lucmp%cansto(mpx),lucmp%ghflux(mpx),lucmp%sghflux(mpx),lucmp%ga(mpx),                &
              lucmp%dgdtg(mpx),lucmp%fev(mpx),lucmp%fes(mpx),lucmp%fhs(mpx),lucmp%wbtot0(mpx),                       &
              lucmp%osnowd0(mpx),lucmp%trad(mpx),lucmp%za_uv(mpx),lucmp%za_tq(mpx),                  &
              lucmp%ssdnn(mpx),lucmp%snowd(mpx),lucmp%snage(mpx),lucmp%osnowd(mpx))
     allocate(lucmp%albedo(mpx,nrb),lucmp%albsoilsn(mpx,nrb))
     allocate(lucmp%tggsn(mpx,msn),lucmp%ssdn(mpx,msn),lucmp%smass(mpx,msn),lucmp%sdepth(mpx,msn))
     allocate(lucmp%tgg(mpx,ms),lucmp%wb(mpx,ms),lucmp%wbice(mpx,ms),lucmp%gammzz(mpx,ms),lucmp%froot(mpx,ms)) 
   
     allocate(lucmp%sumcbal(mpx),lucmp%sumnbal(mpx),lucmp%sumpbal(mpx))
     allocate(lucmp%clabile(mpx))
     allocate(lucmp%cplant(mpx,mplant),lucmp%nplant(mpx,mplant),lucmp%pplant(mpx,mplant))
     allocate(lucmp%clitter(mpx,mlitter),lucmp%nlitter(mpx,mlitter),lucmp%plitter(mpx,mlitter))
     allocate(lucmp%csoil(mpx,msoil),lucmp%nsoil(mpx,msoil),lucmp%psoil(mpx,msoil))
     allocate(lucmp%nsoilmin(mpx))
     allocate(lucmp%psoillab(mpx),lucmp%psoilsorb(mpx),lucmp%psoilocc(mpx))
     allocate(lucmp%cwoodprod(mpx,mwood),lucmp%nwoodprod(mpx,mwood),lucmp%pwoodprod(mpx,mwood))

     ! initialization
     lucmp%iveg=-1;lucmp%isoil=-1;lucmp%soilorder=-1;lucmp%phase=0;lucmp%meth=0;lucmp%isflag=0
     lucmp%zse=0.0
     lucmp%patchfrac=0.0;lucmp%areacell=0.0;lucmp%lai=0.0;lucmp%slai=0.0
     lucmp%hc=0.0;lucmp%xfang=0.0;lucmp%dleaf=0.0;lucmp%frac4=0.0;lucmp%ejmax=0.0;lucmp%vcmax=0.0;lucmp%rp20=0.0;lucmp%rpcoef=0.0
     lucmp%wai=0.0;lucmp%shelrb=0.0;lucmp%vegcf=0.0;lucmp%extkn=0.0;lucmp%tminvj=0.0;lucmp%tmaxvj=0.0;lucmp%vbeta=0.0
     lucmp%canst1=0.0;lucmp%xalbnir=0.0
     lucmp%albsoil=0.0
     lucmp%clay=0.0;lucmp%sand=0.0;lucmp%silt=0.0;lucmp%ssat=0.0;lucmp%sfc=0.0;lucmp%swilt=0.0
     lucmp%bch=0.0;lucmp%hyds=0.0;lucmp%sucs=0.0;lucmp%css=0.0;lucmp%rhosoil=0.0

     lucmp%rtsoil=0.0;lucmp%wetfac=0.0;lucmp%tss=0.0;lucmp%runoff=0.0;lucmp%rnof1=0.0;lucmp%rnof2=0.0
     lucmp%cansto=0.0;lucmp%ghflux=0.0;lucmp%sghflux=0.0;lucmp%ga=0.0
     lucmp%dgdtg=0.0;lucmp%fev=0.0;lucmp%fes=0.0;lucmp%fhs=0.0;lucmp%wbtot0=0.0
     lucmp%osnowd0=0.0;lucmp%trad=0.0;lucmp%za_uv=0.0;lucmp%za_tq=0.0
     lucmp%ssdnn=0.0;lucmp%snowd=0.0;lucmp%snage=0.0;lucmp%osnowd=0.0
     lucmp%albedo=0.0;lucmp%albsoilsn=0.0
     lucmp%tggsn=0.0;lucmp%ssdn=0.0;lucmp%smass=0.0;lucmp%sdepth=0.0
     lucmp%tgg=0.0;lucmp%wb=0.0;lucmp%wbice=0.0;lucmp%gammzz=0.0;lucmp%froot=0.0
     lucmp%sumcbal=0.0;lucmp%sumnbal=0.0;lucmp%sumpbal=0.0

     lucmp%clabile = 0.0
     lucmp%cplant=0.0;  lucmp%nplant=0.0;  lucmp%pplant=0.0
     lucmp%clitter=0.0; lucmp%nlitter=0.0; lucmp%plitter=0.0
     lucmp%csoil=0.0;   lucmp%nsoil=0.0;   lucmp%psoil=0.0
     lucmp%nsoilmin=0.0
     lucmp%psoillab=0.0;lucmp%psoilsorb=0.0;lucmp%psoilocc=0.0
     lucmp%cwoodprod=0.0;lucmp%nwoodprod=0.0;lucmp%pwoodprod=0.0
   END SUBROUTINE landuse_allocate_mp


   SUBROUTINE landuse_deallocate_mp(mpx,lucmp)
   integer     mpx
   TYPE(landuse_mp), INTENT(INOUT)  :: lucmp

     deallocate(lucmp%iveg,lucmp%isoil,lucmp%soilorder,lucmp%phase,lucmp%meth,lucmp%meth)
     deallocate(lucmp%zse)
     deallocate(lucmp%patchfrac,lucmp%areacell,lucmp%lai,lucmp%slai)
     deallocate(lucmp%hc,lucmp%xfang,lucmp%dleaf,lucmp%frac4,lucmp%ejmax,lucmp%vcmax,lucmp%rp20,lucmp%rpcoef, &
                lucmp%wai,lucmp%shelrb,lucmp%vegcf,lucmp%extkn,lucmp%tminvj,lucmp%tmaxvj,lucmp%vbeta, &
                lucmp%canst1,lucmp%xalbnir)                              
     deallocate(lucmp%albsoil)
     deallocate(lucmp%clay,lucmp%sand,lucmp%silt,lucmp%ssat,lucmp%sfc,lucmp%swilt, &
                lucmp%bch,lucmp%hyds,lucmp%sucs,lucmp%css,lucmp%rhosoil)

     deallocate(lucmp%rtsoil,lucmp%wetfac,lucmp%tss,lucmp%runoff,lucmp%rnof1,lucmp%rnof2, &
                lucmp%cansto,lucmp%ghflux,lucmp%sghflux,lucmp%ga,                &
                lucmp%dgdtg,lucmp%fev,lucmp%fes,lucmp%fhs,lucmp%wbtot0,                       &
                lucmp%osnowd0,lucmp%trad,lucmp%za_uv,lucmp%za_tq,                  &
                lucmp%ssdnn,lucmp%snowd,lucmp%snage,lucmp%osnowd)
     deallocate(lucmp%albedo,lucmp%albsoilsn)
     deallocate(lucmp%tggsn,lucmp%ssdn,lucmp%smass,lucmp%sdepth)
     deallocate(lucmp%tgg,lucmp%wb,lucmp%wbice,lucmp%gammzz,lucmp%froot) 
     deallocate(lucmp%sumcbal,lucmp%sumnbal,lucmp%sumpbal)

     deallocate(lucmp%clabile)
     deallocate(lucmp%cplant,lucmp%nplant,lucmp%pplant)
     deallocate(lucmp%clitter,lucmp%nlitter,lucmp%plitter)
     deallocate(lucmp%csoil,lucmp%nsoil,lucmp%psoil)
     deallocate(lucmp%nsoilmin)
     deallocate(lucmp%psoillab,lucmp%psoilsorb,lucmp%psoilocc)
     deallocate(lucmp%cwoodprod,lucmp%nwoodprod,lucmp%pwoodprod)

   END SUBROUTINE landuse_deallocate_mp

END MODULE landuse_patch

  subroutine landuse_driver(ssnow,soil,veg,bal,canopy,rough,rad,phen,casapool,casabal,casamet)
! Re-allocate casacnp restart pools for new land patch array.
! Used for LAND USE CHANGE SIMULATION
! Call by casa_init
! Q.Zhang @ 04/05/2011
  USE cable_IO_vars_module, ONLY: mask,patch,landpt
  USE cable_def_types_mod,  ONLY: mp,mvtype,mstype,mland,r_2,ms,msn,nrb,     &
                                  soil_parameter_type, soil_snow_type,       &
                                  veg_parameter_type, roughness_type,        &
                                  balances_type, canopy_type, radiation_type
  USE casadimension,        ONLY: icycle,mplant,mlitter,msoil,mwood,mso
  USE casavariable,         ONLY: casa_pool,casa_balance,casa_met,casa_biome
  USE phenvariable,         ONLY: phen_variable
  USE landuse_variable
  USE landuse_patch
  IMPLICIT NONE
  TYPE (soil_snow_type)          :: ssnow   ! soil and snow variables
  TYPE (soil_parameter_type)     :: soil    ! soil parameters
  TYPE (veg_parameter_type)      :: veg     ! vegetation parameters
  TYPE (balances_type)           :: bal
  TYPE (canopy_type)             :: canopy
  TYPE (roughness_type)          :: rough   ! roughness varibles
  TYPE (radiation_type)          :: rad     ! radiation variables
  TYPE (phen_variable)           :: phen
  TYPE (casa_pool)               :: casapool
  TYPE (casa_balance)            :: casabal
  TYPE (casa_met)                :: casamet
  TYPE (casa_biome)              :: casabiome

  TYPE (landuse_type),save       :: luc
  TYPE (landuse_mp),save         :: lucmp

  integer     mlon,mlat
  integer,    dimension(:,:), allocatable            :: landmask
  ! "mland" variables
  integer,        dimension(:),        allocatable   :: cstart,cend,nap
  real(kind=dp),  dimension(:,:),      allocatable   :: fracpry
  real(kind=dp),  dimension(:),        allocatable   :: arealand
!  character*120   fxpft,fxluh2cable
  integer ivt, mm, np, np1, p, q, mpx

    ! the following variables are available from "CABLE"
    ! mlon,mlat,landmask  
     mlon     = SIZE(mask,DIM=1)  ! not landpt(:)%ilon
     mlat     = SIZE(mask,DIM=2)  ! not landpt(:)%ilat
     allocate(landmask(mlon,mlat))
     allocate(fracpry(mland,mvtype))
     allocate(cstart(mland),cend(mland),nap(mland))
     allocate(arealand(mland))
     landmask = mask
     DO mm = 1, mland
       cstart(mm)= landpt(mm)%cstart
       cend(mm)  = landpt(mm)%cend
       nap(mm)   = landpt(mm)%nap
       arealand(mm) = SUM(casamet%areacell(cstart(mm):cend(mm)))
     ENDDO
     ! arealand = ????? ! total area of land cell, not patch area (casamet%areacell)
     call landuse_allocate_mp(mp,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)
     call landuse_setup(mland,luc)            

     ! get the mapping matrix from state to PFT
     call landuse_getxluh2(mlat,mlon,landmask,luc)
     call landuse_getdata(mlat,mlon,landmask,luc)
     ! get pool sizes and other states in the "restart" file
     lucmp%zse(:)       = soil%zse(:)          
     lucmp%iveg(:)      = veg%iveg(:)
     lucmp%isoil(:)     = soil%isoilm(:)          
     lucmp%soilorder(:) = casamet%isorder(:)          
     lucmp%phase(:)     = phen%phase(:)          
     lucmp%isflag(:)    = ssnow%isflag(:)          
     lucmp%meth(:)      = veg%meth(:)           
     DO mm = 1, mland
       lucmp%patchfrac(mm) = patch(mm)%frac       ! maybe we should create another variable for "primary%patch"
     ENDDO
     lucmp%areacell(:)  = casamet%areacell(:)
     lucmp%lai(:)       = veg%vlai(:)
     lucmp%slai(:)      = casabiome%sla(veg%iveg(:)) 
           
     lucmp%rtsoil(:)    = ssnow%rtsoil(:)
     lucmp%wetfac(:)    = ssnow%wetfac(:)
     lucmp%tss(:)       = ssnow%tss(:)
     lucmp%runoff(:)    = ssnow%runoff(:)
     lucmp%rnof1(:)     = ssnow%rnof1(:)
     lucmp%rnof2(:)     = ssnow%rnof2(:)
     lucmp%cansto(:)    = canopy%cansto(:)
     lucmp%ghflux(:)    = canopy%ghflux(:)
     lucmp%sghflux(:)   = canopy%sghflux(:)
     lucmp%ga(:)        = canopy%ga(:)
     lucmp%dgdtg(:)     = canopy%dgdtg(:)
     lucmp%fev(:)       = canopy%fev(:)
     lucmp%fes(:)       = canopy%fes(:)
     lucmp%fhs(:)       = canopy%fhs(:)
     lucmp%wbtot0(:)    = bal%wbtot0(:)
     lucmp%osnowd0(:)   = bal%osnowd0(:)
     lucmp%trad(:)      = rad%trad(:)
     lucmp%za_uv(:)     = rough%za_uv(:)
     lucmp%za_tq(:)     = rough%za_tq(:)
     lucmp%ssdnn(:)     = ssnow%ssdnn(:)
     lucmp%snowd(:)     = ssnow%snowd(:)
     lucmp%snage(:)     = ssnow%snage(:)
     lucmp%osnowd(:)    = ssnow%osnowd(:)

     lucmp%albedo(:,:)  = rad%albedo(:,:)
     lucmp%albsoilsn(:,:) = ssnow%albsoilsn(:,:)

     lucmp%tggsn(:,:)    = ssnow%tggsn(:,:)
     lucmp%ssdn(:,:)     = ssnow%ssdn(:,:)
     lucmp%smass(:,:)    = ssnow%smass(:,:)
     lucmp%sdepth(:,:)   = ssnow%sdepth(:,:)
     lucmp%tgg(:,:)      = ssnow%tgg(:,:)
     lucmp%wb(:,:)       = ssnow%wb(:,:)
     lucmp%wbice(:,:)    = ssnow%wbice(:,:)
     lucmp%gammzz(:,:)   = ssnow%gammzz(:,:)
     lucmp%froot(:,:)    = veg%froot(:,:)
    
     do mm=1,mland
       do np=cstart(mm),cend(mm)
          ivt = lucmp%iveg(np)
          if(ivt <=mvtype) then
             luc%primary(mm,ivt) = patch(np)%frac
          endif
       enddo
     enddo

     if(icycle>0) then 
        lucmp%phase(:)       = phen%phase(:)
        lucmp%clabile(:)     = casapool%clabile(:)
        lucmp%cplant(:,:)    = casapool%cplant(:,:)
        lucmp%clitter(:,:)   = casapool%clitter(:,:)
        lucmp%csoil(:,:)     = casapool%csoil(:,:)
        lucmp%cwoodprod(:,:) = casapool%cwoodprod(:,:)
        lucmp%sumcbal(:)     = casabal%sumcbal(:)
     endif
     if(icycle>1) then
        lucmp%nplant(:,:)    = casapool%nplant(:,:)
        lucmp%nlitter(:,:)   = casapool%nlitter(:,:)
        lucmp%nsoil(:,:)     = casapool%nsoil(:,:)
        lucmp%nwoodprod(:,:) = casapool%nwoodprod(:,:)
        lucmp%nsoilmin(:)    = casapool%nsoilmin(:)
        lucmp%sumnbal(:)     = casabal%sumnbal(:)
     endif
     if(icycle >2) then
        lucmp%pplant(:,:)    = casapool%pplant(:,:)
        lucmp%plitter(:,:)   = casapool%plitter(:,:)
        lucmp%psoil(:,:)     = casapool%psoil(:,:)
        lucmp%pwoodprod(:,:) = casapool%pwoodprod(:,:)
        lucmp%psoillab(:)    = casapool%psoillab(:)
        lucmp%psoilsorb(:)   = casapool%psoilsorb(:)
        lucmp%psoilocc(:)    = casapool%psoilocc(:)
        lucmp%sumpbal(:)     = casabal%sumpbal(:)
     endif
      
     ! assign variables var(mp,:) to luc%var_x(mland,mvmax,:)
     call landuse_mp2land(luc,mp,cstart,cend,lucmp)
     ! we need to deallocate "lucmp" because "mp" will be updated after land use change
!     call landuse_deallocate_mp(mp,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)
     call landuse_deallocate_mp(mp,lucmp)
     call landuse_transitx(icycle,luc)
     call landuse_checks(mlon,mlat,landmask,luc,arealand)
     call landuse_update(luc)                    ! assign "var_y" to "var_x"
     ! update "cstart", "cend","nap" and "mp=>mpx"
      cstart=0;cend=0;nap=0
      np =0; cstart(:) = 0; cend(:) =0; nap(:) = 0
      do p=1,mland
         np1 =0
         if(sum(luc%patchfrac_y(p,:))<fracmin) then
            print *, 'WARNING! patch area sum too low',p,luc%patchfrac_y(p,:)
         else
            do q=1,mvmax
               if(luc%patchfrac_y(p,q) >fracmin) then
                  np  = np + 1
                  np1 = np1 + 1
                  if(np1==1) cstart(p) = np
              endif
            enddo
            cend(p) = np
            nap(p)  = max(0,cend(p)-cstart(p) +1)
         endif
      enddo
      mpx = np
      call landuse_allocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)

     ! allocate "lucmp" with "mpx"
     call landuse_allocate_mp(mpx,ms,msn,nrb,mplant,mlitter,msoil,mwood,lucmp)
     call landuse_land2mpx(luc,lucmp,mpx,cstart,cend,nap)
     call landuse_deallocate(luc)

     deallocate(fracpry)
     deallocate(cstart,cend,nap)
     deallocate(arealand)

     close(21)
211  format(i4,a120)
end subroutine landuse_driver

 SUBROUTINE landuse_mp2land(luc,mp,cstart,cend,lucmp)
 use landuse_variable
 USE landuse_patch
 IMPLICIT NONE
 type(landuse_type)    :: luc
 type(landuse_mp)      :: lucmp
 integer mp
 integer g,np,ivt
 integer,        dimension(mland)        :: cstart,cend

   !initilize all old state variables
   !biogeochemical states
   luc%phase_x(:,:) = -1;     luc%patchfrac_x = 0.0
   luc%cplant_x(:,:,:) = 0.0; luc%clitter_x(:,:,:) =0.0; luc%csoil_x(:,:,:) = 0.0
   luc%clabile_x(:,:)  = 0.0; luc%cwoodprod_x(:,:,:)=0.0
   luc%nplant_x(:,:,:) = 0.0; luc%nlitter_x(:,:,:) =0.0; luc%nsoil_x(:,:,:) = 0.0
   luc%nsoilmin_x(:,:) = 0.0; luc%nwoodprod_x(:,:,:)=0.0
   luc%pplant_x(:,:,:) = 0.0; luc%plitter_x(:,:,:) =0.0; luc%psoil_x(:,:,:) = 0.0
   luc%psoillab_x(:,:) = 0.0; luc%pwoodprod_x(:,:,:)=0.0
   luc%psoilsorb_x(:,:)= 0.0; luc%psoilocc_x(:,:)  =0.0

 ! read the CABLE restart file to get the state variables (both biophysical and biogeochemical), and mland, nap(mland)

 ! locate variables to mvtype tiles per grid. 


  do g=1,mland
  do np= cstart(g),cend(g)
     ivt = lucmp%iveg(np) 
     luc%patchfrac_x(g,ivt)   = lucmp%patchfrac(np)
         
     luc%phase_x(g,ivt)       = lucmp%phase(np)
     luc%clabile_x(g,ivt)     = lucmp%clabile(np)
     luc%cplant_x(g,ivt,:)    = lucmp%cplant(np,:)
     luc%clitter_x(g,ivt,:)   = lucmp%clitter(np,:)
     luc%csoil_x(g,ivt,:)     = lucmp%csoil(np,:)
     luc%cwoodprod_x(g,ivt,:) = lucmp%cwoodprod(np,:)

    IF(icycle>1) THEN
      luc%nplant_x(g,ivt,:)    = lucmp%nplant(np,:)
      luc%nlitter_x(g,ivt,:)   = lucmp%nlitter(np,:)
      luc%nsoil_x(g,ivt,:)     = lucmp%nsoil(np,:)
      luc%nsoilmin_x(g,ivt)    = lucmp%nsoilmin(np)
      luc%nwoodprod_x(g,ivt,:) = lucmp%nwoodprod(np,:)
    END IF
    IF(icycle>2) THEN
      luc%pplant_x(g,ivt,:)    = lucmp%pplant(np,:)
      luc%plitter_x(g,ivt,:)   = lucmp%plitter(np,:)
      luc%psoil_x(g,ivt,:)     = lucmp%psoil(np,:)
      luc%psoillab_x(g,ivt)    = lucmp%psoillab(np)
      luc%psoilsorb_x(g,ivt)   = lucmp%psoilsorb(np)
      luc%psoilocc_x(g,ivt)    = lucmp%psoilocc(np)
      luc%pwoodprod_x(g,ivt,:) = lucmp%pwoodprod(np,:)
    END IF
    
  enddo
  enddo

  luc%patchfrac_y = luc%patchfrac_x
  luc%phase_y     = luc%phase_x
  luc%clabile_y   = luc%clabile_x
  luc%cplant_y    = luc%cplant_x
  luc%clitter_y   = luc%clitter_x
  luc%csoil_y     = luc%csoil_x
  luc%cwoodprod_y = luc%cwoodprod_x

  IF(icycle>1) THEN
     luc%nplant_y    = luc%nplant_x
     luc%nlitter_y   = luc%nlitter_x
     luc%nsoil_y     = luc%nsoil_x
     luc%nsoilmin_y  = luc%nsoilmin_x
     luc%nwoodprod_y = luc%nwoodprod_x
  END IF
  IF(icycle>2) THEN
     luc%pplant_y    = luc%pplant_x
     luc%plitter_y   = luc%plitter_x
     luc%psoil_y     = luc%psoil_x
     luc%psoillab_y  = luc%psoillab_x
     luc%psoilsorb_y = luc%psoilsorb_x
     luc%psoilocc_y  = luc%psoilocc_x
     luc%pwoodprod_y = luc%pwoodprod_x
  END IF

END SUBROUTINE landuse_mp2land
  

SUBROUTINE landuse_dims(fpft0,mlon,mlat,mp,landmask,lon,lat)
!  calculate "mland", "landmask",
!  and initial patchfrac
use netcdf
use landuse_constant
implicit none
character*120 fpft0
integer    mlon,mlat,mp
integer,   dimension(mlon,mlat)              :: landmask
real(kind=dp),  dimension(mlon)              :: lon
real(kind=dp),  dimension(mlat)              :: lat
! local variables
real(kind=dp),   dimension(:,:,:),   allocatable   :: primary
real(kind=dp),   dimension(:,:,:),   allocatable   :: secondary
integer ok,ncid1,varxid,i,j,k,np,nland

    allocate(primary(mlon,mlat,mvtype))
    allocate(secondary(mlon,mlat,mvtype))

    ok = nf90_open(fpft0,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"longitude",varxid)
    ok = nf90_get_var(ncid1,varxid,lon)
    ok = nf90_inq_varid(ncid1,"latitude",varxid)
    ok = nf90_get_var(ncid1,varxid,lat)
    ok = nf90_inq_varid(ncid1,"primary",varxid)
    ok = nf90_get_var(ncid1,varxid,primary)
    ok = nf90_inq_varid(ncid1,"secondary",varxid)
    ok = nf90_get_var(ncid1,varxid,secondary)
    ok = nf90_close(ncid1)

    landmask(:,:) = 0
    nland = 0
    np    = 0
    do i=1,mlon
    do j=1,mlat
       
       if(sum(primary(i,j,1:mvtype)+secondary(i,j,1:mvtype))>fracmin) then
          landmask(i,j)=1
          nland = nland +1
          do k=1,4
             if(primary(i,j,k)>fracmin) then
                np = np +1
             endif
          enddo
          do k=5,mvtype
             if((primary(i,j,k)+secondary(i,j,k))>fracmin) then
                np = np +1
             endif
          enddo
          do k=1,4
             if(secondary(i,j,k)>fracmin) then
                np = np +1
             endif
          enddo
       endif
    enddo
    enddo

    mland = nland
    mp = np

    print *, 'total land cells and pathes are ', mland,mp
    deallocate(primary)
    deallocate(secondary)

END SUBROUTINE landuse_dims

SUBROUTINE landuse_init(fpft0,mlon,mlat,mp,landmask,fracpry,cstart,cend,lucmp) 
  !  calculate "mland", "landmask", 
  !  and initial patchfrac
  use netcdf
  use landuse_patch
  implicit none
  TYPE(landuse_mp) :: lucmp
  character*120 fpft0
  integer mlon,mlat,mp
  ! output 
  integer,   dimension(mlon,mlat)                    :: landmask
  real(kind=dp), dimension(mland,mvtype)             :: fracpry  
  integer,   dimension(mland)                        :: cstart    
  integer,   dimension(mland)                        :: cend     
  ! local variables
  real(kind=dp),   dimension(:,:,:),   allocatable   :: primary  
  real(kind=dp),   dimension(:,:,:),   allocatable   :: secondary
  integer ok,ncid1,varxid,i,j,k,np,np1,nland

    allocate(primary(mlon,mlat,mvtype))
    allocate(secondary(mlon,mlat,mvtype))

    ok = nf90_open(fpft0,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"primary",varxid)
    ok = nf90_get_var(ncid1,varxid,primary)
    ok = nf90_inq_varid(ncid1,"secondary",varxid)
    ok = nf90_get_var(ncid1,varxid,secondary)
    ok = nf90_close(ncid1)

    nland = 0
    np    = 0
    fracpry = 0.0; lucmp%patchfrac = 0.0; lucmp%iveg = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
         nland = nland + 1
         do k=1,mvtype
            fracpry(nland,k)   = primary(i,j,k)
         enddo
         
         np1 = 0
         do k=1,4
            if(primary(i,j,k) > fracmin) then
               np = np +1
               np1 = np1 + 1
               if(np1==1) cstart(nland) = np
               lucmp%patchfrac(np) = primary(i,j,k) 
               lucmp%iveg(np) = k
            endif
         enddo
         do k=5,mvtype
            if((primary(i,j,k)+secondary(i,j,k))>fracmin) then
               np = np +1
               np1 = np1 +1 
               if(np1==1) cstart(nland) = np
               lucmp%patchfrac(np) = primary(i,j,k) + secondary(i,j,k)
               lucmp%iveg(np) = k
            endif
         enddo
         do k=1,4
            if(secondary(i,j,k)>fracmin) then
               np = np +1
               np1 = np1 +1 
               if(np1==1) cstart(nland) = np
               lucmp%patchfrac(np) = secondary(i,j,k)
               lucmp%iveg(np) = k+mvtype
            endif
         enddo
         cend(nland) = np
        endif
     enddo
     enddo
 
    deallocate(primary)
    deallocate(secondary)

end SUBROUTINE landuse_init

SUBROUTINE landuse_getxluh2(mlat,mlon,landmask,luc)
! get data: luc%fprimary; luc%fsecondary
  USE netcdf
  use landuse_variable
  IMPLICIT NONE
  TYPE(landuse_type)   :: luc
!  character*120 fxluh2cable
  integer mp,mlat,mlon
  integer,   dimension(mlon,mlat)                  :: landmask
  ! local variables
  real(kind=dp), dimension(:,:,:,:), allocatable   :: xluh2cable
  integer ok,ncid2,varxid
  integer i,j,m,v,s

    allocate(xluh2cable(mlon,mlat,mvmax,mstate))
    ok = nf90_open(fxluh2cable,nf90_nowrite,ncid2)
    ok = nf90_inq_varid(ncid2,"xluh2cable",varxid)
    ok = nf90_get_var(ncid2,varxid,xluh2cable)
    ok = nf90_close(ncid2)
    ! assig the values of luc%variables
    luc%xluh2cable(:,:,:) = 0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          luc%xluh2cable(m,:,:) = xluh2cable(i,j,:,:)
          do s=1,mstate
             do v=1,mvmax
                luc%xluh2cable(m,v,s) = luc%xluh2cable(m,v,s)/sum(luc%xluh2cable(m,1:mvmax,s))
             enddo
          enddo
       endif
    enddo
    enddo

    deallocate(xluh2cable)

 END SUBROUTINE landuse_getxluh2


SUBROUTINE landuse_getdata(mlat,mlon,landmask,luc)
! get data: luc%fprimary; luc%fsecondary
  USE netcdf
  use landuse_variable
  IMPLICIT NONE
  TYPE(landuse_type)   :: luc
!  character*120 fxpft01
  integer mlat,mlon
  integer,   dimension(mlon,mlat)                  :: landmask 
  ! local variables
  real(kind=dp), dimension(:,:,:),   allocatable   :: fracharvw
  real(kind=dp), dimension(:,:,:,:), allocatable   :: transitx
  integer ok,ncid1,varxid
  integer i,j,m,k

    allocate(fracharvw(mlon,mlat,mharvw))
    allocate(transitx(mlon,mlat,mvmax,mvmax))

    ok = nf90_open(fxpft,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"harvest",varxid)
    ok = nf90_get_var(ncid1,varxid,fracharvw)
    ok = nf90_inq_varid(ncid1,"transition",varxid)
    ok = nf90_get_var(ncid1,varxid,transitx)
    ok = nf90_close(ncid1)

    ! assig the values of luc%variables
    luc%fharvw(:,:) =0.0; luc%atransit(:,:,:)=0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          luc%atransit(m,:,:)   = transitx(i,j,:,:)
          luc%fharvw(m,:)       = fracharvw(i,j,:)
       endif
    enddo
    enddo
    
    deallocate(fracharvw)
    deallocate(transitx)

END SUBROUTINE landuse_getdata

SUBROUTINE landuse_transitx(luc)
  
   USE casaparm
   USE casavariable,  ONLY: casa_biome
   USE landuse_variable
   IMPLICIT NONE
   TYPE(casa_biome)    :: casabiome
   TYPE(landuse_type)  :: luc
   !INTEGER             :: mvtype1  ! not used
   !integer, parameter        :: mvtype1 = mvtype+1 ! wrong, mvtype has no number yet
!!  integer, parameter        :: leaf =1
!!  integer, parameter        :: wood =2
!!  integer, parameter        :: froot =3
!!  integer, parameter        :: metb =1
!!  integer, parameter        :: str =2
!!  integer, parameter        :: cwd =3
!!  integer, parameter        :: mic =1
!!  integer, parameter        :: slow=2
!!   ! local variables
!!   real, dimension(mvtype,3)  :: ratioNCplantmax=(/0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02, &
!!                                                   0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005,0.005, &
!!                                                   0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01/)
!!
!!   real, dimension(mvtype,3)  :: ratioNPplantmin=(/15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0, &
!!                                                  15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0, &
!!                                                  15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0,15.0/)
!!
!!
!!   real, dimension(mvtype,3)  :: ftransNPtoL =(/0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,                  &
!!                                                0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95,0.95, &
!!                                                0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9/)
!!
!!   real, dimension(mvtype,3)  :: fracLigninplant =(/0.25,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15,0.25,0.1, &
!!                                                    0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,        &
!!                                                    0.25,0.2,0.2,0.2,0.2,0.1,0.1,0.1,0.1,0.1,0.15,0.15,0.15,0.15,0.15,0.25,0.1/)
!!
!   integer,       dimension(mvtype)               :: ivt2=(/3,3,3,3,2,1,1,2,1,1,0,0,0,1,0,0,0/)
   INTEGER,       DIMENSION(mvtype)               :: ivt2
!   DATA  ivt2 /3,3,3,3,2,1,1,2,1,1,0,0,0,1,0,0,0/
   real(kind=dp), dimension(mland,mvmax)          :: dclabile
   real(kind=dp), dimension(mland,mvmax,mplant)   :: dcplant,dnplant,dpplant
   real(kind=dp), dimension(mland,mvmax,mlitter)  :: dclitter,dnlitter,dplitter
   real(kind=dp), dimension(mland,mvmax,msoil)    :: dcsoil,dnsoil,dpsoil
   real(kind=dp), dimension(mland,mvmax)          :: dnsoilmin,dpsoillab,dpsoilsorb,dpsoilocc
   real(kind=dp), dimension(mland,mvmax,mwood)    :: dcwoodprod,dnwoodprod,dpwoodprod
   real(kind=dp), DIMENSION(3)                    :: ratioLignintoN
   real(kind=dp), DIMENSION(mvmax,3,3)            :: fromPtoL
   real(kind=dp), DIMENSION(mland,mvmax)          :: delarea
   real(kind=dp), DIMENSION(mvmax)                :: delareax
   real(kind=dp), DIMENSION(3)                    :: totcwoodprod, totclitter, totcsoil
   real(kind=dp), DIMENSION(3)                    :: totnwoodprod, totnlitter, totnsoil
   real(kind=dp), DIMENSION(3)                    :: totpwoodprod, totplitter, totpsoil
   real(kind=dp)                                     totclabile, totnsoilmin, totpsoillab, totpsoilsorb, totpsoilocc
   real(kind=dp)                                     tempx
   integer p,d,r,q,r1,r2,r3,r4,ierror,ivt,k

   ivt2 = (/3,3,3,3,2,1,1,2,1,1,0,0,0,1,0,0,0/)

   delarea          = 0.0
   dcplant(:,:,:)   = 0.0; dnplant(:,:,:)   = 0.0; dpplant(:,:,:)    = 0.0; dclabile(:,:) = 0.0
   dclitter(:,:,:)  = 0.0; dnlitter(:,:,:)  = 0.0; dplitter(:,:,:)   = 0.0
   dcsoil(:,:,:)    = 0.0; dnsoil(:,:,:)    = 0.0; dpsoil(:,:,:)     = 0.0
   dnsoilmin(:,:)   = 0.0; dpsoillab(:,:)   = 0.0; dpsoilsorb(:,:)   = 0.0; dpsoilocc(:,:) = 0.0
   dcwoodprod(:,:,:) =0.0; dnwoodprod(:,:,:)=0.0;  dpwoodprod(:,:,:) = 0.0

   do p = 1,mland
      fromPtoL(:,1,:)=1.0; fromPtoL(:,2,:) = 0.0
      do d=1,mvmax
         if(luc%cplant_x(p,d,leaf) > 0.001) then
            ! calculate the fraction of litter or root litter into metabolic litter pool
            if(d>mvtype) then
               ivt = d-mvtype
            else
               ivt=mvtype
            endif

            ratioLignintoN(leaf) = (luc%cplant_x(p,d,leaf) &
                                 /(max(1.0e-10,luc%nplant_x(p,d,leaf)) * casabiome%ftransNPtoL(ivt,leaf))) &
                                 * casabiome%fracLigninplant(ivt,leaf)
            ratioLignintoN(froot)= (luc%cplant_x(p,d,froot)&
                                 /(max(1.0e-10,luc%nplant_x(p,d,froot))* casabiome%ftransNPtoL(ivt,froot))) &
                                 * casabiome%fracLigninplant(ivt,froot)
            fromPtoL(d,metb,leaf)  = max(0.001, 0.85 - 0.018 *ratioLignintoN(leaf))
            fromPtoL(d,metb,froot) = max(0.001, 0.85 - 0.018 *ratioLignintoN(froot))
            fromPtoL(d,str,leaf)   = 1.0 - fromPtoL(d,metb,leaf)
            fromPtoL(d,str,froot)  = 1.0 - fromPtoL(d,metb,froot)
         endif
      enddo

      do d = 1,mvmax
      do r = 1,mvmax
         ! transfer leaf and root into litter (metabolic and structural litter),
         ! transfer wood into wood prodoct pool (three wood product pools)
         ! then calculate the delpool(p,d,:), and delpool(p,r,:) for C, N and P

         if(luc%atransit(p,r,d) > 0.0.and.d/=r.and.luc%patchfrac_x(p,d)>0.0) then
            ! transfer the area from donor (d) to receiver (r)
            delarea(p,d) = delarea(p,d) - luc%atransit(p,r,d)
            delarea(p,r) = delarea(p,r) + luc%atransit(p,r,d)
            ! donor pool changes
            dcplant(p,d,:)    = dcplant(p,d,:)    - luc%atransit(p,r,d) * luc%cplant_x(p,d,:)
            dclitter(p,d,:)   = dclitter(p,d,:)   - luc%atransit(p,r,d) * luc%clitter_x(p,d,:)
            dcsoil(p,d,:)     = dcsoil(p,d,:)     - luc%atransit(p,r,d) * luc%csoil_x(p,d,:)
            dclabile(p,d)     = dclabile(p,d)     - luc%atransit(p,r,d) * luc%clabile_x(p,d)
            dcwoodprod(p,d,:) = dcwoodprod(p,d,:) - luc%atransit(p,r,d) * luc%cwoodprod_x(p,d,:)
            ! receiver pool changes
            ! move the donor leaf and root biomass to receiver structural litter pool
            !(to avoid C:N imbalance)
            ! using max function to avoid dividing by zero, ypw 14/may/2008

            ! calculate the fraction of litter or root litter into metabolic litter pool

            dclitter(p,r,1) = dclitter(p,r,1) + luc%atransit(p,r,d)  &
                            * (luc%clitter_x(p,d,1) + luc%cplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                    + luc%cplant_x(p,d,froot)* fromPtoL(d,metb,froot))
            dclitter(p,r,2) = dclitter(p,r,2) + luc%atransit(p,r,d)  &
                            * (luc%clitter_x(p,d,2) + luc%cplant_x(p,d,leaf) * fromPtoL(d,str,leaf)  &
                                                    + luc%cplant_x(p,d,froot) *fromPtoL(d,str,froot))
            dcsoil(p,r,:)   = dcsoil(p,r,:)   + luc%atransit(p,r,d) * luc%csoil_x(p,d,:)
            !move the labile carbon to receiving tile (maybe better to fast-decomposing wood product pool)
            dclabile(p,r)   = dclabile(p,r)   + luc%atransit(p,r,d) * luc%clabile_x(p,d)
            !move the donor wood to receiver wood product pool including labile C
            dcwoodprod(p,r,1) = dcwoodprod(p,r,1) + luc%atransit(p,r,d)*(luc%cplant_x(p,d,wood)*fwoodprod(1) + luc%cwoodprod_x(p,d,1)) &
                                                  + dclabile(p,d)
            dcwoodprod(p,r,2) = dcwoodprod(p,r,2) + luc%atransit(p,r,d)*(luc%cplant_x(p,d,wood)*fwoodprod(2) + luc%cwoodprod_x(p,d,2))
            dcwoodprod(p,r,3) = dcwoodprod(p,r,3) + luc%atransit(p,r,d)*(luc%cplant_x(p,d,wood)*fwoodprod(3) + luc%cwoodprod_x(p,d,3))

             if(icycle >1) then
                dnplant(p,d,:)  = dnplant(p,d,:)  - luc%atransit(p,r,d) * luc%nplant_x(p,d,:)
                dnlitter(p,d,:) = dnlitter(p,d,:) - luc%atransit(p,r,d) * luc%nlitter_x(p,d,:)
                dnsoil(p,d,:)   = dnsoil(p,d,:)   - luc%atransit(p,r,d) * luc%nsoil_x(p,d,:)
                dnsoilmin(p,d)  = dnsoilmin(p,d)  - luc%atransit(p,r,d) * luc%nsoilmin_x(p,d)

                dnlitter(p,r,1) = dnlitter(p,r,1) + luc%atransit(p,r,d)  &
                                * (luc%nlitter_x(p,d,1) + luc%nplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                        + luc%nplant_x(p,d,froot)* fromPtoL(d,metb,froot))
                dnlitter(p,r,2) = dnlitter(p,r,2) + luc%atransit(p,r,d)  &
                                * (luc%nlitter_x(p,d,2) + luc%nplant_x(p,d,leaf) * fromPtoL(d,str,leaf) &
                                                        + luc%nplant_x(p,d,froot) *fromPtoL(d,str,froot))
                dnsoil(p,r,:)   = dnsoil(p,r,:)   + luc%atransit(p,r,d) * luc%nsoil_x(p,d,:)
                dnsoilmin(p,r)  = dnsoilmin(p,r)  + luc%atransit(p,r,d) * luc%nsoilmin_x(p,d)
                dnwoodprod(p,r,:) = dnwoodprod(p,r,:) + luc%atransit(p,r,d) *(fwoodprod(:)*luc%nplant_x(p,d,wood) + luc%nwoodprod_x(p,d,:))
             endif
             if(icycle >2) then
                dpplant(p,d,:)  = dpplant(p,d,:)  - luc%atransit(p,r,d) * luc%pplant_x(p,d,:)
                dplitter(p,d,:) = dplitter(p,d,:) - luc%atransit(p,r,d) * luc%nlitter_x(p,d,:)
                dpsoil(p,d,:)   = dpsoil(p,d,:)   - luc%atransit(p,r,d) * luc%psoil_x(p,d,:)
                dpsoillab(p,d)  = dpsoillab(p,d)  - luc%atransit(p,r,d) * luc%psoillab_x(p,d)
                dpsoilsorb(p,d) = dpsoilsorb(p,d) - luc%atransit(p,r,d) * luc%psoilsorb_x(p,d)
                dpsoilocc(p,d)  = dpsoilocc(p,d)  - luc%atransit(p,r,d) * luc%psoilocc_x(p,d)

                dplitter(p,r,1) = dplitter(p,r,1) + luc%atransit(p,r,d)  &
                                * (luc%plitter_x(p,d,1) + luc%pplant_x(p,d,leaf) * fromPtoL(d,metb,leaf) &
                                                        + luc%pplant_x(p,d,froot)* fromPtoL(d,metb,froot))
                dplitter(p,r,2) = dplitter(p,r,2) + luc%atransit(p,r,d)  &
                                * (luc%plitter_x(p,d,2) + luc%pplant_x(p,d,leaf) * fromPtoL(d,str,leaf) &
                                                        + luc%pplant_x(p,d,froot) *fromPtoL(d,str,froot))
                dpsoil(p,r,:)   = dpsoil(p,r,:)   + luc%atransit(p,r,d) * luc%psoil_x(p,d,:)
                dpsoillab(p,r)  = dpsoillab(p,r)  + luc%atransit(p,r,d) * luc%psoillab_x(p,d)
                dpsoilsorb(p,r) = dpsoilsorb(p,r) + luc%atransit(p,r,d) * luc%psoilsorb_x(p,d)
                dpsoilocc(p,r)  = dpsoilocc(p,r)  + luc%atransit(p,r,d) * luc%psoilocc_x(p,d)
                dpwoodprod(p,r,:) = dpwoodprod(p,r,:) + luc%atransit(p,r,d) * (fwoodprod(:) * luc%pplant_x(p,d,wood) +luc%pwoodprod_x(p,d,:))
             endif
         endif  ! of "atransit" >0.0, "luctype" >0 etc.
      enddo   ! of "d_tile"
      enddo   ! of "r_tile"
      ! here we deal with wood harvest from primary forest, as a result primart forest into secondary forest
      delareax(:) =0.0;     totclabile = 0.0;  totcwoodprod(:) =0.0

      if(luc%fharvw(p,1) >0.0) then
         do d=1,mvtype
            delareax(d)  = luc%fharvw(p,1) * luc%xluh2cable(p,d,1)
            delarea(p,d) = delarea(p,d) -  delareax(d)
         enddo  ! of "d"

         ! all plant, litter and soil pools are amalgamated
         ! donor pool changes
         dcplant(p,d,:)    = dcplant(p,d,:)    - delareax(d) * luc%cplant_x(p,d,:)
         dclabile(p,d)     = dclabile(p,d)     - delareax(d) * luc%clabile_x(p,d)
         dcwoodprod(p,d,:) = dcwoodprod(p,d,:) - delareax(d) * luc%cwoodprod_x(p,d,:)
         dclitter(p,d,:)   = dclitter(p,d,:)   - delareax(d) * luc%clitter_x(p,d,:)
         dcsoil(p,d,:)     = dcsoil(p,d,:)     - delareax(d) * luc%csoil_x(p,d,:)

         if(icycle>1) then
            dnplant(p,d,:)    = dnplant(p,d,:)    - delareax(d) * luc%nplant_x(p,d,:)
            dnwoodprod(p,d,:) = dnwoodprod(p,d,:) - delareax(d) * luc%nwoodprod_x(p,d,:)
            dnlitter(p,d,:)   = dnlitter(p,d,:)   - delareax(d) * luc%nlitter_x(p,d,:)
            dnsoil(p,d,:)     = dnsoil(p,d,:)     - delareax(d) * luc%nsoil_x(p,d,:)
            dnsoilmin(p,d)    = dnsoilmin(p,d)    - delareax(d) * luc%nsoilmin_x(p,d)
         endif
         if(icycle>2) then
            dpplant(p,d,:)    = dpplant(p,d,:)    - delareax(d) * luc%pplant_x(p,d,:)
            dpwoodprod(p,d,:) = dpwoodprod(p,d,:) - delareax(d) * luc%pwoodprod_x(p,d,:)
            dplitter(p,d,:)   = dplitter(p,d,:)   - delareax(d) * luc%nlitter_x(p,d,:)
            dpsoil(p,d,:)     = dpsoil(p,d,:)     - delareax(d) * luc%psoil_x(p,d,:)
            dpsoillab(p,d)    = dpsoillab(p,d)    - delareax(d) * luc%psoillab_x(p,d)
            dpsoilsorb(p,d)   = dpsoilsorb(p,d)   - delareax(d) * luc%psoilsorb_x(p,d)
            dpsoilocc(p,d)    = dpsoilocc(p,d)    - delareax(d) * luc%psoilocc_x(p,d)
         endif

         do k=1,mwood
            totcwoodprod(k)    = totcwoodprod(k) + delareax(d) * (fwoodprod(k)* luc%cplant_x(p,d,2) + luc%cwoodprod_x(p,d,k))
         enddo

         !  add labile C to fast-decomposing wood product pool
         totcwoodprod(1)    = totcwoodprod(1) + delareax(d) * luc%clabile_x(p,d)

         totclitter(1)      = totclitter(1)   + delareax(d) * (luc%clitter_x(p,d,1)  &
                            + luc%cplant_x(p,d,leaf)*fromPtoL(d,metb,leaf)+luc%cplant_x(p,d,froot)*fromPtoL(d,metb,froot))
         totclitter(2)      = totclitter(2)   + delareax(d) * (luc%clitter_x(p,d,2)  &
                            + luc%cplant_x(p,d,leaf)*fromPtoL(d,str,leaf)+luc%cplant_x(p,d,froot)*fromPtoL(d,str,froot))
         do k=1,msoil
            totcsoil(k)     = totcsoil(k)     + delareax(d) * luc%csoil_x(p,d,k)
         enddo

         if(icycle>1) then
            totnwoodprod(:)    = totnwoodprod(:) + delareax(d) * (fwoodprod(:)* luc%nplant_x(p,d,2)+ luc%nwoodprod_x(p,d,:))
            totnlitter(1)      = totnlitter(1)   + delareax(d) * (luc%nlitter_x(p,d,1)  &
                               + luc%nplant_x(p,d,leaf)*fromPtoL(d,metb,leaf)+luc%nplant_x(p,d,froot)*fromPtoL(d,metb,froot))
            totnlitter(2)      = totnlitter(2)   + delareax(d) * (luc%nlitter_x(p,d,2)  &
                               + luc%nplant_x(p,d,leaf)*fromPtoL(d,str,leaf)+luc%nplant_x(p,d,froot)*fromPtoL(d,str,froot))
            totnsoil(:)        = totnsoil(:)     + delareax(d) * luc%nsoil_x(p,d,:)
            totnsoilmin        = totnsoilmin     + delareax(d) * luc%nsoilmin_x(p,d)
         endif
         if(icycle>2) then
            totpwoodprod(:)    = totpwoodprod(:) + delareax(d) * (fwoodprod(:)* luc%pplant_x(p,d,2)+luc%pwoodprod_x(p,d,:))
            totplitter(1)      = totplitter(1)   + delareax(d) * (luc%plitter_x(p,d,1)  &
                               + luc%pplant_x(p,d,leaf)*fromPtoL(d,metb,leaf)+luc%pplant_x(p,d,froot)*fromPtoL(d,metb,froot))
            totplitter(2)      = totplitter(2)   + delareax(d) * (luc%plitter_x(p,d,2)  &
                               + luc%pplant_x(p,d,leaf)*fromPtoL(d,str,leaf)+luc%pplant_x(p,d,froot)*fromPtoL(d,str,froot))
            totpsoil(:)        = totpsoil(:)     + delareax(d) * luc%psoil_x(p,d,:)
            totpsoillab        = totpsoillab     + delareax(d) * luc%psoillab_x(p,d)
            totpsoilsorb       = totpsoilsorb    + delareax(d) * luc%psoilsorb_x(p,d)
            totpsoilocc        = totpsoilocc     + delareax(d) * luc%psoilocc_x(p,d)
         endif


         do r=mvtype+1,mvmax  !18,21
            !! tempx             = luc%xluh2cable(p,r,3)/sum(luc%xluh2cable(p,18:21,3)+1.0e-10) 
            tempx             = luc%xluh2cable(p,r,3)/sum(luc%xluh2cable(p,(mvtype+1):mvmax,3)+1.0e-10) 
            delarea(p,r)      = delarea(p,r)      + luc%fharvw(p,1) * luc%xluh2cable(p,r,3)

            dcwoodprod(p,r,:) = dcwoodprod(p,r,:) + tempx * totcwoodprod(:)

    !     dclabile(p,r)     = dclabile(p,r)     + tempx * totclabile

            dclitter(p,r,:)   = dclitter(p,r,:)   + tempx * totclitter(:)
            dcsoil(p,r,:)     = dcsoil(p,r,:)     + tempx * totcsoil(:)
            if(icycle >1) then
               dnwoodprod(p,r,:) = dnwoodprod(p,r,:) + tempx * totnwoodprod(:)
               dnlitter(p,r,:)   = dnlitter(p,r,:)   + tempx * totnlitter(:)
               dnsoil(p,r,:)     = dnsoil(p,r,:)     + tempx * totnsoil(:)
               dnsoilmin(p,r)    = dnsoilmin(p,r)    + tempx * totnsoilmin
            endif
            if(icycle >2) then
               dpwoodprod(p,r,:) = dpwoodprod(p,r,:) + tempx * totpwoodprod(:)
               dplitter(p,r,:)   = dplitter(p,r,:)   + tempx * totplitter(:)
               dpsoil(p,r,:)     = dpsoil(p,r,:)     + tempx * totpsoil(:)
               dpsoillab(p,r)    = dpsoillab(p,r)    + tempx * totpsoillab
               dpsoilsorb(p,r)   = dpsoilsorb(p,r)   + tempx * totpsoilsorb
               dpsoilocc(p,r)    = dpsoilocc(p,r)    + tempx * totpsoilocc
            endif
          enddo   ! of "r"
      endif   ! of "fharvw >0.0

      luc%patchfrac_y(p,:)   = luc%patchfrac_x(p,:) + delarea(p,:)

      do d=1,mvmax



      ! here we need to do all biophysical variables and "sumcbal", "sumnbal" and "sumpbal"



      !

      if(luc%patchfrac_x(p,d)+delarea(p,d)>0.0) then
         luc%cplant_y(p,d,leaf)  = (luc%cplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dcplant(p,d,leaf))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%cplant_y(p,d,wood)  = (luc%cplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dcplant(p,d,wood))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%cplant_y(p,d,froot) = (luc%cplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dcplant(p,d,froot)) &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%clabile_y(p,d)      = (luc%clabile_x(p,d)     * luc%patchfrac_x(p,d) + dclabile(p,d))      &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%clitter_y(p,d,metb)  = (luc%clitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dclitter(p,d,metb))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%clitter_y(p,d,str)  = (luc%clitter_x(p,d,str) * luc%patchfrac_x(p,d) + dclitter(p,d,str))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%clitter_y(p,d,cwd)  = (luc%clitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dclitter(p,d,cwd))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%csoil_y(p,d,mic)    = (luc%csoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dcsoil(p,d,mic))    &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%csoil_y(p,d,slow)   = (luc%csoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dcsoil(p,d,slow))   &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%csoil_y(p,d,3)      = (luc%csoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dcsoil(p,d,3))      &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
         luc%cwoodprod_y(p,d,1)  = (luc%cwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,1))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
         luc%cwoodprod_y(p,d,2)  = (luc%cwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,2))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
         luc%cwoodprod_y(p,d,3)  = (luc%cwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dcwoodprod(p,d,3))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))

         if(icycle >1) then
            luc%nplant_y(p,d,leaf)  = (luc%nplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dnplant(p,d,leaf))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nplant_y(p,d,wood)  = (luc%nplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dnplant(p,d,wood))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nplant_y(p,d,froot) = (luc%nplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dnplant(p,d,froot)) &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nlitter_y(p,d,metb) = (luc%nlitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dnlitter(p,d,metb))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nlitter_y(p,d,str)  = (luc%nlitter_x(p,d,str) * luc%patchfrac_x(p,d) + dnlitter(p,d,str))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nlitter_y(p,d,cwd)  = (luc%nlitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dnlitter(p,d,cwd))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nsoil_y(p,d,mic)    = (luc%nsoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dnsoil(p,d,mic))    &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nsoil_y(p,d,slow)   = (luc%nsoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dnsoil(p,d,slow))   &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nsoil_y(p,d,3)      = (luc%nsoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dnsoil(p,d,3))      &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%nsoilmin_y(p,d)     = (luc%nsoilmin_x(p,d)    * luc%patchfrac_x(p,d) + dnsoilmin(p,d))      &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))

            luc%nwoodprod_y(p,d,1)  = (luc%nwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,1))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
            luc%nwoodprod_y(p,d,2)  = (luc%nwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,2))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
            luc%nwoodprod_y(p,d,3)  = (luc%nwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dnwoodprod(p,d,3))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
         endif
         if(icycle >2) then
            luc%pplant_y(p,d,leaf)  = (luc%pplant_x(p,d,leaf) * luc%patchfrac_x(p,d) + dpplant(p,d,leaf))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%pplant_y(p,d,wood)  = (luc%pplant_x(p,d,wood) * luc%patchfrac_x(p,d) + dpplant(p,d,wood))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%pplant_y(p,d,froot) = (luc%pplant_x(p,d,froot)* luc%patchfrac_x(p,d) + dpplant(p,d,froot)) &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%plitter_y(p,d,metb)  = (luc%plitter_x(p,d,metb) * luc%patchfrac_x(p,d) + dplitter(p,d,metb))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%plitter_y(p,d,str)  = (luc%plitter_x(p,d,str) * luc%patchfrac_x(p,d) + dplitter(p,d,str))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%plitter_y(p,d,cwd)  = (luc%plitter_x(p,d,cwd) * luc%patchfrac_x(p,d) + dplitter(p,d,cwd))  &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoil_y(p,d,mic)    = (luc%psoil_x(p,d,mic)   * luc%patchfrac_x(p,d) + dpsoil(p,d,mic))    &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoil_y(p,d,slow)   = (luc%psoil_x(p,d,slow)  * luc%patchfrac_x(p,d) + dpsoil(p,d,slow))   &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoil_y(p,d,3)      = (luc%psoil_x(p,d,3)     * luc%patchfrac_x(p,d) + dpsoil(p,d,3))      &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoillab_y(p,d)     = (luc%psoillab_x(p,d)    * luc%patchfrac_x(p,d) + dpsoillab(p,d))     &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoilsorb_y(p,d)    = (luc%psoilsorb_x(p,d)   * luc%patchfrac_x(p,d) + dpsoilsorb(p,d))    &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))
            luc%psoilocc_y(p,d)     = (luc%psoilocc_x(p,d)    * luc%patchfrac_x(p,d) + dpsoilocc(p,d))     &
                                    /(luc%patchfrac_x(p,d) + delarea(p,d))

            luc%pwoodprod_y(p,d,1)  = (luc%pwoodprod_x(p,d,1) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,1))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
            luc%pwoodprod_y(p,d,2)  = (luc%pwoodprod_x(p,d,2) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,2))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
            luc%pwoodprod_y(p,d,3)  = (luc%pwoodprod_x(p,d,3) * luc%patchfrac_x(p,d) + dpwoodprod(p,d,3))  &
                                    /(luc%patchfrac_x(p,d)   + delarea(p,d))
         endif

      endif   ! luc%patchfrac_x(p,d)+delarea(p,d)>0.0
      enddo   ! of "d"


      ! adding seeding biomass if area > critical value but biomass is too low
      do d=1,mvmax
         if(d>mvtype) then
            ivt = d-mvtype
         else
            ivt = d
         endif

         if(luc%patchfrac_x(p,d)<fracmin.and.luc%patchfrac_y(p,d)>fracmin) then  !newly-born patch
            if(d==1.or.d==2) then
               luc%phase_y(p,d) = 2
            else
               luc%phase_y(p,d) = 0
            endif
         else
            luc%phase_y(p,d) = luc%phase_x(p,d)
         endif

         if(luc%patchfrac_y(p,d) > 0.001.and.sum(luc%cplant_y(p,d,1:3))<fseedling) then
            if(ivt2(ivt)==1) then
               luc%cplant_y(p,d,:) = fseedling * fracgrassseed(:)
            else
               luc%cplant_y(p,d,:) = fseedling * fracwoodseed(:)
            endif
            if(icycle >1) then
               luc%nplant_y(p,d,:) = luc%cplant_y(p,d,:) * casabiome%ratioNCplantmax(ivt,:)
            endif
            if(icycle >2) then
               luc%pplant_y(p,d,:) = luc%nplant_y(p,d,:) /casabiome%ratioNPplantmin(ivt,:)
            endif
         endif
      enddo  ! of "d"

    enddo    ! of "p"


END SUBROUTINE landuse_transitx

 SUBROUTINE landuse_update(luc)                     ! assign "var_y" to "var_x"
 USE landuse_variable
 IMPLICIT NONE
 TYPE(landuse_type) :: luc

    ! general patch variables
    luc%phase_x     = luc%phase_y
    luc%patchfrac_x = luc%patchfrac_y
 
    ! biophysical variables
    ! biogeochemical variables
    luc%cplant_x    = luc%cplant_y
    luc%nplant_x    = luc%nplant_y
    luc%pplant_x    = luc%pplant_y
    luc%clitter_x   = luc%clitter_y
    luc%nlitter_x   = luc%nlitter_y
    luc%plitter_x   = luc%plitter_y
    luc%csoil_x     = luc%csoil_y
    luc%nsoil_x     = luc%nsoil_y
    luc%psoil_x     = luc%psoil_y
    luc%clabile_x   = luc%clabile_y
    luc%nsoilmin_x  = luc%nsoilmin_y
    luc%psoillab_x  = luc%psoillab_y
    luc%psoilsorb_x = luc%psoilsorb_y
    luc%psoilocc_x  = luc%psoilocc_y
    luc%cwoodprod_x = luc%cwoodprod_y
    luc%nwoodprod_x = luc%nwoodprod_y
    luc%pwoodprod_x = luc%pwoodprod_y

 END SUBROUTINE landuse_update

 SUBROUTINE landuse_land2mpx(luc,lucmp,mpx,cstart,cend,nap)
 USE landuse_constant
 USE landuse_variable
 USE landuse_patch
 IMPLICIT NONE
 TYPE(landuse_type)          :: luc
 TYPE(landuse_mp)            :: lucmp
 integer, dimension(mland)   :: cstart,cend,nap
 integer mpx
 integer np,np1,p,q,n,npnew,npold

    npnew=0; npold=0
    do p=1,mland
       do q=1,mvmax
          if(luc%patchfrac_x(p,q)>fracmin.and.luc%patchfrac_y(p,q)>fracmin) then
             npold=npold +1
          endif
          if(luc%patchfrac_y(p,q)>fracmin) then
             npnew = npnew +1
             lucmp%iveg(npnew)      = q
             lucmp%phase(npnew)     = luc%phase_y(p,q)
             lucmp%patchfrac(npnew) = luc%patchfrac_y(p,q)
             ! assign the new biogeochemocal state variables
             lucmp%cplant(npnew,:)    = luc%cplant_y(p,q,:)
             lucmp%clitter(npnew,:)   = luc%clitter_y(p,q,:)
             lucmp%csoil(npnew,:)     = luc%csoil_y(p,q,:)
             lucmp%clabile(npnew)     = luc%clabile_y(p,q)
             lucmp%cwoodprod(npnew,:) = luc%cwoodprod_y(p,q,:)
             if(icycle >1) then
                lucmp%nplant(npnew,:)    = luc%nplant_y(p,q,:)
                lucmp%nlitter(npnew,:)   = luc%nlitter_y(p,q,:)
                lucmp%nsoil(npnew,:)     = luc%nsoil_y(p,q,:)
                lucmp%nsoilmin(npnew)    = luc%nsoilmin_y(p,q)
                lucmp%nwoodprod(npnew,:) = luc%nwoodprod_y(p,q,:)
             endif
             if(icycle >2) then
                lucmp%pplant(npnew,:)    = luc%pplant_y(p,q,:)
                lucmp%plitter(npnew,:)   = luc%plitter_y(p,q,:)
                lucmp%psoil(npnew,:)     = luc%psoil_y(p,q,:)
                lucmp%psoillab(npnew)    = luc%psoillab_y(p,q)
                lucmp%psoilsorb(npnew)   = luc%psoilsorb_y(p,q)
                lucmp%psoilocc(npnew)    = luc%psoilocc_y(p,q)
                lucmp%pwoodprod(npnew,:) = luc%pwoodprod_y(p,q,:)
             endif
          endif
          !update patch_type
       enddo   ! end of "q"
    enddo  ! end of "p"

    print *, 'npnew npold', npnew,npold
    ! creat restart file here
     
    CONTAINS
      real(kind=dp)  function avgpatch(p,q,areax,x2y,x)
      integer p, q, k
      real(kind=dp)                         :: areax
      real(kind=dp), dimension(mvmax,mvmax) :: x2y
      real(kind=dp), dimension(mvmax)       :: x

        avgpatch = 0.0
        if(areax-sum(x2y(1:mvmax,q))+sum(x2y(q,1:mvmax))> 0.001) then
           do k=1,mvmax
              avgpatch = avgpatch + x2y(k,q) * x(k)
           enddo
           avgpatch = (avgpatch + (areax-sum(x2y(1:mvmax,q))* x(q))) &
               /(areax-sum(x2y(1:mvmax,q))+sum(x2y(q,1:mvmax)))
        else
          avgpatch = x(q)   ! check if this works
        endif

         end function avgpatch

 END SUBROUTINE landuse_land2mpx

 SUBROUTINE landuse_checks(mlon,mlat,landmask,luc,arealand)
 ! check mass balance and write output CNP pool sizes for each PFT
 use netcdf
 use landuse_variable
 IMPLICIT NONE
 integer mlon,mlat
 real, parameter      :: xunit = 1.0e-15
 TYPE(landuse_type)   :: luc
! character*120  fpftx
 integer,       dimension(mlon,mlat) :: landmask
 real(kind=dp), dimension(mland)     :: arealand
 real(kind=dp), dimension(mvmax)     :: areapft,areapftx                     
 real(kind=dp), dimension(mvmax)     :: cpland,npland,ppland
 real(kind=dp), dimension(mvmax)     :: clland,nlland,plland
 real(kind=dp), dimension(mvmax)     :: csland,nsland,psland
 real(kind=dp), dimension(mvmax)     :: clabland,nsminland,pslabland,pssorbland,psoccland
 real(kind=dp), dimension(mvmax)     :: cwoodland,nwoodland,pwoodland
 integer n,v
 real(kind=dp)  totalc,totaln,totalp,totarea

 ! local variables
 real(kind=dp),   dimension(:,:),     allocatable   :: fracpft
 real(kind=dp),   dimension(:,:,:),   allocatable   :: primary
 real(kind=dp),   dimension(:,:,:),   allocatable   :: secondary
 integer ok,ncid1,varxid,i,j,k,np,nland

    allocate(primary(mlon,mlat,mvtype))
    allocate(secondary(mlon,mlat,mvtype))
    allocate(fracpft(mland,mvmax))

    ok = nf90_open(fpft,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"primary",varxid)
    ok = nf90_get_var(ncid1,varxid,primary)
    ok = nf90_inq_varid(ncid1,"secondary",varxid)
    ok = nf90_get_var(ncid1,varxid,secondary)
    ok = nf90_close(ncid1)

    nland = 0
    np    = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          nland = nland +1
          do k=1,4
             fracpft(nland,k)        = primary(i,j,k)
             fracpft(nland,mvtype+k) = secondary(i,j,k)
          enddo
          do k=5,mvtype
             fracpft(nland,k) = primary(i,j,k) + secondary(i,j,k)
          enddo
       endif
    enddo
    enddo
91  format(3(i5,2x),21(f7.5,1x))

    areapft=0.0;     areapftx=0.0
    totalc=0.0;      totaln=0.0;      totalp=0.0;     totarea=0.0
    cpland = 0.0;    npland = 0.0;    ppland=0.0
    clland = 0.0;    nlland = 0.0;    plland = 0.0
    csland = 0.0;    nsland = 0.0;    psland = 0.0
    nsminland=0.0;   pslabland=0.0;   pssorbland=0.0; psoccland = 0.0
    cwoodland = 0.0; nwoodland = 0.0; pwoodland = 0.0

    do n=1,mland
    do v=1,mvmax
    if(luc%patchfrac_y(n,v) > 0.0) then

       areapftx(v)= areapftx(v) + arealand(n) * fracpft(n,v)
       areapft(v)= areapft(v)+ arealand(n) *luc%patchfrac_y(n,v)
       cpland(v) = cpland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%cplant_y(n,v,1:3))
       npland(v) = npland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%nplant_y(n,v,1:3))
       ppland(v) = ppland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%pplant_y(n,v,1:3))
       clland(v) = clland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%clitter_y(n,v,1:3))
       nlland(v) = nlland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%nlitter_y(n,v,1:3))
       plland(v) = plland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%plitter_y(n,v,1:3))
       csland(v) = csland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%csoil_y(n,v,1:3))
       nsland(v) = nsland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%nsoil_y(n,v,1:3))
       psland(v) = psland(v) + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%psoil_y(n,v,1:3))
 
       clabland(v)   = clabland(v)   + arealand(n) *luc%patchfrac_y(n,v) * luc%clabile_y(n,v)
       nsminland(v)  = nsminland(v)  + arealand(n) *luc%patchfrac_y(n,v) * luc%nsoilmin_y(n,v)
       pslabland(v)  = pslabland(v)  + arealand(n) *luc%patchfrac_y(n,v) * luc%psoillab_y(n,v)
       pssorbland(v) = pssorbland(v) + arealand(n) *luc%patchfrac_y(n,v) * luc%psoilsorb_y(n,v)
       psoccland(v)  = psoccland(v)  + arealand(n) *luc%patchfrac_y(n,v) * luc%psoilocc_y(n,v)
       cwoodland(v)  = cwoodland(v)  + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%cwoodprod_y(n,v,1:3))
       nwoodland(v)  = nwoodland(v)  + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%nwoodprod_y(n,v,1:3))
       pwoodland(v)  = pwoodland(v)  + arealand(n) *luc%patchfrac_y(n,v) * sum(luc%pwoodprod_y(n,v,1:3))
    endif
    enddo
    enddo

    areapft= areapft *1000.0*xunit
    areapftx=areapftx*1000.0*xunit
    cpland = cpland * xunit
    npland = npland * xunit
    ppland = ppland * xunit
    clland = clland * xunit
    nlland = nlland * xunit
    plland = plland * xunit
    csland = csland * xunit
    nsland = nsland * xunit
    psland = psland * xunit
    clabland  = clabland * xunit
    nsminland = nsminland * xunit
    pslabland = pslabland * xunit
    pssorbland = pssorbland * xunit
    psoccland = psoccland * xunit
    cwoodland = cwoodland * xunit
    nwoodland = nwoodland * xunit
    pwoodland = pwoodland * xunit
    totarea= sum(areapft)
    totalc = sum(cpland+clland+csland+clabland+cwoodland)
    totaln = sum(npland+nlland+nsland+nsminland+nwoodland)
    totalp = sum(ppland+plland+psland+pslabland+pssorbland+psoccland+pwoodland)

    do v=1,mvmax
       write(21,201) v, areapft(v),areapftx(v),cpland(v)+clland(v)+csland(v)+clabland(v)+cwoodland(v),     &
                        npland(v)+nlland(v)+nsland(v)+nsminland(v)+nwoodland(v),                           &
                        ppland(v)+plland(v)+psland(v)+pslabland(v)+pssorbland(v)+psoccland(v)+pwoodland(v)
    enddo
    write(21,202) totarea,totalc,totaln,totalp

201 format(i3,2x,20(f10.5,2x))
202 format(20(f10.5,2x))
   
    deallocate(fracpft)
    deallocate(primary)
    deallocate(secondary)
 END SUBROUTINE landuse_checks
