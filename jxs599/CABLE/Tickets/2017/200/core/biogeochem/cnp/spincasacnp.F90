SUBROUTINE spincasacnp(fcnpspin,dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
                       casaflux,casamet,casabal,phen)
  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  CHARACTER(LEN=99), INTENT(IN)  :: fcnpspin
  REAL,    INTENT(IN)    :: dels
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN)    :: mloop
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen

  ! local variables
  real,      dimension(:), allocatable, save  :: avg_cleaf2met, avg_cleaf2str, avg_croot2met, avg_croot2str, avg_cwood2cwd
  real,      dimension(:), allocatable, save  :: avg_nleaf2met, avg_nleaf2str, avg_nroot2met, avg_nroot2str, avg_nwood2cwd
  real,      dimension(:), allocatable, save  :: avg_pleaf2met, avg_pleaf2str, avg_proot2met, avg_proot2str, avg_pwood2cwd
  real,      dimension(:), allocatable, save  :: avg_cgpp,      avg_cnpp,      avg_nuptake,   avg_puptake
  real,      dimension(:), allocatable, save  :: avg_nsoilmin,  avg_psoillab,  avg_psoilsorb, avg_psoilocc
  real,      dimension(:), allocatable, save  :: avg_ratioNCsoilmic,  avg_ratioNCsoilslow,  avg_ratioNCsoilpass !chris 12/oct/2012 for spin up casa
  real(r_2), dimension(:), allocatable, save  :: avg_xnplimit,  avg_xkNlimiting,avg_xklitter, avg_xksoil

  ! local variables
  INTEGER                  :: myearspin,nyear, nloop1
  CHARACTER(LEN=99)        :: ncfile
  INTEGER                  :: ktau,ktauday,nday,idoy,ktaux,ktauy,nloop
  INTEGER, save            :: ndays
  real,      dimension(mp)      :: cleaf2met, cleaf2str, croot2met, croot2str, cwood2cwd
  real,      dimension(mp)      :: nleaf2met, nleaf2str, nroot2met, nroot2str, nwood2cwd
  real,      dimension(mp)      :: pleaf2met, pleaf2str, proot2met, proot2str, pwood2cwd
  real,      dimension(mp)      :: xcgpp,     xcnpp,     xnuptake,  xpuptake
  real,      dimension(mp)      :: xnsoilmin, xpsoillab, xpsoilsorb,xpsoilocc
  real(r_2), dimension(mp)      :: xnplimit,  xkNlimiting, xklitter, xksoil,xkleaf, xkleafcold, xkleafdry

  ! more variables to store the spinup pool size over the last 10 loops. Added by Yp Wang 30 Nov 2012
  real,      dimension(10,mvtype,mplant)  :: bmcplant,  bmnplant,  bmpplant
  real,      dimension(10,mvtype,mlitter) :: bmclitter, bmnlitter, bmplitter
  real,      dimension(10,mvtype,msoil)   :: bmcsoil,   bmnsoil,   bmpsoil
  real,      dimension(10,mvtype)         :: bmnsoilmin,bmpsoillab,bmpsoilsorb, bmpsoilocc
  real,      dimension(mvtype)           :: bmarea
  integer nptx,nvt,kloop


  !nptx = 59621

  bmcplant  = 0.0;  bmnplant  = 0.0; bmpplant  = 0.0
  bmclitter = 0.0;  bmnlitter = 0.0; bmplitter = 0.0
  bmcsoil   = 0.0;  bmnsoil   = 0.0; bmpsoil   = 0.0
  bmnsoilmin  = 0.0;  bmpsoillab  = 0.0; bmpsoilsorb = 0.0;  bmpsoilocc = 0.0

  call totcnppools(1,veg,casamet,casapool, &
                bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter,  &
                bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb, &
                bmpsoilocc,bmarea)
!  call spinanalytic(fcnpspin,dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
!                       casaflux,casamet,casabal,phen)

  nloop1= max(1,mloop-4)

  DO nloop=2,mloop

  OPEN(91,file=fcnpspin)
  read(91,*) myearspin
  DO nyear=1,myearspin
     read(91,901) ncfile
     call read_casa_dump(ncfile,casamet,casaflux,ktau,kend)

    DO idoy=1,mdyear
      ktauy=idoy*ktauday
       casamet%tairk(:)       = casamet%Tairkspin(:,idoy)
       casamet%tsoil(:,1)     = casamet%Tsoilspin_1(:,idoy)
       casamet%tsoil(:,2)     = casamet%Tsoilspin_2(:,idoy)
       casamet%tsoil(:,3)     = casamet%Tsoilspin_3(:,idoy)
       casamet%tsoil(:,4)     = casamet%Tsoilspin_4(:,idoy)
       casamet%tsoil(:,5)     = casamet%Tsoilspin_5(:,idoy)
       casamet%tsoil(:,6)     = casamet%Tsoilspin_6(:,idoy)
       casamet%moist(:,1)     = casamet%moistspin_1(:,idoy)
       casamet%moist(:,2)     = casamet%moistspin_2(:,idoy)
       casamet%moist(:,3)     = casamet%moistspin_3(:,idoy)
       casamet%moist(:,4)     = casamet%moistspin_4(:,idoy)
       casamet%moist(:,5)     = casamet%moistspin_5(:,idoy)
       casamet%moist(:,6)     = casamet%moistspin_6(:,idoy)
       casaflux%cgpp(:)       = casamet%cgppspin(:,idoy)
       casaflux%crmplant(:,1) = casamet%crmplantspin_1(:,idoy)
       casaflux%crmplant(:,2) = casamet%crmplantspin_2(:,idoy)
       casaflux%crmplant(:,3) = casamet%crmplantspin_3(:,idoy)
       call biogeochem(ktauy,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                      casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
                      cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                      nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                      pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)
!    write(89,891) idoy,nyear,nloop,nptx,&
!                  casaflux%cgpp(nptx),casaflux%cnpp(nptx),casaflux%crmplant(nptx,:),casaflux%Crgplant(nptx),casaflux%Crsoil(nptx), &
!                  casaflux%clabloss(nptx),casapool%dClabiledt(nptx),casaflux%fracClabile(nptx),                                  &
!                  casapool%Cplant(nptx,:), casapool%clitter(nptx,:), casapool%csoil(nptx,:),casabal%cbalance(nptx)
    ENDDO   ! end of idoy

  ENDDO   ! end of nyear
  close(91)

891 format(4(i6,2x),100(f9.2,1x))

    if(nloop>=nloop1) then
       print *, 'kloop =', 2+nloop-nloop1, nloop,nloop1
       call totcnppools(2+nloop-nloop1,veg,casamet,casapool, &
            bmcplant,bmnplant,bmpplant,bmclitter,bmnlitter,bmplitter, &
            bmcsoil,bmnsoil,bmpsoil,bmnsoilmin,bmpsoillab,bmpsoilsorb, &
            bmpsoilocc,bmarea)
    endif

  ENDDO     ! end of nloop

  ! write the last five loop pool size by PFT type
  open(92,file='cnpspinlast5.txt', position='append')
  write(92,*) 'myearspin =', myearspin
  write(92,921)
921 format('PFT total area in 10**12 m2', f12.4)
  do nvt=1,mvtype
     write(92,*) bmarea(nvt)
  enddo

  do nvt=1,mvtype
  if(bmarea(nvt) >0.0) then
     do kloop=1,5
        write(92,922) nvt, bmcplant(kloop,nvt,:),bmclitter(kloop,nvt,:),bmcsoil(kloop,nvt,:) 
     enddo
     if (icycle >1) then 
        do kloop=1,5
           write(92,922) nvt, bmnplant(kloop,nvt,:),bmnlitter(kloop,nvt,:),bmnsoil(kloop,nvt,:), bmnsoilmin(kloop,nvt) 
        enddo
     endif

     if(icycle >2) then
        do kloop=1,5
           write(92,922) nvt, bmpplant(kloop,nvt,:),bmplitter(kloop,nvt,:),bmpsoil(kloop,nvt,:),  &
                              bmpsoillab(kloop,nvt), bmpsoilsorb(kloop,nvt), bmpsoilocc(kloop,nvt)
        enddo
     endif
  endif 
  enddo

901  format(A99)
922 format(i4,20(f10.4,2x))
    CLOSE(92)   
151 FORMAT(i6,100(f12.5,2x))
END SUBROUTINE spincasacnp

