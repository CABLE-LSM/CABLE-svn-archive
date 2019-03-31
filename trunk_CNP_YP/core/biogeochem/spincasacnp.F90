SUBROUTINE spincasacnp( dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
     casaflux,casamet,casabal,phen,POP,climate,LALLOC )


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

  integer nptx,nvt,kloop

  ktauday=int(24.0*3600.0/dels)
  nday=(kend-kstart+1)/ktauday

     allocate(avg_cleaf2met(mp), avg_cleaf2str(mp), avg_croot2met(mp), avg_croot2str(mp), avg_cwood2cwd(mp), &
              avg_nleaf2met(mp), avg_nleaf2str(mp), avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp), &
              avg_pleaf2met(mp), avg_pleaf2str(mp), avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp), &
               avg_cgpp(mp),      avg_cnpp(mp),      avg_nuptake(mp),   avg_puptake(mp),                     &
              avg_xnplimit(mp),  avg_xkNlimiting(mp), avg_xklitter(mp), avg_xksoil(mp),                      &
              avg_rationcsoilmic(mp),avg_rationcsoilslow(mp),avg_rationcsoilpass(mp),                                 &!chris 12/oct/2012 for spin up casa
              avg_nsoilmin(mp),  avg_psoillab(mp),    avg_psoilsorb(mp), avg_psoilocc(mp))

  OPEN(91, file=fcnpspin)
  read(91,*) myearspin

  ! compute the mean fluxes and residence time of each carbon pool
     avg_cleaf2met=0.0; avg_cleaf2str=0.0; avg_croot2met=0.0; avg_croot2str=0.0; avg_cwood2cwd=0.0
     avg_nleaf2met=0.0; avg_nleaf2str=0.0; avg_nroot2met=0.0; avg_nroot2str=0.0; avg_nwood2cwd=0.0
     avg_pleaf2met=0.0; avg_pleaf2str=0.0; avg_proot2met=0.0; avg_proot2str=0.0; avg_pwood2cwd=0.0
     avg_cgpp=0.0;      avg_cnpp=0.0;      avg_nuptake=0.0;   avg_puptake=0.0
     avg_xnplimit=0.0;  avg_xkNlimiting=0.0; avg_xklitter=0.0; avg_xksoil=0.0
     avg_nsoilmin=0.0;  avg_psoillab=0.0;    avg_psoilsorb=0.0; avg_psoilocc=0.0
     avg_rationcsoilmic=0.0;avg_rationcsoilslow=0.0;avg_rationcsoilpass=0.0
  do nyear=1,myearspin
     read(91,901) ncfile
     call read_casa_dump(ncfile,casamet,casaflux,ktau,kend)
901  format(A99)
    do idoy=1,mdyear
       ktau=(idoy-1)*ktauday +1
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

       CALL biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen,xnplimit,xkNlimiting,xklitter,xksoil,xkleaf,xkleafcold,xkleafdry,&
                    cleaf2met,cleaf2str,croot2met,croot2str,cwood2cwd,         &
                    nleaf2met,nleaf2str,nroot2met,nroot2str,nwood2cwd,         &
                    pleaf2met,pleaf2str,proot2met,proot2str,pwood2cwd)

        WHERE(xkNlimiting .eq. 0)  !Chris Lu 4/June/2012
           xkNlimiting = 0.001
        END WHERE
        !nptx=8173

        avg_cleaf2met = avg_cleaf2met + cleaf2met
        avg_cleaf2str = avg_cleaf2str + cleaf2str
        avg_croot2met = avg_croot2met + croot2met
        avg_croot2str = avg_croot2str + croot2str
        avg_cwood2cwd = avg_cwood2cwd + cwood2cwd

        avg_nleaf2met = avg_nleaf2met + nleaf2met
        avg_nleaf2str = avg_nleaf2str + nleaf2str
        avg_nroot2met = avg_nroot2met + nroot2met
        avg_nroot2str = avg_nroot2str + nroot2str
        avg_nwood2cwd = avg_nwood2cwd + nwood2cwd

        avg_pleaf2met = avg_pleaf2met + pleaf2met
        avg_pleaf2str = avg_pleaf2str + pleaf2str
        avg_proot2met = avg_proot2met + proot2met
        avg_proot2str = avg_proot2str + proot2str
        avg_pwood2cwd = avg_pwood2cwd + pwood2cwd

        avg_cgpp      = avg_cgpp      + casaflux%cgpp
        avg_cnpp      = avg_cnpp      + casaflux%cnpp
        avg_nuptake   = avg_nuptake   + casaflux%Nminuptake
        avg_puptake   = avg_puptake   + casaflux%Plabuptake

        avg_xnplimit    = avg_xnplimit    + xnplimit
        avg_xkNlimiting = avg_xkNlimiting + xkNlimiting
        avg_xklitter    = avg_xklitter    + xklitter
        avg_xksoil      = avg_xksoil      + xksoil

        avg_nsoilmin    = avg_nsoilmin    + casapool%nsoilmin
        avg_psoillab    = avg_psoillab    + casapool%psoillab
        avg_psoilsorb   = avg_psoilsorb   + casapool%psoilsorb
        avg_psoilocc    = avg_psoilocc    + casapool%psoilocc

        avg_rationcsoilmic  = avg_rationcsoilmic  + casapool%ratioNCsoilnew(:,mic)
        avg_rationcsoilslow = avg_rationcsoilslow + casapool%ratioNCsoilnew(:,slow)
        avg_rationcsoilpass = avg_rationcsoilpass + casapool%ratioNCsoilnew(:,pass)
    enddo
    enddo

    CLOSE(91)

    avg_cleaf2met = avg_cleaf2met/real(nday*myearspin)
    avg_cleaf2str = avg_cleaf2str/real(nday*myearspin)
    avg_croot2met = avg_croot2met/real(nday*myearspin)
    avg_croot2str = avg_croot2str/real(nday*myearspin)
    avg_cwood2cwd = avg_cwood2cwd/real(nday*myearspin)

    avg_nleaf2met = avg_nleaf2met/real(nday*myearspin)
    avg_nleaf2str = avg_nleaf2str/real(nday*myearspin)
    avg_nroot2met = avg_nroot2met/real(nday*myearspin)
    avg_nroot2str = avg_nroot2str/real(nday*myearspin)
    avg_nwood2cwd = avg_nwood2cwd/real(nday*myearspin)

    avg_pleaf2met = avg_pleaf2met/real(nday*myearspin)
    avg_pleaf2str = avg_pleaf2str/real(nday*myearspin)
    avg_proot2met = avg_proot2met/real(nday*myearspin)
    avg_proot2str = avg_proot2str/real(nday*myearspin)
    avg_pwood2cwd = avg_pwood2cwd/real(nday*myearspin)

    avg_cgpp      = avg_cgpp/real(nday*myearspin)
    avg_cnpp      = avg_cnpp/real(nday*myearspin)
    avg_nuptake   = avg_nuptake/real(nday*myearspin)
    avg_puptake   = avg_puptake/real(nday*myearspin)

    avg_xnplimit    = avg_xnplimit/real(nday*myearspin)
    avg_xkNlimiting = avg_xkNlimiting/real(nday*myearspin)
    avg_xklitter    = avg_xklitter/real(nday*myearspin)
    avg_xksoil      = avg_xksoil/real(nday*myearspin)

    avg_nsoilmin    = avg_nsoilmin/real(nday*myearspin)
    avg_psoillab    = avg_psoillab/real(nday*myearspin)
    avg_psoilsorb   = avg_psoilsorb/real(nday*myearspin)
    avg_psoilocc    = avg_psoilocc/real(nday*myearspin)

    avg_rationcsoilmic  = avg_rationcsoilmic  /real(nday*myearspin)
    avg_rationcsoilslow = avg_rationcsoilslow /real(nday*myearspin)
    avg_rationcsoilpass = avg_rationcsoilpass /real(nday*myearspin)

    call analyticpool(kend,veg,soil,casabiome,casapool,                                          &
                          casaflux,casamet,casabal,phen,                                         &
                          avg_cleaf2met,avg_cleaf2str,avg_croot2met,avg_croot2str,avg_cwood2cwd, &
                          avg_nleaf2met,avg_nleaf2str,avg_nroot2met,avg_nroot2str,avg_nwood2cwd, &
                          avg_pleaf2met,avg_pleaf2str,avg_proot2met,avg_proot2str,avg_pwood2cwd, &
                          avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake,                          &
                          avg_xnplimit,avg_xkNlimiting,avg_xklitter,avg_xksoil,                  &
                          avg_ratioNCsoilmic,avg_ratioNCsoilslow,avg_ratioNCsoilpass,            &
                          avg_nsoilmin,avg_psoillab,avg_psoilsorb,avg_psoilocc)

END SUBROUTINE spincasacnp
