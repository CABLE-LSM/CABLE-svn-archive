MODULE casa_spin_mod
  USE cable_def_types_mod
  USE casadimension
  USE casavariable

  IMPLICIT NONE

  INTEGER, PARAMETER :: num_vars=7
  INTEGER, PARAMETER :: num_dims=3

  CHARACTER(LEN=*), DIMENSION(num_vars), PARAMETER :: &
        var_name =  (/  "latitude     ", &
                        "longitude    ", &
                        "casamet_tairk", &
                        "tsoil        ", &
                        "moist        ", &
                        "cgpp         ", &
                        "crmplant     " /)
  INTEGER, DIMENSION(num_vars)  :: &
        varID      ! (1) lat (2) lon (3) tair, etc

  CHARACTER(LEN=*), DIMENSION(num_dims), PARAMETER :: &
        dim_name =  (/ "tile", &
                       "soil", &
                       "time" /)
  INTEGER, DIMENSION(num_dims)  :: &
        dimID    & ! (1) tile, (2) soil, (3) time
        dim_len

  INTEGER  :: mloop,        &
              endday,       &
              ncid_spin,    &
              ncok

  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: avg_cleaf2met, avg_cleaf2str,  &
                            avg_croot2met, avg_croot2str, avg_cwood2cwd,  &
                            avg_nleaf2met, avg_nleaf2str,                 &
                            avg_nroot2met, avg_nroot2str, avg_nwood2cwd,  &
                            avg_pleaf2met, avg_pleaf2str,                 &
                            avg_proot2met, avg_proot2str, avg_pwood2cwd,  &
                            avg_cgpp, avg_cnpp, avg_nuptake, avg_puptake, &
                            avg_nsoilmin,  avg_psoillab,                  &
                            avg_psoilsorb, avg_psoilocc
  REAL(r_2), DIMENSION(:), ALLOCATABLE, SAVE :: avg_xnplimit, avg_xkNlimiting, &
                                                avg_xklitter, avg_xksoil

CONTAINS

  SUBROUTINE open_casa_dump( ncfile )
    USE netcdf
    USE cable_diag_module, ONLY : get_var_nc, stderr_nc

    CHARACTER(LEN=99), INTENT(IN) :: ncfile

    ncok = NF90_OPEN(ncfile, NF90_nowrite, ncid_spin)
    IF (ncok /= NF90_noerr ) CALL stderr_nc('Error opening ', ncfile)
    ncok = NF90_INQ_DIMID(ncid_spin,dim_name(3),dimID(3))
    IF (ncok /= NF90_noerr ) CALL stderr_nc('Error inquiring time dimension.')
    ncok = NF90_INQUIRE_DIMENSION(ncid_spin,dimID(3),LEN=endday)
    IF (ncok /= NF90_noerr ) CALL stderr_nc('Error getting time dimension.')

  END SUBROUTINE open_casa_dump


  SUBROUTINE read_casa_dump( casamet, casaflux, iday )
    USE netcdf
    USE casa_cnp_module
    USE cable_diag_module, only : get_var_nc, stderr_nc

    TYPE (casa_flux), INTENT(INOUT) :: casaflux
    TYPE (casa_met),  INTENT(INOUT) :: casamet
    INTEGER, INTENT(IN) :: iday

    CALL get_var_nc(ncid_spin, var_name(3), casamet%tairk,     iday, endday )
    CALL get_var_nc(ncid_spin, var_name(4), casamet%tsoil,     iday, endday )
    CALL get_var_nc(ncid_spin, var_name(5), casamet%moist,     iday, endday )
    CALL get_var_nc(ncid_spin, var_name(6), casaflux%cgpp,     iday, endday )
    CALL get_var_nc(ncid_spin, var_name(7), casaflux%crmplant, iday, endday )

  END SUBROUTINE read_casa_dump


  SUBROUTINE spin_init
    
    ALLOCATE( avg_cleaf2met(mp), avg_cleaf2str(mp),                          &
              avg_croot2met(mp), avg_croot2str(mp), avg_cwood2cwd(mp),       &
              avg_nleaf2met(mp), avg_nleaf2str(mp),                          &
              avg_nroot2met(mp), avg_nroot2str(mp), avg_nwood2cwd(mp),       &
              avg_pleaf2met(mp), avg_pleaf2str(mp),                          &
              avg_proot2met(mp), avg_proot2str(mp), avg_pwood2cwd(mp),       &
              avg_cgpp(mp), avg_cnpp(mp), avg_nuptake(mp), avg_puptake(mp),  &
              avg_xnplimit(mp),  avg_xkNlimiting(mp),                        &
              avg_xklitter(mp),  avg_xksoil(mp),                             &
              avg_nsoilmin(mp),  avg_psoillab(mp),                           &
              avg_psoilsorb(mp), avg_psoilocc(mp) )

    avg_cleaf2met=0.0; avg_cleaf2str=0.0
    avg_croot2met=0.0; avg_croot2str=0.0; avg_cwood2cwd=0.0
    avg_nleaf2met=0.0; avg_nleaf2str=0.0
    avg_nroot2met=0.0; avg_nroot2str=0.0; avg_nwood2cwd=0.0
    avg_pleaf2met=0.0; avg_pleaf2str=0.0
    avg_proot2met=0.0; avg_proot2str=0.0; avg_pwood2cwd=0.0
    avg_cgpp=0.0;      avg_cnpp=0.0;      avg_nuptake=0.0;   avg_puptake=0.0
    avg_xnplimit=0.0;  avg_xkNlimiting=0.0
    avg_xklitter=0.0;  avg_xksoil=0.0
    avg_nsoilmin=0.0;  avg_psoillab=0.0
    avg_psoilsorb=0.0; avg_psoilocc=0.0

  END SUBROUTINE spin_init


  SUBROUTINE spin_update( casaflux, casapool )
    TYPE (casa_flux), INTENT(INOUT) :: casaflux
    TYPE (casa_pool), INTENT(INOUT) :: casapool

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

  END SUBROUTINE spin_update


  SUBROUTINE spin_averaging

    avg_cleaf2met = avg_cleaf2met/REAL(endday)
    avg_cleaf2str = avg_cleaf2str/REAL(endday)
    avg_croot2met = avg_croot2met/REAL(endday)
    avg_croot2str = avg_croot2str/REAL(endday)
    avg_cwood2cwd = avg_cwood2cwd/REAL(endday)

    avg_nleaf2met = avg_nleaf2met/REAL(endday)
    avg_nleaf2str = avg_nleaf2str/REAL(endday)
    avg_nroot2met = avg_nroot2met/REAL(endday)
    avg_nroot2str = avg_nroot2str/REAL(endday)
    avg_nwood2cwd = avg_nwood2cwd/REAL(endday)

    avg_pleaf2met = avg_pleaf2met/REAL(endday)
    avg_pleaf2str = avg_pleaf2str/REAL(endday)
    avg_proot2met = avg_proot2met/REAL(endday)
    avg_proot2str = avg_proot2str/REAL(endday)
    avg_pwood2cwd = avg_pwood2cwd/REAL(endday)

    avg_cgpp      = avg_cgpp/REAL(endday)
    avg_cnpp      = avg_cnpp/REAL(endday)
    avg_nuptake   = avg_nuptake/REAL(endday)
    avg_puptake   = avg_puptake/REAL(endday)

    avg_xnplimit    = avg_xnplimit/REAL(endday)
    avg_xkNlimiting = avg_xkNlimiting/REAL(endday)
    avg_xklitter    = avg_xklitter/REAL(endday)
    avg_xksoil      = avg_xksoil/REAL(endday)

    avg_nsoilmin    = avg_nsoilmin/REAL(endday)
    avg_psoillab    = avg_psoillab/REAL(endday)
    avg_psoilsorb   = avg_psoilsorb/REAL(endday)
    avg_psoilocc    = avg_psoilocc/REAL(endday)

  END SUBROUTINE spin_averaging


  SUBROUTINE analyticPool( veg, soil, casabiome, casapool,       &
                           casaflux, casamet, casabal, phen )
    USE carbon_module
    USE casaparm
    USE phenvariable     !  not used ??

    IMPLICIT NONE
    TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
    TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters  
    TYPE (casa_biome),            INTENT(INOUT) :: casabiome
    TYPE (casa_pool),             INTENT(INOUT) :: casapool
    TYPE (casa_flux),             INTENT(INOUT) :: casaflux
    TYPE (casa_met),              INTENT(INOUT) :: casamet
    TYPE (casa_balance),          INTENT(INOUT) :: casabal
    TYPE (phen_variable),         INTENT(INOUT) :: phen      !  not used ??

    ! local variables
    REAL(r_2), DIMENSION(mso) :: Psorder,Pweasoil,xPsoil50
    REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
    REAL(r_2), DIMENSION(:), ALLOCATABLE :: totPsoil
    INTEGER  npt

    ! Soiltype     soilnumber soil P(g P/m2)
    ! Alfisol     1       61.3
    ! Andisol     2       103.9
    ! Aridisol    3       92.8
    ! Entisol     4       136.9
    ! Gellisol    5       98.2
    ! Histosol    6       107.6
    ! Inceptisol  7       84.1
    ! Mollisol    8       110.1
    ! Oxisol      9       35.4    
    ! Spodosol    10      41.0    
    ! Ultisol     11      51.5    
    ! Vertisol    12      190.6
    DATA Psorder/61.3,103.9,92.8,136.9,98.2,107.6,  &
                 84.1,110.1,35.4,41.0,51.5,190.6/
    DATA Pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,   &
                  0.008,0.007,0.006,0.005,0.004,0.003/   ! not used ??
    DATA fracPlab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
    DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/  ! not used ??
    DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
    DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/  ! not used ??
    DATA xPsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/
  

    ALLOCATE(totPsoil(mp))

    ! To keep klitter non-zero or else divide by zero
    WHERE (avg_xkNlimiting < 0.001)  avg_xkNlimiting = 0.001

    ! compute the mean litter input in g(C, N and P)/day from plant pools
    casaflux%fromLtoS = 0.0
    casaflux%fromStoS = 0.0

    casabal%sumcbal(:)   = 0.0
    casabal%sumnbal(:)   = 0.0
    casabal%sumpbal(:)   = 0.0

    DO npt=1,mp
    IF(casamet%iveg2(npt)/=icewater.and.avg_cnpp(npt) > 0.0) THEN
      casaflux%fromLtoS(npt,mic,metb)   = 0.45
                                          ! metb -> mic
      casaflux%fromLtoS(npt,mic,str)   = 0.45*(1.0-casabiome%fracLigninplant(veg%iveg(npt),leaf))
                                          ! str -> mic
      casaflux%fromLtoS(npt,slow,str)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),leaf)
                                          ! str -> slow
      casaflux%fromLtoS(npt,mic,cwd)   = 0.40*(1.0 - casabiome%fracLigninplant(veg%iveg(npt),wood))
                                          ! CWD -> fmic
      casaflux%fromLtoS(npt,slow,cwd)  = 0.7 * casabiome%fracLigninplant(veg%iveg(npt),wood)
                                          ! CWD -> slow
!! set the following two backflow to set (see Bolker 199x)
!      casaflux%fromStoS(npt,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(npt))
!      casaflux%fromStoS(npt,mic,pass)  = 0.45

      casaflux%fromStoS(npt,slow,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.997 - 0.032*soil%clay(npt))
      casaflux%fromStoS(npt,pass,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.003 + 0.032*soil%clay(npt))
      casaflux%fromStoS(npt,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(npt) )

      casaflux%klitter(npt,metb) = avg_xkNlimiting(npt) * avg_xklitter(npt)*casabiome%litterrate(veg%iveg(npt),metb)
      casaflux%klitter(npt,str)  = avg_xkNlimiting(npt) * avg_xklitter(npt)*casabiome%litterrate(veg%iveg(npt),str)&
                               * exp(-3.0*casabiome%fracLigninplant(veg%iveg(npt),leaf))
      casaflux%klitter(npt,cwd)  = avg_xkNlimiting(npt) * avg_xklitter(npt)*casabiome%litterrate(veg%iveg(npt),cwd)

      casaflux%ksoil(npt,mic)    = avg_xksoil(npt)*casabiome%soilrate(veg%iveg(npt),mic)  &
                               * (1.0 - 0.75 *(soil%silt(npt)+soil%clay(npt)))
      casaflux%ksoil(npt,slow)   = avg_xksoil(npt) * casabiome%soilrate(veg%iveg(npt),slow)
      casaflux%ksoil(npt,pass)   = avg_xksoil(npt) * casabiome%soilrate(veg%iveg(npt),pass)

      IF(veg%iveg(npt)==cropland) THEN     ! for cultivated land type
        casaflux%ksoil(npt,mic)  = casaflux%ksoil(npt,mic) * 1.25
        casaflux%ksoil(npt,slow) = casaflux%ksoil(npt,slow)* 1.5
        casaflux%ksoil(npt,pass) = casaflux%ksoil(npt,pass)* 1.5
      ENDIF

    ENDIF
    ENDDO


    DO npt=1,mp
    IF(casamet%iveg2(npt)/=icewater.and.avg_cnpp(npt) > 0.0) THEN
      ! compute steady-state litter and soil C pool sizes
      casapool%clitter(npt,metb) = (avg_cleaf2met(npt)+avg_croot2met(npt))/casaflux%klitter(npt,metb)
      casapool%clitter(npt,str) = (avg_cleaf2str(npt)+avg_croot2str(npt))/casaflux%klitter(npt,str)
      casapool%clitter(npt,cwd) = (avg_cwood2cwd(npt))/casaflux%klitter(npt,cwd)
      casapool%csoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                 +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                 +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                               /casaflux%ksoil(npt,mic)
      casapool%csoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                 + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                 + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                 + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                /casaflux%ksoil(npt,slow)
      casapool%csoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                  +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                /casaflux%ksoil(npt,pass)
      IF(icycle <=1) THEN
        casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) * casapool%clitter(npt,:)
        casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   * casapool%Csoil(npt,:)
        casapool%nsoilmin(npt) = 2.0
        casabal%sumnbal(npt)   = 0.0
      ELSE
        ! compute steady-state litter and soil N pool sizes
        casapool%nlitter(npt,metb) = (avg_nleaf2met(npt)+avg_nroot2met(npt))/casaflux%klitter(npt,metb)
        casapool%nlitter(npt,str) = (avg_nleaf2str(npt)+avg_nroot2str(npt))/casaflux%klitter(npt,str)
        casapool%nlitter(npt,cwd) = (avg_nwood2cwd(npt))/casaflux%klitter(npt,cwd)

        casapool%nsoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * casapool%ratioNCsoil(npt,mic)/casaflux%ksoil(npt,mic)
        casapool%nsoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * casapool%ratioNCsoil(npt,slow)/casaflux%ksoil(npt,slow)
        casapool%nsoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   * casapool%ratioNCsoil(npt,pass)/casaflux%ksoil(npt,pass)
        casapool%Nsoilmin(npt)    = avg_nsoilmin(npt)

      ENDIF

      IF (icycle<=2) THEN
        totPsoil(npt)          = Psorder(casamet%isorder(npt)) *xPsoil50(casamet%isorder(npt))
        casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:)  * casapool%clitter(npt,:)
        casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)    * casapool%Csoil(npt,:)
        casapool%psoillab(npt) = totPsoil(npt) *fracPlab(casamet%isorder(npt))
        casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                    /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
        casapool%psoilocc(npt) = totPsoil(npt) *fracPocc(casamet%isorder(npt))
      ELSE
        ! compute the steady-state litter and soil P pools
        casapool%plitter(npt,metb) = (avg_pleaf2met(npt)+avg_proot2met(npt))/casaflux%klitter(npt,metb)
        casapool%plitter(npt,str) = (avg_pleaf2str(npt)+avg_proot2str(npt))/casaflux%klitter(npt,str)
        casapool%plitter(npt,cwd) = (avg_pwood2cwd(npt))/casaflux%klitter(npt,cwd)

        casapool%psoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * casapool%ratioPCsoil(npt,mic)/casaflux%ksoil(npt,mic)
        casapool%psoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * casapool%ratioPCsoil(npt,slow)/casaflux%ksoil(npt,slow)
        casapool%psoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   * casapool%ratioPCsoil(npt,pass)/casaflux%ksoil(npt,pass)
        ! assign the mineral pools
        casapool%psoillab(npt)      = avg_psoillab(npt)
        casapool%psoilsorb(npt)     = avg_Psoilsorb(npt)
        casapool%psoilocc(npt)      = avg_Psoilocc(npt)
      ENDIF
    ENDIF
    ENDDO

    DEALLOCATE(totPsoil)

  END SUBROUTINE analyticPool


  SUBROUTINE open_out_casa_dump( ncfile )
    USE netcdf
    USE cable_diag_module, ONLY : def_dims, def_vars, def_var_atts, &
                                  put_var_nc, stderr_nc
    USE io_variables, only : patch

    CHARACTER(LEN=99), INTENT(IN) :: ncfile

    INTEGER, DIMENSION(num_vars)  :: varID      ! (1) lat (2) lon (3) tair, etc

    INTEGER, DIMENSION(num_dims)  :: dim_len    ! (1) tile, (2) soil, (3) time

    dim_len(1) = mp
    dim_len(2) = ms
    dim_len(3) = endday

    ncok = NF90_create(path = ncfile, cmode = NF90_noclobber, ncid = ncid_spin)
    IF (ncok /= NF90_noerr) CALL stderr_nc('ncdf creating ', ncfile)

    ! define dimensions: from name and length
    CALL def_dims(num_dims, ncid_spin, dimID, dim_len, dim_name )

    ! define variables: from name, type, dims
    CALL def_vars(num_vars, ncid_spin,  NF90_float, dimID, var_name, varID )

    ! define variable attributes
    CALL def_var_atts(ncfile, ncid_spin, varID )

    ncok = NF90_enddef(ncid_spin)

    CALL put_var_nc(ncid_spin, var_name(1), patch(:)%latitude )
    CALL put_var_nc(ncid_spin, var_name(2), patch(:)%longitude )

  END SUBROUTINE open_out_casa_dump


  SUBROUTINE write_casa_dump( casamet, casaflux, iday )
    USE netcdf
    USE cable_diag_module, ONLY : put_var_nc, stderr_nc

    TYPE (casa_flux), INTENT(INOUT) :: casaflux
    TYPE (casa_met),  INTENT(INOUT) :: casamet
    INTEGER, INTENT(IN) :: iday

    CALL put_var_nc(ncid_spin, var_name(3), casamet%tairk,     iday )
    CALL put_var_nc(ncid_spin, var_name(4), casamet%tsoil,     iday )
    CALL put_var_nc(ncid_spin, var_name(5), casamet%moist,     iday )
    CALL put_var_nc(ncid_spin, var_name(6), casaflux%cgpp,     iday )
    CALL put_var_nc(ncid_spin, var_name(7), casaflux%crmplant, iday )

  END SUBROUTINE write_casa_dump


  SUBROUTINE close_dump( ncfile )
    ! used for both closing the read file and the write file
    USE netcdf
    USE cable_diag_module, ONLY : stderr_nc

    CHARACTER(LEN=99), INTENT(IN) :: ncfile

    ncok = NF90_close(ncid_spin)
    IF (ncok /= NF90_noerr) CALL stderr_nc('Error closing file ', ncfile)

  END SUBROUTINE close_dump


  SUBROUTINE biogeochem(ktau,dels,idoy,veg,soil,casabiome,casapool,casaflux, &
                    casamet,casabal,phen)
    USE cable_def_types_mod
    USE casadimension
    USE casa_cnp_module
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: ktau
    REAL,    INTENT(IN)    :: dels
    INTEGER, INTENT(IN)    :: idoy
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
    INTEGER  npt,j

    xKNlimiting = 1.0
    call phenology(idoy,veg,phen)
    call avgsoil(veg,soil,casamet)
    call casa_rplant(veg,casabiome,casapool,casaflux,casamet)
    call casa_allocation(veg,soil,casabiome,casaflux,casamet,phen)
    call casa_xrateplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome, &
                         casamet,phen)
    call casa_coeffplant(xkleafcold,xkleafdry,xkleaf,veg,casabiome,casapool, &
                         casaflux,casamet)
    call casa_xnp(xnplimit,xNPuptake,veg,casabiome,casapool,casaflux,casamet)
    call casa_xratesoil(xklitter,xksoil,veg,soil,casamet)
    call casa_coeffsoil(xklitter,xksoil,veg,soil,casabiome,casaflux,casamet)

    IF (icycle>1) THEN
      call casa_xkN(xkNlimiting,casapool,casaflux,casamet,veg)
      DO j=1,mlitter
        casaflux%klitter(:,j) = casaflux%klitter(:,j)* xkNlimiting(:)
      ENDDO
      call casa_nuptake(veg,xkNlimiting,casabiome,casapool,casaflux,casamet)
      IF (icycle >2) call casa_puptake(veg,xkNlimiting,casabiome, &
                                       casapool,casaflux,casamet)
    ENDIF

    call casa_delplant(veg,casabiome,casapool,casaflux,casamet)
    call casa_delsoil(veg,casapool,casaflux,casamet)
    call casa_cnpcycle(veg,casabiome,casapool,casaflux,casamet)

    IF (icycle==1) call casa_ndummy(casapool)

    call casa_cnpbal(casapool,casaflux,casabal)
    call casa_cnpflux(casaflux,casabal)

  END SUBROUTINE biogeochem


END MODULE casa_spin_mod

