MODULE analyticpool_casa_mod

IMPLICIT NONE

#define  UM_BUILD YES

#ifndef UM_BUILD

CONTAINS

SUBROUTINE analyticpool(kend,veg,soil,casabiome,casapool,                                 &
                          casaflux,casamet,casabal,phen,                                    &
                          avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd, &
                          avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd, &
                          avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd, &
                          avgcgpp, avgcnpp, avgnuptake, avgpuptake,                         &
                          avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil,                 &
                          avgratioNCsoilmic,avgratioNCsoilslow,avgratioNCsoilpass,          &
                          avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc,                 &
                          avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr,                   &
                          avg_annual_cnpp)
  USE cable_def_types_mod
  USE cable_carbon_module
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: kend
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_biome),            INTENT(INOUT) :: casabiome
  TYPE (casa_pool),             INTENT(INOUT) :: casapool
  TYPE (casa_flux),             INTENT(INOUT) :: casaflux
  TYPE (casa_met),              INTENT(INOUT) :: casamet
  TYPE (casa_balance),          INTENT(INOUT) :: casabal
  TYPE (phen_variable),         INTENT(INOUT) :: phen


  ! local variables
  real, dimension(mp)      :: avgcleaf2met,avgcleaf2str,avgcroot2met,avgcroot2str,avgcwood2cwd
  real, dimension(mp)      :: avgnleaf2met,avgnleaf2str,avgnroot2met,avgnroot2str,avgnwood2cwd
  real, dimension(mp)      :: avgpleaf2met,avgpleaf2str,avgproot2met,avgproot2str,avgpwood2cwd
  real, dimension(mp)      :: avgcgpp, avgcnpp, avgnuptake, avgpuptake
  real(r_2), dimension(mp) :: avgxnplimit,avgxkNlimiting,avgxklitter,avgxksoil
  real, dimension(mp)      :: avgratioNCsoilmic,  avgratioNCsoilslow,  avgratioNCsoilpass
  real,      dimension(mp) :: avgnsoilmin,avgpsoillab,avgpsoilsorb,avgpsoilocc
  real,      dimension(mp) :: avg_af, avg_aw, avg_ar, avg_lf, avg_lw, avg_lr
  real,      dimension(mp) :: avg_annual_cnpp

  ! local variables
  REAL(r_2), DIMENSION(mso) :: Psorder,pweasoil,xpsoil50
  REAL(r_2), DIMENSION(mso) :: fracPlab,fracPsorb,fracPocc,fracPorg
  REAL(r_2), DIMENSION(mp)  :: totpsoil
  INTEGER  npt,nout,nso

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
  DATA psorder/61.3,103.9,92.8,136.9,98.2,107.6,84.1,110.1,35.4,41.0,51.5,190.6/
  DATA pweasoil/0.05,0.04,0.03,0.02,0.01,0.009,0.008,0.007,0.006,0.005,0.004,0.003/
  DATA fracpLab/0.08,0.08,0.10,0.02,0.08,0.08,0.08,0.06,0.02,0.05,0.09,0.05/
  DATA fracPsorb/0.32,0.37,0.57,0.67,0.37,0.37,0.37,0.32,0.24,0.22,0.21,0.38/
  DATA fracPocc/0.36,0.38,0.25,0.26,0.38,0.38,0.38,0.44,0.38,0.38,0.37,0.45/
  DATA fracPorg/0.25,0.17,0.08,0.05,0.17,0.17,0.17,0.18,0.36,0.35,0.34,0.12/
  DATA xpsoil50/7.6,4.1,4.2,3.4,4.1,4.1,4.8,4.1,6.9,6.9,6.9,1.7/

  INTEGER :: nyear,iyear
  real year

  ! compute the mean litter input in g(C, N and P)/day from plant pools
  casaflux%fromLtoS = 0.0
  casaflux%fromStoS = 0.0

  casabal%sumcbal(:)   = 0.0
  casabal%sumnbal(:)   = 0.0
  casabal%sumpbal(:)   = 0.0

  do npt=1,mp
  if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) THEN
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
!    casaflux%fromStoS(npt,mic,slow)  = 0.45 * (0.997 - 0.009 *soil%clay(npt))
!    casaflux%fromStoS(npt,mic,pass)  = 0.45

    casaflux%fromStoS(npt,slow,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.997 - 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,mic)  = (0.85 - 0.68 * (soil%clay(npt)+soil%silt(npt))) &
                                     * (0.003 + 0.032*soil%clay(npt))
    casaflux%fromStoS(npt,pass,slow) = 0.45 * (0.003 + 0.009 * soil%clay(npt) )

    casaflux%klitter(npt,metb) = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),metb)
    casaflux%klitter(npt,str)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),str)&
                               * exp(-3.0*casabiome%fracLigninplant(veg%iveg(npt),leaf))
    casaflux%klitter(npt,cwd)  = avgxkNlimiting(npt) * avgxklitter(npt)*casabiome%litterrate(veg%iveg(npt),cwd)

    casaflux%ksoil(npt,mic)    = avgxksoil(npt)*casabiome%soilrate(veg%iveg(npt),mic)  &
                               * (1.0 - 0.75 *(soil%silt(npt)+soil%clay(npt)))
    casaflux%ksoil(npt,slow)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),slow)
    casaflux%ksoil(npt,pass)   = avgxksoil(npt) * casabiome%soilrate(veg%iveg(npt),pass)

    if(veg%iveg(npt)==cropland) THEN     ! for cultivated land type
       casaflux%ksoil(npt,mic)  = casaflux%ksoil(npt,mic) * 1.25
       casaflux%ksoil(npt,slow) = casaflux%ksoil(npt,slow)* 1.5
       casaflux%ksoil(npt,pass) = casaflux%ksoil(npt,pass)* 1.5
    endif
  endif
  enddo


  do npt=1,mp
   if(casamet%iveg2(npt)/=icewater.and.avgcnpp(npt) > 0.0) then

      ! Compute steady-state plant C pool sizes
      casapool%cplant(npt,leaf) = (avg_annual_cnpp(npt) * avg_af(npt)) / &
                                   avg_lf(npt)
      casapool%cplant(npt,wood) = (avg_annual_cnpp(npt) * avg_aw(npt)) / &
                                   avg_lw(npt)
      casapool%cplant(npt,froot) = (avg_annual_cnpp(npt) * avg_ar(npt)) / &
                                    avg_lr(npt)

    ! compute steady-state litter and soil C pool sizes
     casapool%clitter(npt,metb) = (avgcleaf2met(npt)+avgcroot2met(npt))/casaflux%klitter(npt,metb)
     casapool%clitter(npt,str) = (avgcleaf2str(npt)+avgcroot2str(npt))/casaflux%klitter(npt,str)
     casapool%clitter(npt,cwd) = (avgcwood2cwd(npt))/casaflux%klitter(npt,cwd)
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

      casabal%clitterlast = casapool%clitter
      casabal%csoillast   = casapool%csoil

      if(icycle <=1) then
         casapool%nlitter(npt,:)= casapool%rationclitter(npt,:) * casapool%clitter(npt,:)
         casapool%nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   * casapool%Csoil(npt,:)
         casapool%nsoilmin(npt) = 2.0
         casabal%sumnbal(npt)   = 0.0
       ELSE
       ! compute steady-state litter and soil N pool sizes
         casapool%nlitter(npt,metb) = (avgnleaf2met(npt)+avgnroot2met(npt))/casaflux%klitter(npt,metb)
         casapool%nlitter(npt,str) = (avgnleaf2str(npt)+avgnroot2str(npt))/casaflux%klitter(npt,str)
         casapool%nlitter(npt,cwd) = (avgnwood2cwd(npt))/casaflux%klitter(npt,cwd)

         casapool%nsoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * avgratioNCsoilmic(npt)/casaflux%ksoil(npt,mic)
         casapool%nsoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * avgratioNCsoilslow(npt)/casaflux%ksoil(npt,slow)
         casapool%nsoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   * avgratioNCsoilpass(npt)/casaflux%ksoil(npt,pass)
          casapool%Nsoilmin(npt)    = avgnsoilmin(npt)

        ENDIF

        IF (icycle<=2) THEN
            totpsoil(npt)          = psorder(casamet%isorder(npt)) *xpsoil50(casamet%isorder(npt))
           casapool%plitter(npt,:)= casapool%Nlitter(npt,:)/casapool%ratioNPlitter(npt,:)
            casapool%psoil(npt,:)  = casapool%Nsoil(npt,:)/casapool%ratioNPsoil(npt,:)
            ! why is this commented here but used in UM
            ! casapool%plitter(npt,:)= casapool%ratiopclitter(npt,:)  * casapool%clitter(npt,:)
            ! casapool%psoil(npt,:)  = casapool%ratioPCsoil(npt,:)    * casapool%Csoil(npt,:)
            casapool%psoillab(npt) = totpsoil(npt) *fracpLab(casamet%isorder(npt))
            casapool%psoilsorb(npt)= casaflux%psorbmax(npt) * casapool%psoillab(npt) &
                                    /(casaflux%kmlabp(npt)+casapool%psoillab(npt))
            casapool%psoilocc(npt) = totpsoil(npt) *fracPocc(casamet%isorder(npt))
        ELSE
        ! compute the steady-state litter and soil P pools
          casapool%plitter(npt,metb) = (avgpleaf2met(npt)+avgproot2met(npt))/casaflux%klitter(npt,metb)
          casapool%plitter(npt,str) = (avgpleaf2str(npt)+avgproot2str(npt))/casaflux%klitter(npt,str)
          casapool%plitter(npt,cwd) = (avgpwood2cwd(npt))/casaflux%klitter(npt,cwd)

          casapool%psoil(npt,mic)   = (casaflux%fromLtoS(npt,mic,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb)   &
                                     +casaflux%fromLtoS(npt,mic,str) *casaflux%klitter(npt,str)*casapool%clitter(npt,str)  &
                                     +casaflux%fromLtoS(npt,mic,cwd) *casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) ) &
                                   * (casapool%ratioNCsoil(npt,mic)/casapool%ratioNPsoil(npt,mic))/casaflux%ksoil(npt,mic)
          casapool%psoil(npt,slow)  = (casaflux%fromLtoS(npt,slow,metb)*casaflux%klitter(npt,metb)*casapool%clitter(npt,metb) &
                                     + casaflux%fromLtoS(npt,slow,str)*casaflux%klitter(npt,str)*casapool%clitter(npt,str) &
                                     + casaflux%fromLtoS(npt,slow,cwd)*casaflux%klitter(npt,cwd)*casapool%clitter(npt,cwd) &
                                     + casaflux%fromStoS(npt,slow,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)  ) &
                                   * (casapool%ratioNCsoil(npt,slow)/casapool%ratioNPsoil(npt,slow))/casaflux%ksoil(npt,slow)
          casapool%psoil(npt,pass)  = (casaflux%fromStoS(npt,pass,mic) *casaflux%ksoil(npt,mic) *casapool%csoil(npt,mic)    &
                                     +casaflux%fromStoS(npt,pass,slow)*casaflux%ksoil(npt,slow)*casapool%csoil(npt,slow) ) &
                                   *  (casapool%ratioNCsoil(npt,pass)/casapool%ratioNPsoil(npt,pass))/casaflux%ksoil(npt,pass)
          ! assign the mineral pools
          casapool%psoillab(npt)      = avgpsoillab(npt)
          casapool%psoilsorb(npt)     = avgPsoilsorb(npt)
          casapool%psoilocc(npt)      = avgPsoilocc(npt)
        ENDIF
  ENDIF
  ENDDO

  END SUBROUTINE analyticpool
#endif

END MODULE analyticpool_casa_mod
