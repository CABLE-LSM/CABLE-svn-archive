SUBROUTINE casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux,casamet, &
                        casabal,phen)
  USE cable_def_types_mod
  USE casadimension
  USE casaparm
  USE casavariable
  USE phenvariable
  USE cable_common_module, only: cable_user
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
  if(casamet%iveg2(npt)>0 ) then
    IF (icycle<2) THEN
      casapool%Nplant(npt,:) = casapool%ratioNCplant(npt,:)  &
                             * casapool%cplant(npt,:)
      casapool%Nlitter(npt,:)= casapool%ratioNClitter(npt,:) &
                             * casapool%clitter(npt,:)
      casapool%Nsoil(npt,:)  = casapool%ratioNCsoil(npt,:)   &
                             * casapool%Csoil(npt,:)
      casapool%nsoilmin(npt) = 2.0
      casabal%sumnbal(npt)   = 0.0
      if(casamet%iveg2(npt)==grass) then
         casapool%nplant(npt,wood) = 0.0
         casapool%nlitter(npt,cwd) = 0.0
      endif
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
      if(casamet%iveg2(npt)==grass) then
         casapool%pplant(npt,wood) = 0.0
         casapool%plitter(npt,cwd) = 0.0
      endif
    ENDIF
  else
     casapool%cplant(npt,:)=0.0; casapool%clitter(npt,:)=0.0; casapool%csoil(npt,:) = 0.0; casapool%clabile(npt) = 0.0
     casapool%nplant(npt,:)=0.0; casapool%nlitter(npt,:)=0.0; casapool%nsoil(npt,:) = 0.0; casapool%nsoilmin(npt) = 0.0
     casapool%pplant(npt,:)=0.0; casapool%plitter(npt,:)=0.0; casapool%psoil(npt,:) = 0.0
     casapool%psoillab(npt) = 0.0; casapool%psoilsorb(npt) = 0.0; casapool%psoilocc(npt) = 0.0
     casabal%sumcbal(npt) =0.0; casabal%sumnbal(npt) =0.0; casabal%sumpbal(npt) = 0.0
  endif

!! vh_js  !! 
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

92    format(5(i6,',',2x),5(f18.6,',',2x),2(i6,',',2x),100(f18.6,',',2x))
END SUBROUTINE casa_poolout


