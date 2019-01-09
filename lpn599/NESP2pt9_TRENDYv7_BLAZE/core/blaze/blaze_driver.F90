SUBROUTINE BLAZE_DRIVER ( BLAZE, SF, met, casapool, casaflux, shootfrac, idoy, curyear, CTLFLAG )
!CLNSUBROUTINE BLAZE_DRIVER ( casapool, casaflux, lat, lon, shootfrac, ddvp09, ddvp15, ddprec, &
!CLN     ddTmin, ddTmax, ddwind,AvgAnnMaxFAPAR, modis_igbp, AvgAnnRainf, idoy, curyear, FLI, DFLI, FFDI, AB, &
!CLN     POPFLAG, CTLFLAG, BLAZEFLX, POP_TO, POP_CWD,POP_STR, IAC, popd, mnest, BLAZE_FSTEP &
!CLN     , AGL_wo1,AGL_wo2,AGL_wo3 )

  USE CABLE_COMMON_MODULE, ONLY: IS_LEAPYEAR, DOYSOD2MDHMS
  USE CASAVARIABLE,        ONLY: casa_pool, casa_flux
  USE BLAZE,               ONLY: RUN_BLAZE, TYPE_TURNOVER, BLAZE_TURNOVER, NTO, &
       METB, STR, CWD, LEAF, WOOD, FROOT, TYPE_BLAZE
  USE SIMFIRE_MOD,         ONLY: TYPE_SIMFIRE
  USE UTILS,               ONLY: Esatf
  USE BGCMODULE,           ONLY: tile_index

  IMPLICIT NONE

  ! CTLFLAG  : Only, when POP is on. Control whether only FLI (or FLI AND POP-related TO) 
  ! POPFLAG  : Check, how POP is set
  ! MODE     : 

  TYPE (casa_pool), INTENT(INOUT)      :: casapool
  REAL,DIMENSION(NCELLS,3),INTENT(OUT) :: IAC
  TYPE (casa_flux), INTENT(IN)         :: casaflux
  INTEGER,          INTENT(IN)         :: NCELLS,idoy, CurYear, POPFLAG, CTLFLAG
  REAL,DIMENSION(NCELLS),INTENT(IN)    :: ddvp09, ddvp15, ddprec, shootfrac
  REAL,DIMENSION(NCELLS),INTENT(IN)    :: ddTmin, ddTmax, ddwind, lon, lat
  REAL,DIMENSION(NCELLS),INTENT(IN)    :: POP_TO, POP_CWD, POP_STR
  REAL,DIMENSION(NCELLS,13), INTENT(INOUT) :: BLAZEFLX

  INTEGER, PARAMETER        :: NPOOLS = 3
!CLN  REAL,DIMENSION(NCELLS,3),OPTIONAL :: IAC
!CLN  TYPE (POP_TYPE), OPTIONAL         :: pop

  TYPE(TYPE_TURNOVER)   ,ALLOCATABLE,SAVE :: TO(:,:)
  REAL,   DIMENSION(:,:),ALLOCATABLE,SAVE :: AGL_w, AGL_g      ! Above-Ground Carbon
  REAL,   DIMENSION(:,:),ALLOCATABLE,SAVE :: BGL_w, BGL_g      ! Below-Ground Carbon
  REAL,   DIMENSION(NCELLS) :: AGL_wo1,AGL_wo2,AGL_wo3

  INTEGER,DIMENSION(NCELLS) :: modis_igbp
  INTEGER,DIMENSION(NCELLS) :: t1, t2
  REAL,   DIMENSION(NCELLS) :: AB, relhum, U10, FLI, DFLI, FFDI, popd, mnest
  REAL,   DIMENSION(NCELLS) :: AvgAnnMaxFAPAR, AvgAnnRainf, ag_lit, tot_lit

  INTEGER       :: MM, DD, i, np
  REAL          :: TSTP, C_CHKSUM
  REAL          :: CPLANT_g (ncells,3),CPLANT_w (ncells,3)
  REAL          :: CLITTER_g(ncells,3),CLITTER_w(ncells,3)
  LOGICAL       :: EOY
  LOGICAL, SAVE :: CALL1 = .TRUE.
  LOGICAL, SAVE :: YEAR1 = .TRUE.

  TYPE (TYPE_BLAZE)   :: BLAZE
  TYPE (TYPE_SIMFIRE) :: SF

  REAL    :: C_BIOMASS, C_FIRE
  LOGICAL,PARAMETER :: CLOSURE_TEST = .FALSE.
  CHARACTER         :: BLAZE_FSTEP*7

  ! INITIALISATION ============================================================
 
  IF ( BURNMODE .EQ. 0 ) RETURN

  BLAZEFLX = 0.

  !CLN ???
  t1 = tile_index(:,1)
  t2 = tile_index(:,2)

  CPLANT_g  = REAL(casapool%cplant (t1,:))
  CPLANT_w  = REAL(casapool%cplant (t2,:))
  CLITTER_g = REAL(casapool%clitter(t1,:))
  CLITTER_w = REAL(casapool%clitter(t2,:))

  ! CLN needs to be altered for beginning of spinup only!!!
  
  IF ( CALL1 ) THEN
     ALLOCATE( AGL_g(NCELLS,NPOOLS),AGL_w(NCELLS,NPOOLS) )
     ALLOCATE( BGL_g(NCELLS,NPOOLS),BGL_w(NCELLS,NPOOLS) )
     ! Initialise above / below-ground partitioning by using fluxes
     ! Grass
     ag_lit  = casaflux%fromPtoL(t1,METB,LEAF) * casaflux%kplant(t1,LEAF) * &
          MAX(casapool%cplant(t1,LEAF),1.e-5) 
     tot_lit = ag_lit + casaflux%fromPtoL(t1,METB,FROOT) * casaflux%kplant(t1,FROOT) * &
          MAX(casapool%cplant(t1,FROOT),1.e-5) 
     AGL_g(:,METB) = CLITTER_g(:,METB) * ag_lit / tot_lit

     ag_lit  = casaflux%fromPtoL(t1, STR,LEAF) * casaflux%kplant(t1,LEAF) * &
          MAX(casapool%cplant(t1,LEAF),1.e-5) 
     tot_lit = ag_lit + casaflux%fromPtoL(t1, STR,FROOT) * casaflux%kplant(t1,FROOT) * &
          MAX(casapool%cplant(t1,FROOT),1.e-5) 
     AGL_g(:,STR ) = CLITTER_g(:, STR) * ag_lit / tot_lit

     AGL_g(:,CWD)  = 0.
     ! Wood
     ag_lit  = casaflux%fromPtoL(t2,METB,LEAF) * casaflux%kplant(t2,LEAF) * &
          MAX(casapool%cplant(t2,LEAF),1.e-5) 
     tot_lit = ag_lit + casaflux%fromPtoL(t2,METB,FROOT) * casaflux%kplant(t2,FROOT) * &
          MAX(casapool%cplant(t2,FROOT),1.e-5) 
     AGL_w(:,METB) = CLITTER_w(:,METB) * ag_lit / tot_lit


     ag_lit  = casaflux%fromPtoL(t2, STR,LEAF) * casaflux%kplant(t2,LEAF) * &
          MAX(casapool%cplant(t2,LEAF),1.e-5) 
     tot_lit = ag_lit + casaflux%fromPtoL(t2, STR,FROOT) * casaflux%kplant(t2,FROOT) * &
          MAX(casapool%cplant(t2,FROOT),1.e-5) 
     AGL_w(:, STR) = CLITTER_w(:, STR) * ag_lit / tot_lit

     AGL_w(:, CWD) = CLITTER_w(:, CWD) * shootfrac

     CALL1 = .FALSE.

  END IF

  CALL DOYSOD2MDHMS( CurYear, idoy, 0, MM=MM, DD=DD )

  IF ( idoy .EQ. 366 .OR. ( idoy .EQ. 365 .AND. .NOT. is_leapyear(CurYear)) ) THEN
     EOY = .TRUE.
  ELSE
     EOY = .FALSE.
  END IF

  ! POPFLAG  
  IF ( POPFLAG .EQ. -1 ) THEN ! POP on -> Run on annual Timestep
     TSTP = 1.
  ELSE IF ( IS_LEAPYEAR(CurYear) ) THEN ! POP off
     TSTP = 1./366.
  ELSE
     TSTP = 1./365.
  END IF

  ! MAIN  ============================================================

!  IF ( CLOSURE_TEST ) &
!       C_CHKSUM = SUM(cplant_g) + SUM(cplant_w) + SUM(clitter_g) + SUM(clitter_w) &
!       + SUM(POP_TO) + SUM(POP_CWD) + SUM(POP_STR) 
!     
 
  ! Update above-ground-partition of metb/str litter pools
  AGL_g(:,METB) = (1. - casaflux%klitter(t1,METB)) * AGL_g(:,METB) + &
       casaflux%fromPtoL(t1,METB,LEAF) * casaflux%kplant(t1,LEAF) * CPLANT_g(:,LEAF)
  AGL_g(:, STR) = (1. - casaflux%klitter(t1,STR )) * AGL_g(:, STR) + &
       casaflux%fromPtoL(t1,STR ,LEAF) * casaflux%kplant(t1,LEAF) * CPLANT_g(:,LEAF) 
  
  AGL_w(:,METB) = (1. - casaflux%klitter(t2,METB)) * AGL_w(:,METB) + &
       casaflux%fromPtoL(t2,METB,LEAF) * casaflux%kplant(t2,LEAF) * CPLANT_w(:,LEAF)
  AGL_w(:, STR) = (1. - casaflux%klitter(t2,STR )) * AGL_w(:, STR) + &
       casaflux%fromPtoL(t2,STR ,LEAF) * casaflux%kplant(t2,LEAF) * CPLANT_w(:,LEAF) 
  AGL_w(:, CWD) = (1. - casaflux%klitter(t2,CWD )) * AGL_w(:, CWD) + &
       casaflux%fromPtoL(t2,CWD ,WOOD) * casaflux%kplant(t2,WOOD) * CPLANT_w(:,WOOD) * shootfrac

  ! If the pools are going to be updated split CLITTER into AGL and BGL 
  ! later add updated AGL at the end of this routine again
  IF ( BURNMODE .EQ. 1 .OR. CTLFLAG .EQ. -1 ) THEN
     DO i = 1, 3
        BGL_g(:,i) = CLITTER_g(:,i) - AGL_g(:,i)
        BGL_w(:,i) = CLITTER_w(:,i) - AGL_w(:,i)
     END DO
  ENDIF

  ! BLAZE used to compute ALL or FLI_ONLY
  IF ( BURNMODE .EQ. 1 .OR.  CTLFLAG .EQ. 1 ) THEN

     ! BURNMode 1: BLAZE computes mortality| BURNMode 2: BLAZE computes FLI only 
     !
     ! Ignition 0: GFED derived BA | Ignition 1: SIMFIRE Simulated BA
     ! Ignition is set in BLAZE at the moment.
     ! 
     ! CTLFLAG  -1: Compute Turnovers (using POP-TO) | 1: Compute FLI 
     
     ! Wind T. McVicar 201?
     ! Conversion to Windmax following S. Matthews, 2014 (pers. comm. so far)
     ! in km/h

     U10 = ddwind * 3.6 ! m/s -> km/h
     U10 = ( 214.7 * ( U10 + 10. ) ** (-1.6968 ) + 1. ) * U10
     
     DO i=1, NCELLS
        relhum(i) = 0.5*(MIN(1.,ddvp09(i) / Esatf(ddTmin(i))) + &
             MIN(1.,ddvp15(i) / Esatf(ddTmax(i)))) * 100.        ! [%]
     END DO

     CALL RUN_BLAZE( NCELLS, lat, lon, shootfrac,CPLANT_g, CPLANT_w, AGL_g, AGL_w, &
          BGL_g, BGL_w, ddprec, ddTMIN, ddTMAX, relhum, U10,AvgAnnMaxFAPAR, &
          modis_igbp, AvgAnnRainf, AB, FLI, DFLI, FFDI, TO, tstp, CurYear, idoy, 1, &
          POPFLAG, popd, mnest,BLAZE_FSTEP )
     IF ( idoy .EQ. 1 ) IAC = 0.
     
!!CRM     IF ( POPFLAG .NE. 0 ) THEN       
!!CRM        WHERE(AB .GT. 0. )
!!CRM           IAC (:,1) = (FLI * AB + IAC(:,1) * IAC(:,2)) / (IAC(:,3)+1.)
!!CRM           IAC (:,2) = IAC(:,2)+AB
!!CRM           IAC (:,3) = IAC(:,3)+1.
!!CRM        END WHERE
!!CRM     END IF

  ! compute c-pool turnovers after POP has provided biomass TO 
  ELSE IF ( CTLFLAG .EQ. -1 ) THEN
!     IF ( .NOT. PRESENT(POP_TO) ) STOP "Provide POP_TO to blaze_casa.f90!"
     DO np = 1, NCELLS
        CALL BLAZE_TURNOVER( AB(np), CPLANT_g(np,:), CPLANT_w(np,:), AGL_g(np,:), &
             AGL_w(np,:), BGL_g(np,:), BGL_w(np,:),shootfrac(np),TO(np,:), &
             POPFLAG, BLAZEFLX(np,:), POP_TO(np) )
     END DO
  ELSE
     STOP "Wrong MODE in blaze_driver.f90!"
  ENDIF

  ! When TURNOVERS have been computed, update LITTER with AGL
  ! then add POP non-fire disturbance litter to AGL
  IF ( BURNMODE .EQ. 1 .OR. CTLFLAG .EQ. -1 ) THEN
     DO i = 1, 3
        CLITTER_g(:,i) = BGL_g(:,i) + AGL_g(:,i)
        CLITTER_w(:,i) = BGL_w(:,i) + AGL_w(:,i)
     END DO
     ! Update AGL only (BGL not saved)
     IF ( POPFLAG .NE. 0 ) THEN
        AGL_w(:, CWD) = AGL_w(:, CWD) + POP_CWD * shootfrac
        AGL_w(:, STR) = AGL_w(:, STR) + POP_STR 
     ENDIF
  ENDIF


  AGL_wo1 = ( AGL_w(:,METB))
  AGL_wo2 = ( AGL_w(:,STR ))
  AGL_wo3 = ( AGL_w(:,CWD ))


  ! C - CLOSURE check
  IF ( CLOSURE_TEST ) THEN
     C_BIOMASS = (SUM( casapool%cplant ) + SUM(casapool%clitter) ) - &
          (SUM(cplant_g) + SUM(cplant_w) + SUM(clitter_g) + SUM(clitter_w))
          
     C_FIRE = SUM(BLAZEFLX(:,1:2)) + SUM(BLAZEFLX(:,6:9)) + SUM(BLAZEFLX(:,11:13))
     IF ( ABS(C_BIOMASS - C_FIRE) .GT. .5 ) THEN
        PRINT*," IMBALANCE in BLAZE ", C_BIOMASS , C_FIRE, C_BIOMASS - C_FIRE
        PRINT*," GRASS : "
        PRINT*," Live 1 : ", SUM(cplant_g(:,1)),SUM(casapool%cplant(t1,1)), -SUM(BLAZEFLX(:,11)),&
             SUM(cplant_g(:,1))+SUM(BLAZEFLX(:,11))-SUM(casapool%cplant(t1,1))
        PRINT*," Live 2 : ", SUM(cplant_g(:,2)),SUM(casapool%cplant(t1,2)),SUM(cplant_g(:,2))-SUM(casapool%cplant(t1,2))
        PRINT*," Live 3 : ", SUM(cplant_g(:,3)),SUM(casapool%cplant(t1,3)),SUM(cplant_g(:,3))-SUM(casapool%cplant(t1,3))

        PRINT*," Litter1: ", SUM(clitter_g(:,1)),SUM(casapool%clitter(t1,1)),-SUM(BLAZEFLX(:,12)),&
             (SUM(clitter_g(:,1))-SUM(casapool%clitter(t1,1)))+SUM(BLAZEFLX(:,12))
        PRINT*," Litter2: ", SUM(clitter_g(:,2)),SUM(casapool%clitter(t1,2)),-SUM(BLAZEFLX(:,13)),&
             (SUM(clitter_g(:,2))-SUM(casapool%clitter(t1,2)))+SUM(BLAZEFLX(:,13))
        PRINT*," Litter3: ", SUM(clitter_g(:,3)),SUM(casapool%clitter(t1,3)),(SUM(clitter_g(:,3))-SUM(casapool%clitter(t1,3)))
        PRINT*," WOOD : "
        PRINT*," Live  LEAF:  ", SUM(cplant_w(:,LEAF)),SUM(casapool%cplant(t2,LEAF)),-(SUM(BLAZEFLX(:,1))+SUM(BLAZEFLX(:,3))),&
             SUM(cplant_w(:,LEAF))-SUM(casapool%cplant(t2,LEAF))+(SUM(BLAZEFLX(:,1))+SUM(BLAZEFLX(:,3))) 
        PRINT*," Live  WOOD:  ", SUM(cplant_w(:,WOOD)),SUM(casapool%cplant(t2,WOOD)),-(SUM(BLAZEFLX(:,2))+SUM(BLAZEFLX(:,4:5))),&
             SUM(cplant_w(:,WOOD))-SUM(casapool%cplant(t2,WOOD))+(SUM(BLAZEFLX(:,2))+SUM(BLAZEFLX(:,4:5)))  
        
        PRINT*," Live  FROOT: ", SUM(cplant_w(:,FROOT)),SUM(casapool%cplant(t2,FROOT)),-SUM(BLAZEFLX(:,9:10)), &
             SUM(cplant_w(:,FROOT))-SUM(casapool%cplant(t2,FROOT))+SUM(BLAZEFLX(:,9:10)) 
        PRINT*," Litter METB: ", SUM(clitter_w(:,METB)),SUM(casapool%clitter(t2,METB)),-SUM(BLAZEFLX(:,6)),&
             SUM(clitter_w(:,METB))-SUM(casapool%clitter(t2,METB))+SUM(BLAZEFLX(:,6))
        PRINT*," Litter STR:  ", SUM(clitter_w(:,STR )),SUM(casapool%clitter(t2,STR )),-SUM(BLAZEFLX(:,7))+SUM(BLAZEFLX(:,3:4))+SUM(BLAZEFLX(:,9)),&
             SUM(clitter_w(:,STR ))-SUM(casapool%clitter(t2,STR ))+SUM(BLAZEFLX(:,7))-(SUM(BLAZEFLX(:,3:4))+SUM(BLAZEFLX(:,9)))
        PRINT*," Litter CWD : ", SUM(clitter_w(:,CWD)),SUM(casapool%clitter(t2,CWD)),-SUM(BLAZEFLX(:,8))+SUM(BLAZEFLX(:,5)),&
             SUM(clitter_w(:,CWD))-SUM(casapool%clitter(t2,CWD))+SUM(BLAZEFLX(:,8))-SUM(BLAZEFLX(:,5))

        PRINT*,"MTO(LEAF )%TO_STR : ", SUM(BLAZEFLX(:,3))
        PRINT*,"MTO(WOOD )%TO_STR : ", SUM(BLAZEFLX(:,4))
        PRINT*,"MTO(WOOD )%TO_CWD : ", SUM(BLAZEFLX(:,5))
        PRINT*,"MTO(FROOT)%TO_STR : ", SUM(BLAZEFLX(:,9))
        PRINT*,"MTO(MLIT )%TO_ATM : ", -SUM(BLAZEFLX(:,6))
        PRINT*,"MTO(SLIT )%TO_ATM : ", -SUM(BLAZEFLX(:,7))
        PRINT*,"MTO(CLIT )%TO_ATM : ", -SUM(BLAZEFLX(:,8))
        PRINT*,"AB * AGL_g(METB)  : ", -SUM(BLAZEFLX(:,12))
        PRINT*,"AB * AGL_g(STR )  : ", -SUM(BLAZEFLX(:,13))
     ENDIF
  ENDIF

  casapool%cplant (t1,:) = DBLE(CPLANT_g )
  casapool%cplant (t2,:) = DBLE(CPLANT_w )
  casapool%clitter(t1,:) = DBLE(CLITTER_g)
  casapool%clitter(t2,:) = DBLE(CLITTER_w)

END SUBROUTINE BLAZE_DRIVER
