MODULE POPLUC_CONSTANTS
  USE TYPEdef, ONLY: dp, i4b

  INTEGER(i4b),PARAMETER ::  LENGTH_SECDF_HISTORY = 4000
  INTEGER(i4b),PARAMETER :: AGE_MAX = 1000  
  INTEGER(i4b),PARAMETER :: disturbance_interval = 100 ! N.B. needs to be the same as veg%disturbacne_interval
  LOGICAL, PARAMETER :: IFHARVEST=.FALSE.
  INTEGER(i4b), PARAMETER :: ROTATION=70
  INTEGER(i4b), PARAMETER :: nLU=3 ! number of land-use tiles (pf, sf, grass)
  INTEGER(i4b), PARAMETER :: nTrans=4 ! number of possible gross transition types (ptog, ptos, stog, gtos)

END MODULE POPLUC_CONSTANTS

!*******************************************************************************
MODULE POPLUC_Types
  USE TYPEdef, ONLY: dp, i4b
  USE POPLUC_Constants, ONLY: LENGTH_SECDF_HISTORY, AGE_MAX

  TYPE POPLUC_TYPE
     INTEGER(i4b),POINTER :: it
     INTEGER(i4b),POINTER :: np
     INTEGER(i4b),POINTER :: firstyear
     INTEGER(i4b),POINTER :: thisyear
     INTEGER(i4b), DIMENSION(:),POINTER :: n_event ! number of secondary forest transitions
     REAL(dp), DIMENSION(:),POINTER :: latitude, longitude
     REAL(dp), DIMENSION(:),POINTER :: primf, secdf, grass,       &  ! land cover types
          ptos, ptog, stop, stog, gtop, gtos,    & ! transitions
          frac_primf, frac_forest
     REAL(dp), DIMENSION(:,:),POINTER ::  freq_age_primary, freq_age_secondary, &
          biomass_age_primary, biomass_age_secondary   
     REAL(dp), DIMENSION(:,:),POINTER :: age_history_secdf, area_history_secdf
     REAL(dp), DIMENSION(:,:),POINTER :: FNEP, Clitt, Csoil, Cbiomass
     REAL(dp), DIMENSION(:,:),POINTER :: FHarvest, FClearance, FTransferNet
     REAL(dp), DIMENSION(:,:),POINTER :: FTransferGross
     REAL(dp), DIMENSION(:),POINTER :: pharv, smharv, syharv
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProd, ClearProd
     REAL(dp), DIMENSION(:,:),POINTER :: fracHarvProd, fracClearProd
     REAL(dp), DIMENSION(:,:),POINTER :: HarvProdLoss, ClearProdLoss
     REAL(dp), DIMENSION(:),POINTER :: fracHarvResid, fracHarvSecResid, fracClearResid

  END TYPE POPLUC_TYPE

END MODULE POPLUC_Types
!*******************************************************************************

MODULE POPLUC_Module

  !-------------------------------------------------------------------------------
  ! * This module contains all subroutines for POPLUC calcs at a single time step.
  !-------------------------------------------------------------------------------
  USE TYPEdef, ONLY: sp, i4b
  USE POPLUC_Types
  USE POPLUC_Constants
  USE casavariable, ONLY: casa_pool, casa_balance
  USE POP_Types, ONLY: POP_TYPE
  USE cable_IO_vars_module, ONLY: landpt, patch
  USE CABLE_LUC_EXPT, ONLY: LUC_EXPT_TYPE

CONTAINS

  !*******************************************************************************
  SUBROUTINE ZeroPOPLUC(POPLUC)
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER:: g,np

    np = popluc%np

    POPLUC%firstyear = 0
    POPLUC%thisyear = 0
    POPLUC%primf = 0
    POPLUC%secdf = 0
    POPLUC%grass = 0
    POPLUC%ptos = 0
    POPLUC%ptog = 0
    POPLUC%stop = 0
    POPLUC%stog = 0
    POPLUC%gtop = 0
    POPLUC%gtos = 0
    POPLUC%frac_forest = 0
    POPLUC%frac_primf = 0
    POPLUC%area_history_secdf = 0
    POPLUC%age_history_secdf = 0
    POPLUC%n_event = 0
    POPLUC%freq_age_secondary = 0
    POPLUC%freq_age_primary = 0
    POPLUC%biomass_age_primary = 0
    POPLUC%biomass_age_secondary = 0
    POPLUC%FNEP = 0
    POPLUC%Clitt = 0
    POPLUC%Csoil = 0
    POPLUC%Cbiomass = 0
    POPLUC%FHarvest = 0
    POPLUC%FClearance = 0
    POPLUC%FTransferNet = 0
    POPLUC%FTransferGross = 0
    POPLUC%pharv = 0
    POPLUC%smharv = 0
    POPLUC%syharv = 0
    POPLUC%HarvProd = 0
    POPLUC%HarvProdLoss = 0
    POPLUC%ClearProd = 0
    POPLUC%ClearProdLoss = 0


  END SUBROUTINE ZeroPOPLUC
  !*******************************************************************************

  SUBROUTINE execute_luc_event(from_state,to_state,frac_change_grid,g,POPLUC) 
    ! Execute a transition between land use types (states)
    !  frac_change_grid = fractional change in unit fraction of grid cell this year

    IMPLICIT NONE

    CHARACTER(5), INTENT(IN) :: from_state, to_state
    REAL(dp), INTENT(INOUT):: frac_change_grid
    TYPE(POPLUC_TYPE),INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: g  ! grid cell index
    REAL :: frac_open_grid, remaining
    INTEGER(i4b) :: n  ! position of new element in POPLUC%SecFor array
    INTEGER(i4b) :: i
    frac_open_grid = 1.0 -POPLUC%frac_forest(g)
    n = POPLUC%n_event(g)

    IF (from_state=='PRIMF') THEN

       IF(frac_change_grid.GT.POPLUC%primf(g)) THEN
          !PRINT*, "Warning: requested reduction in primary forest area &
          !     exceeds primary forest area"
          IF (to_state=='SECDF') POPLUC%ptos(g) = POPLUC%primf(g)
          IF (to_state=='C3ANN') POPLUC%ptog(g) = POPLUC%primf(g)
          frac_change_grid = POPLUC%primf(g)
          POPLUC%primf(g) = 0.0
       ELSE
          POPLUC%primf(g) = POPLUC%primf(g) &
               -frac_change_grid
       ENDIF

       IF (to_state=='SECDF') THEN 
          ! Transition from primary -> secondary forest(ptos)
          POPLUC%n_event(g) = POPLUC%n_event(g)+1
          n = POPLUC%n_event(g)
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid

       ENDIF

    ELSEIF (from_state=='SECDF') THEN

       IF (to_state=='PRIMF') THEN
          print*, "Error: cannot create primary forest from secondary forest"
          STOP

       ELSE
          ! Transition from secondary -> non forest (stog)
          ! Assumption: youngest stands cleared first
!!$          remaining = frac_change_grid
!!$          i = n
!!$
!!$          DO WHILE (remaining > 0.0 .and. i >0 )
!!$
!!$             IF (POPLUC%area_history_secdf(g,i).GE.remaining) THEN
!!$                POPLUC%area_history_secdf(g,i) =POPLUC%area_history_secdf(g,i)&
!!$                     - remaining
!!$                remaining = 0.0
!!$             ELSE
!!$                remaining = remaining - POPLUC%area_history_secdf(g,i)
!!$                i = i-1
!!$                POPLUC%n_event(g) = POPLUC%n_event(g)-1
!!$             ENDIF
!!$
!!$          ENDDO
   
          remaining = frac_change_grid
          i = 1
          DO WHILE (remaining > 0.0 .and. i <= age_max )
             IF (POPLUC%freq_age_secondary(g,i).GE.remaining) THEN
                POPLUC%freq_age_secondary(g,i) =POPLUC%freq_age_secondary(g,i) &
                     - remaining
                remaining = 0.0
             ELSE
                remaining = remaining - POPLUC%freq_age_secondary(g,i)
                POPLUC%freq_age_secondary(g,i) = 0.0
                i = i+1
             ENDIF

          ENDDO
          !if (remaining.gt.frac_change_grid) POPLUC%stog(g) = POPLUC%stog(g)-remaining
          if (remaining.gt.0.0) POPLUC%stog(g) = POPLUC%stog(g)-remaining
       ENDIF

    ELSEIF (to_state=='SECDF') THEN

       POPLUC%n_event(g) = POPLUC%n_event(g)+1
       n = POPLUC%n_event(g) 
       ! Transition from non-forest to secondary forest (gtos)

       if (frac_change_grid.LE.frac_open_grid) THEN
          POPLUC%area_history_secdf(g,n) = frac_change_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_change_grid
       ELSE
          POPLUC%area_history_secdf(g,n) = frac_open_grid
          POPLUC%age_history_secdf(g,n) = 0
          POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + frac_open_grid
          ! gtos to frac_change_grid here!!
          POPLUC%gtos(g) =  frac_open_grid


       ENDIF

    ELSEIF (to_state=='PRIMF') THEN

       print*, "Error: cannot create primary forest from non-forest"
       STOP

    ENDIF


  ENDSUBROUTINE execute_luc_event


  !*******************************************************************************
  SUBROUTINE CALCULATE_WEIGHTS(POPLUC, g)
    ! Calculates weights (fraction of total forest area on grid cell)
    !for primary and secondary forest stands up to specified maximum stand age


    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: g
    REAL(dp), PARAMETER :: EPS= 1.e-9
    INTEGER(i4b) :: age , i, iage
    REAL(dp) :: fac, disturbance_freq

    ! First get relative weights for primary forest
    disturbance_freq=1.0/REAL(disturbance_interval,dp)

    !fac = POPLUC%frac_primf/POPLUC%frac_forest
    fac = 1.0
    DO iage = 1, age_max
       POPLUC%freq_age_primary(g,iage) =  REALExponential(disturbance_freq,REAL(iage-1,dp))
       POPLUC%freq_age_secondary(g,iage) = 0.0
    END DO
    POPLUC%freq_age_primary(g,:) =  POPLUC%freq_age_primary(g,:) &
         / sum(POPLUC%freq_age_primary(g,:))*fac


    !  Loop through secondary forest stands to transfer weights
    !fac = POPLUC%frac_forest
    fac = 1.0
    DO i = 1, POPLUC%n_event(g)

       age = POPLUC%age_history_secdf(g,i)
       POPLUC%freq_age_secondary(g,age+1) = POPLUC%area_history_secdf(g,i)/fac
       ! write(*,*) 'calc_wt', age,  POPLUC%SecFor%area(i), fac
    ENDDO
    !STOP

  END SUBROUTINE CALCULATE_WEIGHTS
  !*******************************************************************************
  SUBROUTINE INCREMENT_AGE(POPLUC,g)

    IMPLICIT NONE

    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b) :: n_event, i, j
    INTEGER(i4b), INTENT(IN) :: g
    REAL(dp):: area, remaining

    n_event =  POPLUC%n_event(g)

    POPLUC%freq_age_secondary(g,2:age_max)=POPLUC%freq_age_secondary(g,1:age_max-1)
    POPLUC%freq_age_secondary(g,1) = 0.0

    ! adjust secondary age distribution for secondary mature forest harvest area
    if (POPLUC%smharv(g).gt.0) then

       remaining = POPLUC%smharv(g)
       i = age_max
       do while (remaining.gt.0.0)
          if (POPLUC%freq_age_secondary(g,i) .gt. remaining) then

             POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i) - remaining
             POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + remaining
             remaining = 0.0
          else
             POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + POPLUC%freq_age_secondary(g,i)
             POPLUC%freq_age_secondary(g,i) = 0.0
             i = i-1;
          end if
          if (i.lt.2) then
             POPLUC%smharv(g) = POPLUC%smharv(g) - remaining
             remaining = 0.0
          endif
       enddo
    endif

    IF (IFHARVEST .and. POPLUC%freq_age_secondary(g,ROTATION+1).gt.0.0) THEN
       POPLUC%freq_age_secondary(g,1) = POPLUC%freq_age_secondary(g,1) + POPLUC%freq_age_secondary(g,ROTATION+1)
       POPLUC%freq_age_secondary(g,ROTATION+1) = 0.0
    ENDIF

    ! adjust secondary age distribution for natural disturbance
    i = age_max
    DO i = age_max, 2 , -1
       POPLUC%freq_age_secondary(g,1) =  POPLUC%freq_age_secondary(g,1) +  POPLUC%freq_age_secondary(g,i)/disturbance_interval
       POPLUC%freq_age_secondary(g,i) = POPLUC%freq_age_secondary(g,i)*(1. - 1./disturbance_interval)
    ENDDO


!!$    i = 1
!!$    DO WHILE (i.LE.n_event)
!!$       IF (POPLUC%area_history_secdf(g,i)>0) THEN
!!$          POPLUC%age_history_secdf(g,i) = POPLUC%age_history_secdf(g,i)+1
!!$          IF (IFHARVEST) THEN
!!$
!!$             IF ( POPLUC%age_history_secdf(g,i) == ROTATION) THEN
!!$                ! Harvest: remove this stand from record and create\
!!$                ! a new secondary forest stand of age 0
!!$                area = POPLUC%area_history_secdf(g,i)
!!$                DO j=i,n_event-1
!!$                   POPLUC%age_history_secdf(g,j) = POPLUC%age_history_secdf(g,j+1)
!!$                   POPLUC%area_history_secdf(g,j) = POPLUC%area_history_secdf(g,j+1)
!!$                ENDDO
!!$                POPLUC%area_history_secdf(g,n_event) = area
!!$                POPLUC%age_history_secdf(g,n_event) = -1 ! will be incremented to 
!!$                ! 0 at end of 'i' loop
!!$                i = i-1 ! repeat this run through i loop as we have shifted list downward
!!$             ENDIF
!!$
!!$          ENDIF  ! HARVEST
!!$       ENDIF  ! area>0
!!$       i = i+ 1
!!$    ENDDO




  END SUBROUTINE INCREMENT_AGE
  !*******************************************************************************
  SUBROUTINE POPLUCStep(POPLUC,year)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    INTEGER(i4b), INTENT(IN) :: year
    INTEGER(i4b) :: g,j

    POPLUC%it = POPLUC%it + 1


    IF (year.lt.POPLUC%thisyear) THEN

       DO g = 1,POPLUC%np
          ! CALL calculate_weights(POPLUC, g)
          CALL increment_age(POPLUC,g)
       ENDDO

    ELSE

       DO g = 1,POPLUC%np

        
          if (POPLUC%ptos(g) .gt. 0.0) &
               CALL execute_luc_event('PRIMF','SECDF',POPLUC%ptos(g),g,POPLUC)

          if (POPLUC%ptog(g) .gt. 0.0) &
               CALL execute_luc_event('PRIMF','C3ANN',POPLUC%ptog(g),g,POPLUC)

          if (POPLUC%stop(g) .gt.0.0) &
               CALL execute_luc_event('SECDF','PRIMF',POPLUC%stop(g),g,POPLUC)
          if (POPLUC%stog(g) .gt.0.0) &
               CALL execute_luc_event('SECDF','C3ANN',POPLUC%stog(g),g,POPLUC)

          if (POPLUC%gtop(g) .gt.0.0) &
               CALL execute_luc_event('C3ANN','PRIMF',POPLUC%gtop(g),g,POPLUC)
          if (POPLUC%gtos(g) .gt.0.0) &
               CALL execute_luc_event('C3ANN','SECDF',POPLUC%gtos(g),g,POPLUC)

          POPLUC%frac_forest(g) =  POPLUC%primf(g)+ SUM(POPLUC%freq_age_secondary(g,:))
!!$          DO j=1,POPLUC%n_event(g)
!!$             POPLUC%frac_forest(g) =  POPLUC%frac_forest(g) + &
!!$                  POPLUC%area_history_secdf(g,j)
!!$
!!$          ENDDO

          ! CALL calculate_weights(POPLUC, g)
          CALL increment_age(POPLUC,g)

       ENDDO

    ENDIF

  END SUBROUTINE POPLUCStep
  !*******************************************************************************
  SUBROUTINE POPLUC_weights_transfer(POPLUC,POP,LUC_EXPT)

    !-------------------------------------------------------------------------------
    ! This subroutine transfers age distributions to POP and redestributes carbon amongst 
    ! pools according to land-use transtions
    !-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    TYPE(POP_TYPE), INTENT(INOUT) :: POP
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    integer:: g, k, j, l
    REAL(dp), DIMENSION(:), ALLOCATABLE:: freq_age

    DO l = 1, POP%np
       pop%pop_grid(l)%freq_age = 0.0
    ENDDO
    ALLOCATE (freq_age(age_max))

    DO g = 1,POPLUC%np
       IF (.NOT.LUC_EXPT%prim_only(g) .and.sum(POPLUC%freq_age_secondary(g,:)).gt.1e-12) THEN
          j = landpt(g)%cstart + 1
          freq_age = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
          DO l = 1, POP%np
             if (POP%Iwood(l).eq.j) then
                ! freq_age seems to be a vector, right side is constant and can
                ! be computed outside the loop
!                pop%pop_grid(l)%freq_age = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
                pop%pop_grid(l)%freq_age = freq_age
             endif
             !write(71,*) POPLUC%it,l,  sum(pop%pop_grid(l)%freq_age)
          ENDDO
       ENDIF
    ENDDO
    DEALLOCATE (freq_age)

!!$    DO g = 1,POPLUC%np
!!$       IF (.NOT.LUC_EXPT%prim_only(g)) THEN
!!$          j = landpt(g)%cstart + 1
!!$          DO l = 1, POP%np
!!$             if (POP%Iwood(l).eq.j .and.sum(POPLUC%freq_age_secondary(g,:)).gt.1e-12 ) then
!!$                pop%pop_grid(l)%freq_age = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
!!$             endif
!!$             !write(71,*) POPLUC%it,l,  sum(pop%pop_grid(l)%freq_age)
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO

!!$write(69,991) POPLUC%POPLUC_Grid(1)%freq_age_secondary 
!!$991  format(1166(e12.4,2x)) 

  END SUBROUTINE POPLUC_weights_transfer
  !*******************************************************************************
  SUBROUTINE POP_LUC_CASA_transfer(POPLUC,POP,LUC_EXPT,casapool,casabal)

    !-------------------------------------------------------------------------------
    ! This subroutine transfers age distributions to POP and redestributes carbon amongst 
    ! pools according to land-use transtions
    !-------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE(POP_TYPE), INTENT(IN) :: POP
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    TYPE (casa_pool),           INTENT(INOUT) :: casapool
    TYPE (casa_balance),        INTENT(INOUT) :: casabal
    integer:: g, k, j, l, idp, irp, idlu, irlu, ilu
    INTEGER,  PARAMETER :: &
         ptos         =  1, &
         ptog         =  2, &
         stog         =  3, &
         gtos         =  4
    INTEGER,  PARAMETER :: &
         p         =  1, &
         s         =  2, &
         gr         =  3
    REAL(dp) :: dcplant(nLU,nLU,3), dclitter(nLU,nLU,3), dcsoil(nLU,nLU,3)
    REAL(dp) :: dcplant_r(nLU,3), dclitter_r(nLU,3), dcsoil_r(nLU,3)
    REAL(dp) :: dcplant_d(nLU,3), dclitter_d(nLU,3), dcsoil_d(nLU,3)
    REAL(dp) :: dnplant(nLU,nLU,3), dnlitter(nLU,nLU,3), dnsoil(nLU,nLU,3)
    REAL(dp) :: dnplant_r(nLU,3), dnlitter_r(nLU,3), dnsoil_r(nLU,3)
    REAL(dp) :: dnplant_d(nLU,3), dnlitter_d(nLU,3), dnsoil_d(nLU,3)
    REAL(dp) :: dpplant(nLU,nLU,3), dplitter(nLU,nLU,3), dpsoil(nLU,nLU,3)
    REAL(dp) :: dpplant_r(nLU,3), dplitter_r(nLU,3), dpsoil_r(nLU,3)
    REAL(dp) :: dpplant_d(nLU,3), dplitter_d(nLU,3), dpsoil_d(nLU,3)
    REAL(dp) :: dnsoilmin_r(nLU), dclabile_r(nLU)
    REAL(dp) :: dnsoilmin_d(nLU), dclabile_d(nLU)
    REAL(dp) :: dA_r(nLU), dA_d(nLU), dA(nLU), deltaA, dwood_transfer
    REAL(dp), ALLOCATABLE :: dcHarvCLear(:), dcHarv(:), dcClear(:), dcExpand(:), FHarvClear(:)
    REAL(dp) :: kHarvProd(3), kClearProd(3)
    

    ! turnover rates for harvest and clearance products (y-1)
    kHarvProd(1) = 1.0/1.0
    kHarvProd(2) = 1.0/10.0
    kHarvProd(3) = 1.0/100.0
    kClearProd = kHarvProd

    ! zero POPLUC fluxes
    popluc%FtransferNet = 0.0
    popluc%FtransferGross = 0.0
    popluc%FHarvest = 0.0
    popluc%FClearance = 0.0

    ! local variable for storing sum of biomass change due to 
    !secondary harvest, clearance and expansion, and secondary forest
    !  harvest and clearance fluxes
    Allocate(FHarvClear(POPLUC%np))
    Allocate(dcHarvClear(POPLUC%np))
    Allocate(dcHarv(POPLUC%np))
    Allocate(dcClear(POPLUC%np))
    Allocate(dcExpand(POPLUC%np))
    FHarvClear = 0.0
    dcHarvClear = 0.0
    dcHarv = 0.0
    dcClear = 0.0
    dcExpand = 0.0 
    POPLUC%FHarvest = 0.0
    POPLUC%FClearance = 0.0
    
    ! Transfer POP age-dependent biomass to POLUC variables (diagnostic only) and
    ! catastrophic mortality in secondary forest tiles to changes in biomass associated
    ! with secondary forest harvest and clearance and expansion.  
    DO g = 1,POPLUC%np

       j = landpt(g)%cstart
       DO l = 1, POP%np
          IF (.NOT.LUC_EXPT%prim_only(g)) THEN
             if (POP%Iwood(l).eq.j+1 .and. & 
                  (patch(j+1)%frac+ POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g)) .gt. 0.0   ) then
                POPLUC%biomass_age_secondary(g,:)=POP%pop_grid(l)%biomass_age
                ! change in biomass density due to secondary forest expansion 
                dcExpand(g) = -(POPLUC%gtos(g)+POPLUC%ptos(g))*casapool%cplant(j+1,2)/ &
                     (patch(j+1)%frac + POPLUC%gtos(g)+POPLUC%ptos(g)-POPLUC%stog(g))
                if (POP%pop_grid(l)%cmass_sum_old .gt. 0.001) then
                   dcHarvClear(g) = -max(POP%pop_grid(l)%cat_mortality/ &
                        (POP%pop_grid(l)%cmass_sum_old + POP%pop_grid(l)%growth) &
                        - 1.0/ disturbance_interval ,0.0) &
                        *casapool%cplant(j+1,2) - dcExpand(g)
                  if (( POPLUC%smharv(g) + POPLUC%stog(g)).gt.1e-10) then
                     dcHarv(g) = dcHarvClear(g) * POPLUC%smharv(g) & 
                          /( POPLUC%smharv(g) + POPLUC%stog(g))
                     dcClear(g) = dcHarvClear(g) * POPLUC%stog(g) & 
                          /( POPLUC%smharv(g) + POPLUC%stog(g))

                     if ((casapool%cplant(j+1,2) + dcExpand(g)+dcClear(g)+dcHarv(g)).gt.0.0) then
                        FHarvClear(g) = -(patch(j+1)%frac *  (dcHarvClear(g) + dcExpand(g)) + &
                             (POPLUC%ptos(g)+POPLUC%gtos(g)-POPLUC%stog(g))* &
                             ( casapool%cplant(j+1,2) + dcHarvClear(g) + dcExpand(g) ) ) 
                        POPLUC%FHarvest(g,2) =  (1.-POPLUC%fracHarvSecResid(g))*FHarvClear(g) &
                             * POPLUC%smharv(g) &
                             /( POPLUC%smharv(g) + POPLUC%stog(g))
                        POPLUC%FClearance(g,2) =  (1.-POPLUC%fracClearResid(g))* FHarvClear(g) * &
                             POPLUC%stog(g) & 
                             /( POPLUC%smharv(g) + POPLUC%stog(g))
                     else
                        dcHarv(g) = 0.0
                        dcClear(g) = 0.0
                        POPLUC%FHarvest(g,2) = 0.0
                        POPLUC%FClearance(g,2) = 0.0
                     endif
                        
                  endif
                endif

             endif
          ENDIF
          if (POP%Iwood(l).eq.j ) then
             POPLUC%biomass_age_primary(g,:)=POP%pop_grid(l)%biomass_age
             POPLUC%freq_age_primary(g,:)=POP%pop_grid(l)%freq_age
          endif
       ENDDO
    ENDDO


    
    ! Calculate Carbon Pool Transfers
    DO g = 1,POPLUC%np
       j = landpt(g)%cstart
       l = landpt(g)%cend
       dA_r = 0.0
       dA_d = 0.0
       dA = 0.0
       dcsoil = 0.0
       dcplant = 0.0
       dclitter = 0.0
       dcsoil_r = 0.0
       dcplant_r = 0.0
       dclitter_r = 0.0
       dcsoil_r = 0.0
       dcplant_r = 0.0
       dclitter_r = 0.0
       dclabile_r = 0.0

       dcsoil_d = 0.0
       dcplant_d = 0.0
       dclitter_d = 0.0
       dcsoil_d = 0.0
       dcplant_d = 0.0
       dclitter_d = 0.0
       dclabile_d = 0.0

       dnsoil = 0.0
       dnplant = 0.0
       dnlitter = 0.0
       dnsoil_r = 0.0
       dnplant_r = 0.0
       dnlitter_r = 0.0
       dnsoil_r = 0.0
       dnplant_r = 0.0
       dnlitter_r = 0.0
       dnsoilmin_r = 0.0
    
       dnlitter = 0.0
       dnsoil_d = 0.0
       dnplant_d = 0.0
       dnlitter_d = 0.0
       dnsoil_d = 0.0
       dnplant_d = 0.0
       dnlitter_d = 0.0
       dnsoilmin_d = 0.0

       dpsoil = 0.0
       dpplant = 0.0
       dplitter = 0.0
       dpsoil_r = 0.0
       dpplant_r = 0.0
       dplitter_r = 0.0
       dpsoil_r = 0.0
       dpplant_r = 0.0
       dplitter_r = 0.0
       dwood_transfer = 0.0

       dpsoil_d = 0.0
       dpplant_d = 0.0
       dplitter_d = 0.0
       dpsoil_d = 0.0
       dpplant_d = 0.0
       dplitter_d = 0.0

       IF (.NOT.LUC_EXPT%prim_only(g)) THEN
          DO k = 1,nTrans
             if (k==1) then
                deltaA = POPLUC%ptos(g)        ! transition area
                idp = j                        ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = p                       ! donor land use index
                irlu = s                       ! receiver land use index
             elseif (k==2) then
                deltaA = POPLUC%ptog(g)
                idp = j                        ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = p                       ! donor land use index
                irlu = gr                       ! receiver land use index
             elseif (k==3) then
                deltaA = POPLUC%stog(g)
                idp = j+1                      ! donor patch index
                irp = j+2                      ! receiver patch index
                idlu = s                       ! donor land use index
                irlu = gr                       ! receiver land use index
             elseif(k==4) then
                deltaA = POPLUC%gtos(g)
                idp = j+2                      ! donor patch index
                irp = j+1                      ! receiver patch index
                idlu = gr                       ! donor land use index
                irlu = s                       ! receiver land use index
             endif
             dcsoil   = 0
             dclitter = 0
             dcplant  = 0
             dnsoil   = 0
             dnlitter = 0
             dnplant  = 0
             dpsoil   = 0
             dplitter = 0
             dpplant  = 0
             dwood_transfer = 0.0

             ! transfer fluxes
             if (deltaA.gt.0.0) then
                ! all soil pools
                dcsoil(irlu,idlu,:) = deltaA*casapool%csoil(idp,:)
                dcsoil_r(irlu,:) =  dcsoil_r(irlu,:) + deltaA*casapool%csoil(idp,:)
                dcsoil_d(idlu,:) =  dcsoil_d(irlu,:) - deltaA*casapool%csoil(idp,:)

                dnsoil(irlu,idlu,:) = deltaA*casapool%nsoil(idp,:)
                dnsoil_r(irlu,:) =  dnsoil_r(irlu,:) + deltaA*casapool%nsoil(idp,:)
                dnsoil_d(idlu,:) =  dnsoil_d(irlu,:) - deltaA*casapool%nsoil(idp,:)

                dnsoilmin_r(irlu) = dnsoilmin_r(irlu) + deltaA*casapool%nsoilmin(idp)
                dnsoilmin_d(irlu) = dnsoilmin_d(irlu) - deltaA*casapool%nsoilmin(idp)

                dpsoil(irlu,idlu,:) = deltaA*casapool%psoil(idp,:)
                dpsoil_r(irlu,:) =  dpsoil_r(irlu,:) + deltaA*casapool%psoil(idp,:)
                dpsoil_d(idlu,:) =  dpsoil_d(irlu,:) - deltaA*casapool%psoil(idp,:)

                ! microbial litter 
                dclitter(irlu,idlu,1) = deltaA*casapool%clitter(idp,1)
                dclitter_r(irlu,1) =  dclitter_r(irlu,1) + dclitter(irlu,idlu,1)
                dclitter_d(idlu,1) =  dclitter_d(irlu,1) - dclitter(irlu,idlu,1)

                dnlitter(irlu,idlu,1) = deltaA*casapool%nlitter(idp,1)
                dnlitter_r(irlu,1) =  dnlitter_r(irlu,1) + dnlitter(irlu,idlu,1)
                dnlitter_d(idlu,1) =  dnlitter_d(irlu,1) - dnlitter(irlu,idlu,1)

                dplitter(irlu,idlu,1) = deltaA*casapool%plitter(idp,1)
                dplitter_r(irlu,1) =  dplitter_r(irlu,1) + dplitter(irlu,idlu,1)
                dplitter_d(idlu,1) =  dplitter_d(irlu,1) - dplitter(irlu,idlu,1)

                ! CWD : donor pool inherits below-ground woody biomass
                if (idlu == s .and. casapool%cplant(idp,2).gt.1e-5 ) then
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) + &
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3) + & 
                        POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))  &
                        *casapool%nplant(idp,2)/ casapool%cplant(idp,2)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3) + & 
                       POPLUC%fracClearResid(g)*POPLUC%Fclearance(g,2)/ &
                        (1.-POPLUC%fracClearResid(g))  &
                          *casapool%pplant(idp,2)/ casapool%cplant(idp,2)
                elseif (idlu == p .and. irlu == s) then
                   dclitter(irlu,idlu,3) = deltaA*(casapool%clitter(idp,3) + &
                         POPLUC%fracHarvResid(g)*casapool%cplant(idp,2) )
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracHarvResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + & 
                       POPLUC%fracHarvResid(g) *casapool%pplant(idp,2) )
                elseif (idlu == p .and. irlu == gr) then
                   dclitter(irlu,idlu,3) = deltaA*(casapool%clitter(idp,3) + &
                         POPLUC%fracClearResid(g)*casapool%cplant(idp,2) )
                   dnlitter(irlu,idlu,3) = deltaA*(casapool%nlitter(idp,3) + &
                        POPLUC%fracClearResid(g)*casapool%nplant(idp,2) )
                   dplitter(irlu,idlu,3) = deltaA*(casapool%plitter(idp,3) + & 
                       POPLUC%fracClearResid(g) *casapool%pplant(idp,2) )
                else
                   dclitter(irlu,idlu,3) = deltaA*casapool%clitter(idp,3) 
                   dnlitter(irlu,idlu,3) = deltaA*casapool%nlitter(idp,3)
                   dplitter(irlu,idlu,3) = deltaA*casapool%plitter(idp,3)
                endif
                dclitter_r(irlu,3) =  dclitter_r(irlu,3) + dclitter(irlu,idlu,3)
                dclitter_d(idlu,3) =  dclitter_d(irlu,3) - dclitter(irlu,idlu,3)

                dnlitter_r(irlu,3) =  dnlitter_r(irlu,3) + dnlitter(irlu,idlu,3)
                dnlitter_d(idlu,3) =  dnlitter_d(irlu,3) - dnlitter(irlu,idlu,3)

                dplitter_r(irlu,3) =  dplitter_r(irlu,3) + dplitter(irlu,idlu,3)
                dplitter_d(idlu,3) =  dplitter_d(irlu,3) - dplitter(irlu,idlu,3)

                ! fine structural litter: donor pool inherits leaves and fine roots
                dclitter(irlu,idlu,2) = deltaA*(casapool%clitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dclitter_r(irlu,2) =  dclitter_r(irlu,2) + dclitter(irlu,idlu,2)
                dclitter_d(idlu,2) =  dclitter_d(irlu,2) - dclitter(irlu,idlu,2)

                dnlitter(irlu,idlu,2) = deltaA*(casapool%nlitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dnlitter_r(irlu,2) =  dnlitter_r(irlu,2) + dnlitter(irlu,idlu,2)
                dnlitter_d(idlu,2) =  dnlitter_d(irlu,2) - dnlitter(irlu,idlu,2)

                dplitter(irlu,idlu,2) = deltaA*(casapool%plitter(idp,2) + casapool%cplant(idp,1) + &
                     casapool%cplant(idp,3))
                dplitter_r(irlu,2) =  dplitter_r(irlu,2) + dplitter(irlu,idlu,2)
                dplitter_d(idlu,2) =  dplitter_d(irlu,2) - dplitter(irlu,idlu,2)

                ! biomass: no biomass inherited
                dclabile_r(irlu) =  dclabile_r(irlu) + deltaA*casapool%clabile(idp)
                dclabile_d(irlu) =  dclabile_d(irlu) - deltaA*casapool%clabile(idp)

                dcplant(irlu,idlu,:) = 0.0
                dcplant_r(irlu,:) =  dcplant_r(irlu,:) + dcplant(irlu,idlu,:)
                dcplant_d(idlu,:) =  dcplant_d(irlu,:) - dcplant(irlu,idlu,:)

                dnplant(irlu,idlu,:) = 0.0
                dnplant_r(irlu,:) =  dnplant_r(irlu,:) + dnplant(irlu,idlu,:)
                dnplant_d(idlu,:) =  dnplant_d(irlu,:) - dnplant(irlu,idlu,:)

                dpplant(irlu,idlu,:) = 0.0
                dpplant_r(irlu,:) =  dpplant_r(irlu,:) + dpplant(irlu,idlu,:)
                dpplant_d(idlu,:) =  dpplant_d(irlu,:) - dpplant(irlu,idlu,:)

                ! Gross Transfer Flux
                popluc%FtransferGross(g,k) = &
                     sum(dcsoil(irlu,idlu,:) + dclitter(irlu,idlu,:) +  dcplant(irlu,idlu,:))

                ! Harvest Flux
                ! primary to secondary forest
                if (idlu==p .and. irlu ==s) then
                   popluc%FHarvest(g,idlu) =  (1.0 -POPLUC%fracHarvResid(g)) &
                        *casapool%cplant(idp,2)*deltaA
                endif
                ! Clearance Flux
                ! primary forest to grass
                if ((idlu==p) .and. irlu ==gr) then
                   popluc%FClearance(g,idlu) =  (1.0 -POPLUC%fracClearResid(g)) &
                        * casapool%cplant(idp,2)*deltaA
                endif
                  
                ! transition area
                dA_r(irlu) = dA_r(irlu) + deltaA
                dA_d(idlu) = dA_d(idlu) - deltaA
             endif
          ENDDO  ! ntrans

          DO ilu=1,nLU
             ! update pools
             irp = ilu + j -1
             dwood_transfer = 0.0
             dA(ilu) = dA_r(ilu) + dA_d(ilu)
             ! Net Transfer Flux
             popluc%FtransferNet(g,ilu) = sum(dcsoil_r(ilu,:) + dcsoil_d(ilu,:)) + &
                  sum(dclitter_r(ilu,:) + dclitter_d(ilu,:))  + &
                  sum(dcplant_r(ilu,:)  + dcplant_d(ilu,:))
             
             if (ilu ==2) then
                dclitter_r(ilu,3) = dclitter_r(irlu,3) + &
                     POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g)) 
                if (casapool%cplant(irp,2) .gt. 1.e-5) then
                   dnlitter_r(ilu,3) = dnlitter_r(irlu,3) +  &
                     POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g)) &
                     *casapool%nplant(irp,2)/ casapool%cplant(irp,2)
                   dplitter_r(ilu,3) = dplitter_r(irlu,3) +  &
                     POPLUC%fracHarvSecResid(g)*POPLUC%FHarvest(g,2)/(1.0 -POPLUC%fracHarvSecResid(g)) &
                     *casapool%pplant(irp,2)/ casapool%cplant(irp,2)
                endif
             endif

                
            ! if (abs(dA(ilu)).gt.0.0 .and.(patch(irp)%frac+dA(ilu)).gt.1.e-6  ) then
             if ((patch(irp)%frac+dA(ilu)).gt.1.e-6  ) then

                casapool%nsoilmin(irp) = casapool%nsoilmin(irp) +  &
                     (dnsoilmin_r(ilu) - casapool%nsoilmin(irp)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%clabile(irp) = casapool%clabile(irp) +  &
                     (dclabile_r(ilu) - casapool%clabile(irp)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%csoil(irp,:)  = casapool%csoil(irp,:) + &
                     (dcsoil_r(ilu,:)-casapool%csoil(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%nsoil(irp,:)  = casapool%nsoil(irp,:) + &
                     (dnsoil_r(ilu,:)-casapool%nsoil(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%psoil(irp,:)  = casapool%psoil(irp,:) + &
                     (dpsoil_r(ilu,:)-casapool%psoil(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))


                casapool%clitter(irp,:)  = casapool%clitter(irp,:) + &
                     (dclitter_r(ilu,:)-casapool%clitter(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%nlitter(irp,:)  = casapool%nlitter(irp,:) + &
                     (dnlitter_r(ilu,:)-casapool%nlitter(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%plitter(irp,:)  = casapool%plitter(irp,:) + &
                     (dplitter_r(ilu,:)-casapool%plitter(irp,:)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,1)  = casapool%cplant(irp,1) + &
                     (dcplant_r(ilu,1)-casapool%cplant(irp,1)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,1)  = casapool%nplant(irp,1) + &
                     (dnplant_r(ilu,1)-casapool%nplant(irp,1)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,1)  = casapool%pplant(irp,1) + &
                     (dpplant_r(ilu,1)-casapool%pplant(irp,1)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,3)  = casapool%cplant(irp,3) + &
                     (dcplant_r(ilu,3)-casapool%cplant(irp,3)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,3)  = casapool%nplant(irp,3) + &
                     (dnplant_r(ilu,3)-casapool%nplant(irp,3)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,3)  = casapool%pplant(irp,3) + &
                     (dpplant_r(ilu,3)-casapool%pplant(irp,3)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%cplant(irp,2)  = casapool%cplant(irp,2) + &
                     (dcplant_r(ilu,2)-casapool%cplant(irp,2)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%nplant(irp,2)  = casapool%nplant(irp,2) + &
                     (dnplant_r(ilu,2)-casapool%nplant(irp,2)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))

                casapool%pplant(irp,2)  = casapool%pplant(irp,2) + &
                     (dpplant_r(ilu,2)-casapool%pplant(irp,2)*(dA(ilu) - dA_d(ilu)))/(patch(irp)%frac+dA(ilu))




                ! account here for change in secondary forest biomass density due to:
                !  harvest and clearing, as well as increase in (below ground) CWD 
                ! where secondary forest harvest occurs
                if (ilu .eq.s .and. (casapool%cplant(irp,2) +(dcHarv(g)+dCClear(g)) ).gt.0.0  &
                   .and. casapool%cplant(irp,2).gt.1.e-10  ) then
                   
                 
                   casapool%nplant(irp,2) = casapool%nplant(irp,2) + (dcHarv(g)+dcClear(g)) &
                        * casapool%nplant(irp,2)/casapool%cplant(irp,2)
                   casapool%pplant(irp,2) = casapool%pplant(irp,2) + (dcHarv(g)+dcClear(g)) &
                        * casapool%pplant(irp,2)/casapool%cplant(irp,2)
                   casapool%cplant(irp,2) = casapool%cplant(irp,2) + (dcHarv(g)+dcClear(g))

               elseif (ilu .eq.s .and. (casapool%cplant(irp,2) +(dcHarv(g)+dCClear(g)) ).le.0.0  ) then

                  POPLUC%FHarvest(g,2) = 0.0
                  POPLUC%FClearance(g,2) = 0.0
               endif


            endif

          
            casapool%Cclear(irp,1) =  casapool%Cclear(irp,1) + popluc%FClearance(g,ilu)
            casapool%Chwp(irp,1) =  casapool%Chwp(irp,1) + popluc%FHarvest(g,ilu)
         

         ENDDO


          POPLUC%FNEP(g,:) = casabal%Fcneeyear(j:l)*patch(j:l)%frac  ! note NEE = NEP here (+ve into surface)
          DO ilu=1,nLU
             ! update area weights
             irp = ilu + j -1
             ! Net Area Change
             patch(irp)%frac= max(patch(irp)%frac + dA_r(ilu) + dA_d(ilu), 0.0)

          ENDDO
          ! pools in gC per m2 of gridcell
          ! NEP in g C y-1 per m2 of gridcell

          POPLUC%csoil(g,:) = sum(casapool%csoil(j:l,:),2)*patch(j:l)%frac
          POPLUC%clitt(g,:) = sum(casapool%clitter(j:l,:),2)*patch(j:l)%frac
          POPLUC%cbiomass(g,:) = sum(casapool%cplant(j:l,:),2)*patch(j:l)%frac
      
          POPLUC%primf(g) = patch(j)%frac
          POPLUC%secdf(g) = patch(j+1)%frac
          POPLUC%grass(g) = patch(j+2)%frac

       ELSE

          POPLUC%csoil(g,:) = 0.0
          POPLUC%clitt(g,:) = 0.0
          POPLUC%cbiomass(g,:) = 0.0
          POPLUC%FNEP(g,:) = 0.0

!!$          POPLUC%csoil(g,1:l-j+1) = sum(casapool%csoil(j:l,:),2)*patch(j:l)%frac
!!$          POPLUC%clitt(g,1:l-j+1) = sum(casapool%clitter(j:l,:),2)*patch(j:l)%frac
!!$          POPLUC%cbiomass(g,1:l-j+1) = sum(casapool%cplant(j:l,:),2)*patch(j:l)%frac
!!$          POPLUC%FNEP(g,1:l-j+1) =  casabal%Fcneeyear(j:l)*patch(j:l)%frac

          POPLUC%csoil(g,1) = sum(casapool%csoil(j,:))*patch(j)%frac
          POPLUC%clitt(g,1) = sum(casapool%clitter(j,:))*patch(j)%frac
          POPLUC%cbiomass(g,1) = sum(casapool%cplant(j,:))*patch(j)%frac
          POPLUC%FNEP(g,1) =  casabal%Fcneeyear(j)*patch(j)%frac
          POPLUC%primf(g) = patch(j)%frac
          POPLUC%secdf(g) = 0.0

          if (POPLUC%grass(g).gt.0.0 .and. l.eq.j+1) then
             POPLUC%csoil(g,3) = sum(casapool%csoil(l,:))*patch(l)%frac
             POPLUC%clitt(g,3) = sum(casapool%clitter(l,:))*patch(l)%frac
             POPLUC%cbiomass(g,3) = sum(casapool%cplant(l,:))*patch(l)%frac
             POPLUC%FNEP(g,3) = casabal%Fcneeyear(l)*patch(l)%frac
             POPLUC%grass(g) = patch(l)%frac
          endif

       ENDIF

    
       POPLUC%HarvProdLoss(g,:) = kHarvProd * POPLUC%HarvProd(g,:)
       POPLUC%ClearProdLoss(g,:) = kClearProd * POPLUC%ClearProd(g,:)
       
       DO j=1,3
          POPLUC%HarvProd(g,j) = POPLUC%HarvProd(g,j) + &
               POPLUC%fracHarvProd(g,j)*sum(POPLUC%FHarvest(g,:)) - POPLUC%HarvProdLoss(g,j) 

          POPLUC%ClearProd(g,j) = POPLUC%ClearProd(g,j) + &
               POPLUC%fracClearProd(g,j)*sum(POPLUC%FClearance(g,:)) - POPLUC%ClearProdLoss(g,j) 
       
       ENDDO
    

    ENDDO
991 format(1166(e14.7,2x)) 
!!$    write(615,991) casabal%Fcneeyear(2), casabal%dcdtyear(2), casabal%dcdtyear(2) - casabal%Fcneeyear(2)
!!$   ! write(615,991) casabal%dcdtyear
    casapool%ctot = sum(casapool%cplant,2)+sum(casapool%clitter,2)+sum(casapool%csoil,2)
    casabal%cplantlast  = casapool%cplant
    casabal%clabilelast = casapool%clabile
    casabal%clitterlast = casapool%clitter
    casabal%csoillast   = casapool%csoil

  END SUBROUTINE POP_LUC_CASA_transfer
  !*******************************************************************************
  SUBROUTINE POPLUC_Init(POPLUC,LUC_EXPT,np)


    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC
    TYPE (LUC_EXPT_TYPE), INTENT(IN) :: LUC_EXPT
    INTEGER(i4b), INTENT(IN) :: np
    INTEGER(i4b) :: g

    CALL alloc_POPLUC(POPLUC,np)

    POPLUC%it = 0
    POPLUC%np = np

    CALL ZeroPOPLUC(POPLUC)

    DO g = 1,np
       POPLUC%fracharvProd(g, 1) = 0.9
       POPLUC%fracharvProd(g, 2) = 0.04
       POPLUC%fracharvProd(g, 3) = 0.06
       
       POPLUC%fracClearProd(g, 1) = 0.597
       POPLUC%fracClearProd(g, 2) = 0.403
       POPLUC%fracClearProd(g, 3) = 0.00
       POPLUC%fracClearResid(g) = 0.33
       POPLUC%fracHarvResid(g) = 0.79
       POPLUC%fracHarvSecResid(g) = 0.81          

       IF (LUC_EXPT%biome(g)==1 .OR. LUC_EXPT%biome(g)==2) THEN        
          ! Tropical Evergreen and Tropical Deciduous
          POPLUC%fracharvProd(g, 1) = 0.9
          POPLUC%fracharvProd(g, 2) = 0.04
          POPLUC%fracharvProd(g, 3) = 0.06

          POPLUC%fracClearProd(g, 1) = 0.597
          POPLUC%fracClearProd(g, 2) = 0.403
          POPLUC%fracClearProd(g, 3) = 0.00

          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvResid(g) = 0.79
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvResid(g) = 0.86
          IF (LUC_EXPT%biome(g)==1) POPLUC%fracHarvSecResid(g) = 0.71
          IF (LUC_EXPT%biome(g)==2) POPLUC%fracHarvSecResid(g) = 0.81

          POPLUC%fracClearResid(g) = 0.33

       ELSEIF (LUC_EXPT%biome(g).GE.4 .OR. LUC_EXPT%biome(g).LE.10) THEN       

          ! Other Forest
          POPLUC%fracharvProd(g, 1) = 0.4
          POPLUC%fracharvProd(g, 2) = 0.24
          POPLUC%fracharvProd(g, 3) = 0.36

          POPLUC%fracClearProd(g, 1) = 0.597
          POPLUC%fracClearProd(g, 2) = 0.2985
          POPLUC%fracClearProd(g, 3) = 0.1045

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvResid(g) = 0.83
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvResid(g) = 0.87
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvResid(g) = 0.78

          IF (LUC_EXPT%ivegp(g)==2) POPLUC%fracHarvSecResid(g) = 0.75
          IF (LUC_EXPT%ivegp(g)==1 .OR. LUC_EXPT%ivegp(g)==3 ) POPLUC%fracHarvSecResid(g) = 0.82
          IF (LUC_EXPT%ivegp(g)==4) POPLUC%fracHarvSecResid(g) = 0.70


          POPLUC%fracClearResid(g) = 0.33

       ELSEIF (LUC_EXPT%biome(g).EQ.3 .OR. LUC_EXPT%biome(g).GE.11) THEN
          ! savanna and shrub
          POPLUC%fracharvProd(g, 1) = 1.0
          POPLUC%fracharvProd(g, 2) = 0.00
          POPLUC%fracharvProd(g, 3) = 0.00

          POPLUC%fracClearProd(g, 1) = 0.8
          POPLUC%fracClearProd(g, 2) = 0.2
          POPLUC%fracClearProd(g, 3) = 0.0


          IF (LUC_EXPT%biome(g).EQ.13 .OR.(LUC_EXPT%biome(g).EQ.14)  ) THEN
             POPLUC%fracHarvResid(g) = 0.78
             POPLUC%fracHarvSecResid(g) = 0.70
             
             
          ELSE
             POPLUC%fracHarvResid(g) = 0.86
             POPLUC%fracHarvSecResid(g) = 0.81
             
          ENDIF
          
          POPLUC%fracClearResid(g) = 0.5

       ENDIF
       
    ENDDO

  END SUBROUTINE POPLUC_Init
  !*******************************************************************************
  SUBROUTINE alloc_POPLUC(POPLUC, arraysize)
    IMPLICIT NONE
    TYPE(POPLUC_TYPE),INTENT(INOUT) :: POPLUC
    INTEGER,            INTENT(IN) :: arraysize

    ALLOCATE(POPLUC%it)
    ALLOCATE(POPLUC%np)
    ALLOCATE(POPLUC%firstyear)
    ALLOCATE(POPLUC%thisyear)
    ALLOCATE(POPLUC%latitude(arraysize))
    ALLOCATE(POPLUC%longitude(arraysize))
    ALLOCATE(POPLUC%n_event(arraysize))
    ALLOCATE(POPLUC%primf(arraysize))
    ALLOCATE(POPLUC%secdf(arraysize))
    ALLOCATE(POPLUC%grass(arraysize))
    ALLOCATE(POPLUC%ptos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%stop(arraysize))
    ALLOCATE(POPLUC%stog(arraysize))
    ALLOCATE(POPLUC%gtop(arraysize))
    ALLOCATE(POPLUC%gtos(arraysize))
    ALLOCATE(POPLUC%ptog(arraysize))
    ALLOCATE(POPLUC%pharv(arraysize))
    ALLOCATE(POPLUC%smharv(arraysize))
    ALLOCATE(POPLUC%syharv(arraysize))
    ALLOCATE(POPLUC%frac_primf(arraysize))
    ALLOCATE(POPLUC%frac_forest(arraysize))
    ALLOCATE(POPLUC%freq_age_primary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%freq_age_secondary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%biomass_age_primary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%biomass_age_secondary(arraysize,AGE_MAX))
    ALLOCATE(POPLUC%age_history_secdf(arraysize,LENGTH_SECDF_HISTORY))
    ALLOCATE(POPLUC%area_history_secdf(arraysize,LENGTH_SECDF_HISTORY))
    ALLOCATE(POPLUC%FNEP(arraysize,nLU))
    ALLOCATE(POPLUC%Clitt(arraysize,nLU))
    ALLOCATE(POPLUC%Csoil(arraysize,nLU)) 
    ALLOCATE(POPLUC%Cbiomass(arraysize,nLU)) 
    ALLOCATE(POPLUC%FHarvest(arraysize,nLU))
    ALLOCATE(POPLUC%FClearance(arraysize,nLU))
    ALLOCATE(POPLUC%FTransferNet(arraysize,nLU)) 
    ALLOCATE(POPLUC%FTransferGross(arraysize,nTrans))
    ALLOCATE(POPLUC%HarvProd(arraysize,3))
    ALLOCATE(POPLUC%ClearProd(arraysize,3))
    ALLOCATE(POPLUC%HarvProdLoss(arraysize,3))
    ALLOCATE(POPLUC%ClearProdLoss(arraysize,3))
    ALLOCATE(POPLUC%fracHarvProd(arraysize,3))
    ALLOCATE(POPLUC%fracClearProd(arraysize,3))
    ALLOCATE(POPLUC%fracHarvResid(arraysize))
    ALLOCATE(POPLUC%fracHarvSecResid(arraysize))
    ALLOCATE(POPLUC%fracClearResid(arraysize))

  END SUBROUTINE alloc_POPLUC
  !*******************************************************************************
  ! Exponential distribution
  ! Returns probability of a given time-between-events (x)
  ! Given a Poisson process with expected frequency (events per unit time) lambda
  ! Reference: http://en.wikipedia.org/wiki/Exponential_distribution
  ! Use to determine average age (x, years) of patches with a given random disturbance
  ! frequency lambda (disturbances per year)

  REAL(dp) FUNCTION REALExponential(lambda, x)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) ::  x
    REAL(dp), INTENT(IN) ::  lambda

    IF (x.LT.0) THEN ! Shouldn't happen but ...
       REALExponential=0.0
    ELSE
       REALExponential=lambda*EXP(-lambda*x)
    ENDIF

  END FUNCTION REALExponential


  !*******************************************************************************
  SUBROUTINE WRITE_LUC_OUTPUT_NC ( POPLUC, ctime, FINAL )

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    LOGICAL, INTENT(IN)    :: FINAL
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, hist_ID, t_ID, nLU_ID, nTrans_ID, i, mp
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp,t)
    CHARACTER(len=20),DIMENSION(13):: A1
    ! 2 dim integer arrays (mp,t)
    CHARACTER(len=20),DIMENSION(1):: AI1
    ! 3 dim real arrays (mp,age_max,t)
    CHARACTER(len=25),DIMENSION(4) :: A2
    ! 3 dim real arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: A3
    ! 3 dim integer arrays (mp,LENGTH_SECDF_HISTORY,t)
    CHARACTER(len=20),DIMENSION(1) :: AI3
    ! 3 dim real arrays (mp,nLU,t)
    CHARACTER(len=20),DIMENSION(11) :: A4
    ! 3 dim real arrays (mp,nTrans,t)
    CHARACTER(len=20),DIMENSION(1) :: A5

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1)),VIDI1(SIZE(AI1))
    INTEGER, SAVE :: VID2(SIZE(A2)),VID3(SIZE(A3)),VIDI3(SIZE(AI3))
    INTEGER, SAVE :: VID4(SIZE(A4)),VID5(SIZE(A5))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    CHARACTER(len=50) :: RecordDimName
    REAL, ALLOCATABLE :: freq_age_secondary(:,:)
    INTEGER :: g
    LOGICAL :: put_age_vars
    mp = POPLUC%np
    put_age_vars=.FALSE.
    allocate(freq_age_secondary(mp,age_max))

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'ptos'
    A1(5) = 'ptog'
    A1(6) = 'stog'
    A1(7) = 'gtop'
    A1(8) = 'gtos'
    A1(9) = 'frac_primf'
    A1(10) = 'frac_forest'
    A1(11) = 'pharv'
    A1(12) = 'smharv'
    A1(13) = 'syharv'
    

    AI1(1) = 'n_event'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'
    A2(3) = 'biomass_age_primary'
    A2(4) = 'biomass_age_secondary'

    A3(1) = 'area_history_secdf'

    AI3(1) = 'age_history_secdf'

    A4(1) = 'FHarvest'
    A4(2) = 'FClearance'
    A4(3) = 'FNEP'
    A4(4) = 'CLitt'
    A4(5) = 'CSoil'
    A4(6) = 'CBiomass'
    A4(7) = 'FTransferNet'
    A4(8) = 'HarvProd'
    A4(9) = 'ClearProd'
    A4(10) = 'HarvProdLoss'
    A4(11) = 'ClearProdLoss'

    A5(1) = 'FTransferGross'

    DO g=1,mp
       if (sum(POPLUC%freq_age_secondary(g,:)).gt.1e-12 ) then
          freq_age_secondary(g,:) = POPLUC%freq_age_secondary(g,:)/sum(POPLUC%freq_age_secondary(g,:))
       else
          freq_age_secondary(g,:) = 0.0
       endif
    ENDDO


    CNT = CNT + 1

    IF ( CALL1 ) THEN
       ! Get File-Name

       WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
          WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
       ENDIF
       fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//&
            TRIM(dum)//'_LUC_out.nc'

       ! Create NetCDF file:
       STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       ! Put the file in define mode:
       STATUS = NF90_redef(FILE_ID)

       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
       STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

       ! Define dimensions:
       ! Land (number of points)
       STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mage' , age_max , age_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mhist',LENGTH_SECDF_HISTORY , hist_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mLU', nLU , nLU_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_def_dim(FILE_ID, 'mTrans',nTrans , nTrans_ID)
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

       DO i = 1, SIZE(AI1)
          STATUS = NF90_def_var(FILE_ID,TRIM(AI1(i)) ,NF90_INT,(/land_ID,t_ID/),VIDI1(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO
       if(put_age_vars) then
          DO i = 1, SIZE(A2)
             STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,age_ID,t_ID/),VID2(i))
             IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
          END DO
       endif
!!$
!!$       DO i = 1, SIZE(A3)
!!$          STATUS = NF90_def_var(FILE_ID,TRIM(A3(i)) ,NF90_FLOAT,(/land_ID,hist_ID,t_ID/),VID3(i))
!!$          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$       END DO
!!$
!!$       DO i = 1, SIZE(A3)
!!$          STATUS = NF90_def_var(FILE_ID,TRIM(AI3(i)) ,NF90_INT,(/land_ID,hist_ID,t_ID/),VIDI3(i))
!!$          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
!!$       END DO

       DO i = 1, SIZE(A4)
          STATUS = NF90_def_var(FILE_ID,TRIM(A4(i)) ,NF90_FLOAT,(/land_ID,nLU_ID,t_ID/),VID4(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       DO i = 1, SIZE(A5)
          STATUS = NF90_def_var(FILE_ID,TRIM(A5(i)) ,NF90_FLOAT,(/land_ID,ntrans_ID,t_ID/),VID5(i))
          IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       END DO

       ! End define mode:
       STATUS = NF90_enddef(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


       ! PUT LAT / LON ( mp )
       STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), POPLUC%latitude )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), POPLUC%longitude )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

       CALL1 = .FALSE.

    ENDIF ! CALL1



    ! TIME  ( t )
    STATUS = NF90_PUT_VAR(FILE_ID, VIDtime, ctime, start=(/ CNT /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), POPLUC%primf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), POPLUC%secdf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), POPLUC%grass,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 4), POPLUC%ptos,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 5), POPLUC%ptog,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 6), POPLUC%stog,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 7), POPLUC%gtop,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 8), POPLUC%gtos,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 9), POPLUC%frac_primf,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1(10), POPLUC%frac_forest,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 11), POPLUC%pharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 12), POPLUC%smharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 13), POPLUC%syharv,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)



    STATUS = NF90_PUT_VAR(FILE_ID, VIDI1(1), POPLUC%n_event,        start=(/ 1, CNT /), count=(/ mp, 1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)


    if (put_age_vars) then
       ! PUT 3D VARS ( mp, mage, t )
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), POPLUC%freq_age_primary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),freq_age_secondary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(3), POPLUC%biomass_age_primary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
       STATUS = NF90_PUT_VAR(FILE_ID, VID2(4), POPLUC%biomass_age_secondary,   &
            start=(/ 1,1,CNT /), count=(/ mp,age_max,1 /) )
       IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    endif

!!$    ! PUT 3D VARS ( mp, mage, t )
!!$    STATUS = NF90_PUT_VAR(FILE_ID, VID3(1), POPLUC%area_history_secdf,   &
!!$         start=(/ 1,1,CNT /), count=(/ mp,LENGTH_SECDF_HISTORY,1 /) )
!!$    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
!!$    STATUS = NF90_PUT_VAR(FILE_ID, VIDI3(1), POPLUC%age_history_secdf,   &
!!$         start=(/ 1,1,CNT /), count=(/ mp,LENGTH_SECDF_HISTORY,1 /) )
!!$    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nLU, t )

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(1), POPLUC%FHarvest,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(2), POPLUC%FClearance,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(3), POPLUC%FNEP,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(4), POPLUC%Clitt,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(5), POPLUC%CSoil,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(6), POPLUC%CBiomass,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%FTransferNet,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%HarvProd,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%ClearProd,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%HarvProdLoss,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID4(7), POPLUC%ClearProdLoss,   &
         start=(/ 1,1,CNT /), count=(/ mp,nLU,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, nTrans, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID5(1), POPLUC%FTransferGross,   &
         start=(/ 1,1,CNT /), count=(/ mp,nTrans,1 /) )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    IF ( FINAL ) THEN
       ! Close NetCDF file:
       STATUS = NF90_close(FILE_ID)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       WRITE(*,*) " POPLUC Output written to ",fname
    ENDIF

  END SUBROUTINE WRITE_LUC_OUTPUT_NC
  !************************************************************************************************************************************
  SUBROUTINE WRITE_LUC_RESTART_NC ( POPLUC, ctime )

    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf

    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(IN) :: POPLUC
    INTEGER, INTENT(IN)    :: ctime

    INTEGER   :: STATUS
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20),DIMENSION(5):: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25),DIMENSION(2) :: A2
   

    INTEGER, SAVE :: VIDtime, VID0(SIZE(A0)),VID1(SIZE(A1))
    INTEGER, SAVE :: VID2(SIZE(A2))
    INTEGER, SAVE :: FILE_ID, CNT = 0
    CHARACTER(len=50) :: RecordDimName
    INTEGER :: g
   
    mp = POPLUC%np

    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'HarvProd'
    A1(5) = 'ClearProd'

    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'

    

    ! Get File-Name

    WRITE( dum, FMT="(I4,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    IF (CABLE_USER%YEARSTART.lt.1000.and.CABLE_USER%YEAREND.lt.1000) THEN
       WRITE( dum, FMT="(I3,'_',I3)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    ELSEIF (CABLE_USER%YEARSTART.lt.1000) THEN
       WRITE( dum, FMT="(I3,'_',I4)")CABLE_USER%YEARSTART,CABLE_USER%YEAREND
    ENDIF
    fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//'LUC_rst.nc'

    ! Create NetCDF file:
    STATUS = NF90_create(fname, NF90_CLOBBER, FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! Put the file in define mode:
    STATUS = NF90_redef(FILE_ID)

    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "StartYear", CABLE_USER%YEARSTART )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "EndYear"  , CABLE_USER%YEAREND   )
    STATUS = NF90_PUT_ATT( FILE_ID, NF90_GLOBAL, "RunIden"  , CABLE_USER%RunIden   )

    ! Define dimensions:
    ! Land (number of points)
    STATUS = NF90_def_dim(FILE_ID, 'land'   , mp     , land_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    STATUS = NF90_def_dim(FILE_ID, 'mage' , age_max , age_ID)
   

    ! Define variables

    DO i = 1, SIZE(A0)
       STATUS = NF90_def_var(FILE_ID,TRIM(A0(i)) ,NF90_FLOAT,(/land_ID/),VID0(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    DO i = 1, SIZE(A1)
       STATUS = NF90_def_var(FILE_ID,TRIM(A1(i)) ,NF90_FLOAT,(/land_ID/),VID1(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO


    DO i = 1, SIZE(A2)
       STATUS = NF90_def_var(FILE_ID,TRIM(A2(i)) ,NF90_FLOAT,(/land_ID,age_ID/),VID2(i))
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    END DO

    ! End define mode:
    STATUS = NF90_enddef(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

    ! PUT LAT / LON ( mp )
    STATUS = NF90_PUT_VAR(FILE_ID, VID0(1), POPLUC%latitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    STATUS = NF90_PUT_VAR(FILE_ID, VID0(2), POPLUC%longitude )
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 2D VARS ( mp, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 1), POPLUC%primf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 2), POPLUC%secdf)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID1( 3), POPLUC%grass)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)

    ! PUT 3D VARS ( mp, mage, t )
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(1), POPLUC%freq_age_primary)
    IF(STATUS /= NF90_NoErr) CALL handle_err(STATUS)
    STATUS = NF90_PUT_VAR(FILE_ID, VID2(2),POPLUC%freq_age_secondary)

    ! Close NetCDF file:
    STATUS = NF90_close(FILE_ID)
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
    WRITE(*,*) " POPLUC Restart written to ",fname


  END SUBROUTINE WRITE_LUC_RESTART_NC
  !*******************************************************************************
  SUBROUTINE READ_LUC_RESTART_NC (POPLUC)
    USE CABLE_COMMON_MODULE, ONLY:  filename, cable_user, HANDLE_ERR
    USE netcdf
    IMPLICIT NONE
    TYPE(POPLUC_TYPE), INTENT(INOUT) :: POPLUC


    INTEGER   :: STATUS, land_dim, mage_dim
    INTEGER   :: land_ID, age_ID, nLU_ID, nTrans_ID, i, mp, FILE_ID, dID
    LOGICAL   :: CASAONLY
    CHARACTER :: CYEAR*4, FNAME*99,dum*50
    LOGICAL, SAVE :: CALL1 = .TRUE.
    REAL , ALLOCATABLE :: TMP(:), TMP2(:,:)

    ! 1 dim arrays (mp )
    CHARACTER(len=20),DIMENSION(2) :: A0
    ! 2 dim real arrays (mp)
    CHARACTER(len=20),DIMENSION(5):: A1
    ! 2 dim real arrays (mp,age_max)
    CHARACTER(len=25),DIMENSION(2) :: A2

    mp = POPLUC%np
    ALLOCATE(tmp(mp))
    ALLOCATE(tmp2(mp,age_max))
    A0(1) = 'latitude'
    A0(2) = 'longitude'

    A1(1) = 'primf'
    A1(2) = 'secdf'
    A1(3) = 'grass'
    A1(4) = 'HarvProd'
    A1(5) = 'ClearProd'


    A2(1) = 'freq_age_primary'
    A2(2) = 'freq_age_secondary'


    !fname = TRIM(filename%path)//'/'//TRIM( cable_user%RunIden )//&
    !     '_LUC_rst.nc'

 fname = TRIM(filename%path)//'/'//TRIM(cable_user%RunIden)//'_'//'LUC_rst.nc'
    STATUS = NF90_OPEN( TRIM(fname), NF90_NOWRITE, FILE_ID )
    IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

   
 ! DIMS
  STATUS = NF90_INQ_DIMID( FILE_ID, 'land', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=land_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

  STATUS = NF90_INQ_DIMID( FILE_ID, 'mage', dID )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
  STATUS = NF90_INQUIRE_DIMENSION( FILE_ID, dID, LEN=mage_dim )
  IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)


    ! READ 1-dimensional fields
    DO i = 1, SIZE(A1)
       STATUS = NF90_INQ_VARID( FILE_ID, A1(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A1(i)))
       CASE ('primf'      ) ; POPLUC%primf       = TMP
       CASE ('secdf'   ) ; POPLUC%secdf   = TMP
       CASE ('grass'  ) ; POPLUC%grass  = TMP

       END SELECT
    END DO

 991  format(1000(e12.4,2x))
    ! READ 2-dimensional fields (mage)
    DO i = 1, SIZE(A2)
       STATUS = NF90_INQ_VARID( FILE_ID, A2(i), dID )
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)
       STATUS = NF90_GET_VAR( FILE_ID, dID, TMP2)
       IF (STATUS /= NF90_noerr) CALL handle_err(STATUS)

       SELECT CASE ( TRIM(A2(i)))
       CASE ('freq_age_primary' ) ; POPLUC%freq_age_primary = TMP2
       CASE ('freq_age_secondary' ) ; POPLUC%freq_age_secondary = TMP2
       END SELECT
    END DO

    STATUS = NF90_CLOSE( FILE_ID )

  END SUBROUTINE READ_LUC_RESTART_NC
 !*******************************************************************************
END MODULE POPLUC_MODULE

