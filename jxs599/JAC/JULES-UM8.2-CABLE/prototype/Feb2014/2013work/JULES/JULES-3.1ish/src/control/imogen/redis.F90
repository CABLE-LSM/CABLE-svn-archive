#if !defined(UM_JULES)
!*********************************************************************
! Routine to redistribute rainfall if maximum precipitation rate is
! exceeded.
!
! Written by Chris Huntingford (September 2001)
!*********************************************************************
  SUBROUTINE REDIS(                                               &
    NSDMAX,STEP_DAY,MAX_PRECIP_RATE,PREC_LOC,N_EVENT_LOCAL,       &
    N_TALLY                                                       &
  )

    INTEGER ::                                                    &
      STEP_DAY,                                                   &
              !IN Calculated number of timesteps per day
      NSDMAX  !IN Maximum possible number of sub-daily timesteps. 

    REAL ::                                                       &
      MAX_PRECIP_RATE,                                            &
              !IN Maximum allowed precip. rate allowed within
              !   each sub-daily timestep (mm/day). (This only
              !   applies when STEP_DAY.GE.2)
      PREC_LOC(NSDMAX),                                           &
              !INOUT Temporary value of rainfall for each gridbox.
      PREC_TOT,                                                   &
              !WORK Total precip during entire day (mm/day
      PREC_TOT_ADJ,                                               &
              !WORK Temporary adjusted total (mm/day)
      PREC_CHANGE,                                                &
              !WORK The amount of precipitation to be
              !     redistributed within day (mm/day)
      EXTRA_PER_REAL
              !WORK Number of extra periods of rainfall
              !     (expressed as a real number) that are required.

    INTEGER ::                                                    &
      N_EVENT_LOCAL(NSDMAX),                                      &
              ! 1: if rains/snows during timestep period
              ! 0: otherwise
      N_TALLY,                                                    &
              !IN Number of precipitation periods before
              !   redistribution. 
      I,                                                          &
              ! WORK Looping parameter
      FIRST_PERIOD,                                               &
              ! First period of precipitation
      LAST_PERIOD,                                                &
              ! Last period of precipitation
      EXTRA_PER_INT
              ! WORK Integer value of EXTRA_PER_REAL



!  Recalculate the total precipitation for the day (in units of mm) 
    PREC_TOT = 0.0
    DO I = 1,STEP_DAY
      PREC_TOT = PREC_TOT + (1.0 / FLOAT(STEP_DAY)) * PREC_LOC(I)
    ENDDO

! Now limit the possible precipitation rate
    DO I = 1,STEP_DAY
      PREC_LOC(I) = MIN(PREC_LOC(I), MAX_PRECIP_RATE)
    ENDDO

! Calculate the new total precipitation for the day
    PREC_TOT_ADJ = 0.0
    DO I = 1,STEP_DAY
      PREC_TOT_ADJ = PREC_TOT_ADJ                                 &
                   + (1.0 / FLOAT(STEP_DAY)) * PREC_LOC(I)
    ENDDO

! Now scatter remaining amount across other periods.
! First revisit the hours where the precipitation may be placed.
! This includes initially before the rainfall event, and then
! after. Start by finding the LAST period. 
    DO I = 1,STEP_DAY 
      IF(N_EVENT_LOCAL(I) == 1) LAST_PERIOD = I
    ENDDO

    FIRST_PERIOD = LAST_PERIOD - N_TALLY + 1

    PREC_CHANGE = PREC_TOT - PREC_TOT_ADJ 

! Calculate the number of extra periods required. Ideally this
! number is less than STEP_DAY-N_TALLY

! Check that only looking at case when PREC_CHANGE is non-zero
    IF(PREC_CHANGE > 0.0) THEN

      EXTRA_PER_REAL = (PREC_CHANGE * FLOAT(STEP_DAY))            &
                     / MAX_PRECIP_RATE
      EXTRA_PER_INT = INT(EXTRA_PER_REAL)

! First case where it is not possible to distribute all the rainfall.
      IF(EXTRA_PER_INT >= STEP_DAY - N_TALLY) THEN
        DO I = 1,STEP_DAY 
          IF(N_EVENT_LOCAL(I) == 0) PREC_LOC(I) = MAX_PRECIP_RATE
        ENDDO
      ENDIF

! Now fill in the time periods before and after the storm.
! Periods are calculated and are moving away from FIRST_PERIOD
! and LAST_PERIOD. Options are as follows:
!
! 1) No periods available before storm, hence all distributed
! after the storm.
!
! 2) All may be accommodated before the storm and nothing after
!
! 3) Then case 3, whereby some of the rain must fall after the
! storm period.
!

! Do case (1); FIRST_PERIOD = 1, hence all precipitation is
! after the storm. 
      IF(FIRST_PERIOD == 1) THEN
        DO I = LAST_PERIOD + 1,LAST_PERIOD + EXTRA_PER_INT
          PREC_LOC(I) = MAX_PRECIP_RATE
        ENDDO
        PREC_LOC(LAST_PERIOD + EXTRA_PER_INT + 1) =               &
           MAX_PRECIP_RATE*(EXTRA_PER_REAL - FLOAT(EXTRA_PER_INT))

! Now do the case where all rain may be accommodated before storm.
      ELSE IF((EXTRA_PER_INT + 1) <= FIRST_PERIOD - 1) THEN
        DO I = FIRST_PERIOD - 1,FIRST_PERIOD - EXTRA_PER_INT,-1
          PREC_LOC(I) = MAX_PRECIP_RATE
        ENDDO
        PREC_LOC(FIRST_PERIOD - EXTRA_PER_INT - 1) =              &
           MAX_PRECIP_RATE*(EXTRA_PER_REAL - FLOAT(EXTRA_PER_INT))

! Now do case where some rain falls after the storm, and all periods
! before the sotrm it rains.
      ELSE
        DO I = FIRST_PERIOD - 1,1,-1
          PREC_LOC(I) = MAX_PRECIP_RATE
        ENDDO
   
        DO I = LAST_PERIOD + 1,                                   &
               LAST_PERIOD + EXTRA_PER_INT - (FIRST_PERIOD - 1)
          PREC_LOC(I) = MAX_PRECIP_RATE
        ENDDO

        PREC_LOC(LAST_PERIOD+EXTRA_PER_INT-(FIRST_PERIOD-1)+1) =  &
           MAX_PRECIP_RATE*(EXTRA_PER_REAL - FLOAT(EXTRA_PER_INT))
      ENDIF      !End of distribution options. 
    ENDIF        !End of check whether redistribution is required.

    RETURN

  END SUBROUTINE REDIS
#endif
