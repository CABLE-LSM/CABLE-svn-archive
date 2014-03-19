#if !defined(UM_JULES)
  MODULE imogen_time

    IMPLICIT NONE

!-----------------------------------------------------------------
! Define basic timestepping variables
!-----------------------------------------------------------------
    INTEGER ::                                                    &
      YEAR1=1860,                                                 &
                   ! First year of the numerical experiment
      IYSTART=1860,                                               &
                   ! Start year for this submitted run
      IYEND=2100,                                                 &
                   ! Stop year of the ENTIRE run
      IYEAR        ! Year we are at in the run

    INTEGER ::                                                    &
      MM,                                                         &
                ! Number of months in year
      MD,                                                         &
                ! Number of days in (GCM) month
      STEP_DAY,                                                   &
                ! Number of daily timesteps
                ! This is derived from the timestep
      NSDMAX       
                ! Maximum number of possible subdaily increments
    
    PARAMETER(MM=12)
    PARAMETER(MD=30)   !At present hardwired as 30 days per month (as
                       !in GCM)
                       !l_360 is forced to be true at the control level
    PARAMETER(NSDMAX=24)

  END MODULE imogen_time
#endif
