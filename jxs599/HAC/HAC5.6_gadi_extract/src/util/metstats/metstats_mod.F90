! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use
! and distribution under the JULES collaboration agreement, subject
! to the terms and conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************

MODULE metstats_mod

  !No imports

IMPLICIT NONE
!
! Description:
!   Module defining a set of TYPE variables to store data for various
!   meteorological statistics.
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.

  !Switch for whether we're acutally using the metstats
LOGICAL :: l_metstats = .FALSE.

!Template TYPE for cumulative statistics (eg max, min, ave)
TYPE stat_cum
  REAL :: run !Running 'total' since the start of sampling period
  REAL :: fin !Value logged at the end of last sampling period
END TYPE stat_cum

!Template TYPE for stats that are simply the value at a defined timestamp
!(eg at 12h local)
TYPE stat_pnt
  REAL :: fin !Value logged at the end of last sampling period
END TYPE stat_pnt

!----------------------------------------------------------------------------
!Data structure to hold all the metstats. Prognostic.

!Times indicate when the %fin fields are populated and refer to local time
!(as opposed to the model timestamp which is usually GMT+0)
!%fin fields are overwritten at the end of their sampling periods.
!At present, all metstats are daily stats.

!Examples:
!temp_max_00h is the maximum temperature in the 24h to midnight local.
!rhum_pnt_12h is the relative humidity at 1200 local.
TYPE metstats_prog_struct
  TYPE (stat_cum) :: temp_max_00h !Temperature   Units K
  TYPE (stat_cum) :: temp_ave_00h
  TYPE (stat_pnt) :: temp_pnt_12h

  TYPE (stat_cum) :: prec_tot_00h !Precip        Units mm
  TYPE (stat_cum) :: prec_tot_12h

  TYPE (stat_cum) :: rhum_min_00h !Rel hums      Units %
  TYPE (stat_pnt) :: rhum_pnt_12h

  TYPE (stat_cum) :: dewp_ave_00h !Dewpoint      Units K

  TYPE (stat_cum) :: wind_ave_00h !Wind          Units m/s
  TYPE (stat_pnt) :: wind_pnt_12h

  !Difference between model timestamp and local time
  REAL            :: lon_time_diff
END TYPE metstats_prog_struct


!Contains flags for each stat to control whether they are calculated or not
TYPE metstats_flag_struct
  LOGICAL         :: temp_max_00h = .FALSE.,                                  &
                     temp_ave_00h = .FALSE.,                                  &
                     temp_pnt_12h = .FALSE.,                                  &
                     prec_tot_00h = .FALSE.,                                  &
                     prec_tot_12h = .FALSE.,                                  &
                     rhum_min_00h = .FALSE.,                                  &
                     rhum_pnt_12h = .FALSE.,                                  &
                     dewp_ave_00h = .FALSE.,                                  &
                     wind_ave_00h = .FALSE.,                                  &
                     wind_pnt_12h = .FALSE.
END TYPE metstats_flag_struct


!Contains initial values to which prognostics should be set
TYPE metstats_init_struct
  REAL ::            temp_max_00h =   0.0, & !(K)  ie min physically possible
                     temp_ave_00h =   0.0, & !(K)
                     temp_pnt_12h =   0.0, & !(K)
                     prec_tot_00h =   0.0, & !(mm)
                     prec_tot_12h =   0.0, & !(mm)
                     rhum_min_00h = 100.0, & !(%)  ie max physically possible
                     rhum_pnt_12h =   0.0, & !(%)
                     dewp_ave_00h =   0.0, & !(K)
                     wind_ave_00h =   0.0, & !(m/s)
                     wind_pnt_12h =   0.0    !(m/s)
END TYPE metstats_init_struct

!Declare instances for running the module
TYPE (metstats_prog_struct) , ALLOCATABLE :: metstats_prog(:)
TYPE (metstats_flag_struct) , SAVE        :: metstats_flag
TYPE (metstats_init_struct) , SAVE        :: metstats_inis

END MODULE metstats_mod
