! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting Parameters for each plant functional type, for TRIFFID.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE trif

IMPLICIT NONE

INTEGER, ALLOCATABLE ::                                                       &
 crop(:)          !  Flag for crop types: 1 for crop type,
                  !  0 for non-crop.

REAL, ALLOCATABLE ::                                                          &
 g_area(:)                                                                    &
                  !  Disturbance rate (/360days).
,g_grow(:)                                                                    &
                  !  Rate of leaf growth (/360days)
,g_root(:)                                                                    &
                  !  Turnover rate for root biomass (/360days).
,g_wood(:)                                                                    &
                  !  Turnover rate for woody biomass (/360days).
,lai_max(:)                                                                   &
                  !  Maximum projected LAI.
,lai_min(:)                                                                   &
                  !  Minimum projected LAI
,alloc_fast(:)                                                                &
                  ! Fraction of landuse carbon allocated fast
                  ! product pool
,alloc_med(:)                                                                 &
                  ! Fraction of landuse carbon allocated medium
                  ! product pool
,alloc_slow(:)                                                                &
                  ! Fraction of landuse carbon allocated slow
                  ! product pool
,dpm_rpm_ratio(:)                                                             &
                  !  Ratio of each PFTs litter allocated to the
                  !  DPM soil carbon pool versus the RPM soil
                  !  carbon pool
,retran_l(:)                                                                  &
                  ! Leaf N retranslocation factor
,retran_r(:)
                  ! Root N retranslocation factor

END MODULE trif
