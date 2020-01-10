! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holding parameter arrays for non-vegetation surface types.


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE nvegparm

IMPLICIT NONE

REAL, ALLOCATABLE ::                                                          &
 albsnc_nvg(:)                                                                &
                  ! Snow-covered albedo.
,albsnf_nvgu(:)                                                               &
                  ! Max Snow-free albedo, when scaled to obs
,albsnf_nvg(:)                                                                &
                  ! Snow-free albedo.
,albsnf_nvgl(:)                                                               &
                  ! Min Snow-free albedo, when scaled to obs
,catch_nvg(:)                                                                 &
                  ! Canopy capacity for water (kg/m2).
,gs_nvg(:)                                                                    &
                  ! Surface conductance (m/s).
,infil_nvg(:)                                                                 &
                  ! Infiltration enhancement factor.
,z0_nvg(:)                                                                    &
                  ! Roughness length (m).
,ch_nvg(:)                                                                    &
                  ! "Canopy" heat capacity (J/K/m2)
,vf_nvg(:)                                                                    &
                  ! Fractional "canopy" coverage
,emis_nvg(:)

END MODULE nvegparm
