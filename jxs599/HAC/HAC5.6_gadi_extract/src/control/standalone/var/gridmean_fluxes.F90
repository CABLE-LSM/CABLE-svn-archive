#if !defined(UM_JULES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing grid-box mean surface fluxes would be passed 
! back to an atmopheric model
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! Code Owner: Please refer to ModuleLeaders.txt and UM file CodeOwners.txt
! This file belongs in section: Land
!

MODULE gridmean_fluxes

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module is required only in standalone configuration
!-----------------------------------------------------------------------------
REAL, ALLOCATABLE :: fqw_1_ij(:,:)
!   Moisture flux between layers (kg per square metre per sec)
!   FQW(,1) is total water flux from surface, 'E'
REAL, ALLOCATABLE :: ftl_1_ij(:,:)
!   FTL(,K) contains net turbulent sensible heat flux into layer K from below
!   so FTL(,1) is the surface sensible heat, H.(W/m2)
REAL, ALLOCATABLE :: taux_1_ij(:,:)
!   W'ly component of surface wind stress (N/sq m)
REAL, ALLOCATABLE :: tauy_1_ij(:,:)
!   S'ly component of surface wind stress (N/sq m)
!   On V-grid; comments as per TAUX

END MODULE gridmean_fluxes

#endif
