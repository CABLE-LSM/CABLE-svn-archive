! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds parameters for each crop plant functional type

MODULE cropparm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Module containing parameters for each crop PFT
!
! Code Owner: Please refer to ModuleLeaders.txt
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

! Public module variables
REAL, ALLOCATABLE ::                                                          &
  t_bse(:),                                                                   &
      ! Base temp for development (K).
  t_opt(:),                                                                   &
      ! Optimum temp for development (K).
  t_max(:),                                                                   &
      ! Maximum temp for development (K).
  tt_emr(:),                                                                  &
      ! Thermal time for emergence (degree days).
  crit_pp(:),                                                                 &
      ! Critical daylength for photoperiod sensitivity (hours).
  pp_sens(:),                                                                 &
      ! Sensitivity to daylength (hours-1).
  rt_dir(:),                                                                  &
      ! Alpha for root growth direction
  alpha1(:),                                                                  &
  alpha2(:),                                                                  &
  alpha3(:),                                                                  &
  beta1(:),                                                                   &
  beta2(:),                                                                   &
  beta3(:),                                                                   &
      ! Coefficients to calculate partition coefficients
  r_gamma(:),                                                                 &
  delta(:),                                                                   &
      ! Coefficients for sla calculation (m2 kg-1).
  remob(:),                                                                   &
      ! Remobilisation factor
  cfrac_s(:),                                                                 &
      ! Carbon fraction of stems
  cfrac_r(:),                                                                 &
      ! Carbon fraction of roots
  cfrac_l(:),                                                                 &
      ! Carbon fraction of leaves
  allo1(:),                                                                   &
  allo2(:),                                                                   &
      ! Allometric coefficients for stemc <-> canht
  mu(:),                                                                      &
      ! Coefficient for senescence calculation
  nu(:),                                                                      &
      ! Coefficient for senescence calculation
  yield_frac(:),                                                              &
      ! Fraction of the harv carbon pool converted to yield carbon
  initial_carbon(:),                                                          &
      ! Carbon in crops at DVI=initial_c_dvi (kgC/m2).
  initial_c_dvi(:),                                                           &
      ! DVI at which total crop carbon is set to initial_carbon. Should be 
      ! emergence (0.0) or shortly after.
  sen_dvi(:),                                                                 &
      ! DVI at which leaf senescence begins. 
  t_mort(:)
      ! Soil temperature (second level) at which to kill crop if dvi>1 (K).
        
END MODULE cropparm
