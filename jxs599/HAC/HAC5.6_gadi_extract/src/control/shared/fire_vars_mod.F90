! *****************************COPYRIGHT****************************************
! (c) Crown copyright, Met Office. All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
! *****************************COPYRIGHT****************************************
!
! Module containing surface fluxes.
!
! Code Description:
!   Language: FORTRAN 90
!
! Code Owner: Please refer to ModuleLeaders.txt
!

MODULE fire_vars


USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Cloud to ground lightning strikes - only required for standalone
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                        &
  flash_rate(:),                                                              &
    ! The lightning flash rate (flashes/km2), Cloud to Ground
  pop_den(:),                                                                 &
    ! The population density (ppl/km2)
  flammability_ft(:,:),                                                       &
    ! PFT-specific flammability
  burnt_area(:),                                                              &
    ! Gridbox mean burnt area fraction (/s)
  burnt_area_ft(:,:),                                                         &
    ! PFT burnt area fraction (/s)
  emitted_carbon(:),                                                          &
    ! Gridbox mean emitted carbon (kgC/m2/s)
  emitted_carbon_ft(:,:),                                                     &
    ! PFT emitted carbon (kgC/m2/s)
  emitted_carbon_DPM(:),                                                      &
    ! DPM emitted carbon (kgC/m2/s)
  emitted_carbon_RPM(:),                                                      &
    ! RPM emitted carbon (kgC/m2/s)
  fire_em_CO2(:),                                                             &
    ! Gridox mean fire CO2 emission (kg/m2/s)
  fire_em_CO2_ft(:,:),                                                        &
    ! PFT fire CO2 emission (kg/m2/s)
  fire_em_CO2_DPM(:),                                                         &
    ! DPM fire CO2 emission (kg/m2/s)
  fire_em_CO2_RPM(:),                                                         &
    ! RPM fire CO2 emission (kg/m2/s)
  fire_em_CO(:),                                                              &
    ! Gridbox mean fire CO emission (kg/m2/s)
  fire_em_CO_ft(:,:),                                                         &
    ! PFT fire CO emission (kg/m2/s)
  fire_em_CO_DPM(:),                                                          &
    ! DPM fire CO emission (kg/m2/s)
  fire_em_CO_RPM(:),                                                          &
    ! RPM fire CO emission (kg/m2/s)
  fire_em_CH4(:),                                                             &
    ! Gridbox mean fire CH4 emission (kg/m2/s)
  fire_em_CH4_ft(:,:),                                                        &
    ! PFT fire CH4 emission (kg/m2/s)
  fire_em_CH4_DPM(:),                                                         &
    ! DPM fire CH4 emission (kg/m2/s)
  fire_em_CH4_RPM(:),                                                         &
    ! RPM fire CH4 emission (kg/m2/s)
  fire_em_NOx(:),                                                             &
    ! Gridbox mean fire NOx emission (kg/m2/s)
  fire_em_NOx_ft(:,:),                                                        &
    ! PFT fire NOx emission (kg/m2/s)
  fire_em_NOx_DPM(:),                                                         &
    ! DPM fire NOx emission (kg/m2/s)
  fire_em_NOx_RPM(:),                                                         &
    ! RPM fire NOx emission (kg/m2/s)
  fire_em_SO2(:),                                                             &
    ! Gridbox mean fire SO2 emission (kg/m2/s)
  fire_em_SO2_ft(:,:),                                                        &
    ! PFT fire SO2 emission (kg/m2/s)
  fire_em_SO2_DPM(:),                                                         &
    ! DPM fire SO2 emission (kg/m2/s)
  fire_em_SO2_RPM(:),                                                         &
    ! RPM fire SO2 emission (kg/m2/s)
  fire_em_OC(:),                                                              &
    ! Gridbox mean fire Organic Carbon emission (kg/m2/s)
  fire_em_OC_ft(:,:),                                                         &
     ! PFT fire OC emission (kg/m2/s)
  fire_em_OC_DPM(:),                                                          &
     ! DPM fire OC emission (kg/m2/s)
  fire_em_OC_RPM(:),                                                          &
     ! RPM fire OC emission (kg/m2/s)
  fire_em_BC(:),                                                              &
     ! Gridbox mean fire Black Carbon emission (kg/m2/s)
  fire_em_BC_ft(:,:),                                                         &
     ! PFT fire BC emission (kg/m2/s)
  fire_em_BC_DPM(:),                                                          &
     ! DPM fire BC emission (kg/m2/s)
  fire_em_BC_RPM(:)
     ! RPM fire BC emission (kg/m2/s)

END MODULE fire_vars
