!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: subroutines for calculating carbon, nitrogen, phosphorus cycle
!          including plant growth
!
! Called from: biogeochem (mostly) or casa_xnp
!
! Contact: Yingping.Wang@csiro.au
!
! History: Developed by Yingping Wang (Wang et al., BG, 2011)
!          Current version uses fixed phenology.
!
! Sep 2015: option of climate-driven phenology (V. Haverd)
!           search for cable_user%PHENOLOGY_SWITCH (Ticket #110)
! May 2016: option of acclimation of auttrophic respiration (V. Haverd)
!            search for cable_user%CALL_climate (Ticket#110)
!         : fixes to prevent carbon and nitrogen pools from going negative
!           search for Ticket#108 (V.Haverd)
!         : alternative functional form of vcmax, called when cable_user%vcmax=='Walker2014'
!           (V.Haverd)
!         : alternative allocation switch integer: LALLOC=3. (V.Haverd)
!           leaf:wood allocation set to maintain LA:SA ratio
!           below target value (requires casaflux%sapwood_area 
!           inherited from POP demography module. (Ticket#61)
! ==============================================================================
!
! This module contains the following subroutines:
!   casa_xnp
!   casa_allocation
!   casa_rplant
!   casa_xrateplant,        casa_xratesoil
!   casa_coeffplant,        casa_coeffsoil
!   casa_delplant,          casa_delsoil
!   avgsoil
!   casa_xkN
!   casa_nuptake,           casa_puptake
!   casa_Nrequire,          casa_Prequire
!   casa_cnpcycle
!   casa_poolzero
!   casa_cnpbal
!   casa_ndummy
!   phenology

MODULE casa_cnp_module
USE cable_def_types_mod
USE casadimension
USE casaparm
USE casavariable
USE phenvariable
USE cable_common_module, only: cable_user ! Custom soil respiration: Ticket #42
IMPLICIT NONE
  REAL(r_2), PARAMETER :: zero = 0.0_r_2
  REAL(r_2), PARAMETER :: one  = 1.0_r_2
CONTAINS



END MODULE casa_cnp_module
