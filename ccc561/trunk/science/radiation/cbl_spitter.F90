!******************************************************************************
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
!******************************************************************************

MODULE cbl_spitter_module

!-----------------------------------------------------------------------------
! Description: 
!   Calculates the beam fraction from Spitters et al. 1986, agric. 
!   for meteorol., 38:217-229
!
! This MODULE is USEd in:
!     cbl_init_radiation.F90 (JULES)
! 
! This MODULE contains 1 public Function:
!     spitter
!
! Module specific documentation: https://trac.nci.org.au/trac/cable/wiki/TBC
! Where it fits in the model flow: https://trac.nci.org.au/trac/cable/wiki/TBC
!-----------------------------------------------------------------------------
   
   IMPLICIT NONE

   PUBLIC spitter
   PRIVATE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate beam fraction: See Spitters et al. 1986, agric. for meteorol., 38:217-229
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION spitter(mp, cpi, doy, coszen, fsd) RESULT(fbeam)
implicit none
!re-decl input args
integer :: mp
real :: cpi
integer :: doy(mp)       ! day of year !typecast from integer
real :: coszen(mp)    ! cos(zenith angle of sun)
real :: fsd(mp)       ! short wave down (positive) w/m^2

   REAL, DIMENSION(mp) ::                                                      &
      fbeam,      & ! beam fraction (result)
      tmpr,       & !
      tmpk,       & !
      tmprat        !

   REAL, PARAMETER :: solcon = 1370.0

   fbeam = 0.0
   tmpr = 0.847 + coszen * (1.04 * coszen - 1.61)
   tmpk = (1.47 - tmpr) / 1.66

   WHERE (coszen > 1.0e-10 .AND. fsd > 10.0)
      tmprat = fsd / ( solcon * ( 1.0 + 0.033 * COS( 2. * CPI* ( doy-10.0 )&
               / 365.0 ) ) * coszen )
   ELSEWHERE
     tmprat = 0.0
   END WHERE

   WHERE ( tmprat > 0.22 ) fbeam = 6.4 * ( tmprat - 0.22 )**2

   WHERE ( tmprat > 0.35 ) fbeam = MIN( 1.66 * tmprat - 0.4728, 1.0 )

   WHERE ( tmprat > tmpk ) fbeam = MAX( 1.0 - tmpr, 0.0 )

END FUNCTION spitter


END MODULE cbl_spitter_module
