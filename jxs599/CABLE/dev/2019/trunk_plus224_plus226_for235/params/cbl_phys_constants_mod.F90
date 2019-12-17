MODULE cable_phys_constants_mod

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   CABLE physical constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in 
!-----------------------------------------------------------------------------

  REAL, PARAMETER :: tfrz   = 273.16     ! Temp (K) corresp. to 0 C
  REAL, PARAMETER :: sboltz = 5.67e-8    ! Stefan-Boltz. constant (W/m2/K4)
  REAL, PARAMETER :: emsoil = 1.0        ! soil emissivity
  REAL, PARAMETER :: emleaf = 1.0        ! leaf emissivity
  REAL, PARAMETER :: capp   = 1004.64    ! air spec. heat (J/kg/K)

END MODULE cable_phys_constants_mod
