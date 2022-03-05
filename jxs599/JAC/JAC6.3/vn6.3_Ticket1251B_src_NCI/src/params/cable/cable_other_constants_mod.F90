MODULE cable_other_constants_mod

USE grid_constants_mod_cbl, ONLY: nrb

IMPLICIT NONE

PUBLIC

!-----------------------------------------------------------------------------
! Description:
!   Other CABLE constants
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

REAL, PARAMETER :: gauss_w(nrb)=[0.308,0.514,0.178 ] ! Gaussian integ. weights
REAL, PARAMETER :: rad_thresh = 0.001
                        ! minimum zenithal angle for downward SW radiation
REAL, PARAMETER :: lai_thresh = 0.001
                        ! threshold for minimum significant LAI

! minimum (cosine)zenith angle of sun signalling sunrise
REAL, PARAMETER :: coszen_tols = 1.0e-4

REAL, PARAMETER :: z0surf_min = 1.0e-7 ! min. roughness of bare soil surface
!H!REAL, PARAMETER :: z0snow_min = 1.e-7 ! min. roughness of bare snow surface

END MODULE cable_other_constants_mod
