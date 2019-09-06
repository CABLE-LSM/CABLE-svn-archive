#if !defined(UM_JULES)

MODULE init_cable_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE init_cable_grid()

USE cable_types_mod,          ONLY: l_tile_pts
USE ancil_info,               ONLY: frac_surft, land_pts
USE jules_surface_types_mod,  ONLY: ntype

IMPLICIT NONE

!------------------------------------------------------------------------------
! Description:
!   Initialises the JULES/CABLE grid array, which aligns JULES grid points
!   with CABLE land points
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!------------------------------------------------------------------------------

INTEGER :: i, j

ALLOCATE(l_tile_pts(land_pts, ntype))

l_tile_pts(:,:) = .FALSE.

DO j = 1, ntype
  DO i = 1, land_pts
    IF ( frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO

RETURN

END SUBROUTINE init_cable_grid


SUBROUTINE init_cable_veg()

USE cable_types_mod,         ONLY: veg, vegin, mp, l_tile_pts
USE ancil_info,              ONLY: surft_index

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   init_cables veg parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

INTEGER :: h

veg%iveg   = PACK(surft_index, l_tile_pts)

! Prescribe parameters for current gridcell based on veg/soil type
! (which may have loaded from default value file or met file):
DO h = 1, mp          ! over each patch in current grid
  veg%taul(h,1)   = vegin%taul(1,veg%iveg(h))
  veg%taul(h,2)   = vegin%taul(2,veg%iveg(h))
  veg%refl(h,1)   = vegin%refl(1,veg%iveg(h))
  veg%refl(h,2)   = vegin%refl(2,veg%iveg(h))
  veg%hc(h)       = vegin%hc(veg%iveg(h))
END DO ! over each veg patch in land point

END SUBROUTINE init_cable_veg


SUBROUTINE init_cable_soil()

USE cable_types_mod,          ONLY: ssnow, soil, veg, albsoil, mp,            &
                                    perm_ice_veg, perm_ice_soil, non_ice_soil
USE p_s_parms,                ONLY: albsoil_soilt
USE cable_pack_mod,           ONLY: cable_pack_lp
USE ancil_info,               ONLY: land_pts

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises soil parameters using values read from namelist
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

LOGICAL, PARAMETER :: skip =.TRUE.
REAL :: dummy(10)
INTEGER :: i

!--- soil%isoilm defines soiltype.


! set soil type for permanent ice based on where permanent ice
! located in vegetation map
soil%isoilm(:) = non_ice_soil
DO i = 1, mp
  IF (veg%iveg(i) == perm_ice_veg)                                            &
    soil%isoilm(:) = perm_ice_soil
END DO

! set for CABLE which has only one tile per land point
ALLOCATE(albsoil(land_pts))
albsoil(:) = albsoil_soilt(:,1)

! set CABLE-var soil%albsoil from UM/JULES var albsoil
CALL cable_pack_lp( albsoil, dummy, soil%albsoil(:,1), soil%isoilm, skip )

END SUBROUTINE init_cable_soil


SUBROUTINE init_cable_met()

USE cable_types_mod,            ONLY: met, rad_bands
USE cable_other_constants_mod,  ONLY: n_sw_bands
USE ancil_info,                 ONLY: row_length, rows
USE cable_pack_mod,             ONLY: cable_pack_rr, cable_pack_met
USE forcing,                    ONLY: tl_1_ij

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Repacks UM/JULES forcing vars of dimension (row_length, rows) into a
!   single vector of active tiles.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

ALLOCATE(rad_bands%sw_down_dir(row_length, rows))
ALLOCATE(rad_bands%sw_down_dif(row_length, rows))
ALLOCATE(rad_bands%sw_down_vis(row_length, rows))
ALLOCATE(rad_bands%sw_down_nir(row_length, rows))
ALLOCATE(rad_bands%fbeam(row_length, rows, n_sw_bands))

CALL cable_pack_met()

CALL cable_pack_rr(tl_1_ij, met%tk)

met%coszen = MAX(met%coszen, EPSILON(0.0))

END SUBROUTINE init_cable_met


SUBROUTINE init_cable_rad()

USE cable_types_mod,            ONLY: met, rad, mp
USE cable_other_constants_mod,  ONLY: rad_thresh
USE cable_pack_mod,             ONLY: cable_pack_rad

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises and repacks radiation variables
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!
!-----------------------------------------------------------------------------

INTEGER :: i

CALL cable_pack_rad()

DO i = 1, mp
  IF (met%coszen(i) < rad_thresh) THEN
    rad%fbeam(i,1) = 0.0
    rad%fbeam(i,2) = 0.0
    rad%fbeam(i,3) = 0.0
  END IF
END DO

END SUBROUTINE init_cable_rad


SUBROUTINE init_cable_canopy()

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Repacks UM/JULES canopy variables of dimension (row_length, rows) into a
!   single vector of active tiles.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Added for the sake of setting up directory structure.
! Code will be added as required.

END SUBROUTINE init_cable_canopy


SUBROUTINE init_cable_soilsnow()

USE cable_types_mod,           ONLY: ssnow, l_tile_pts
USE cable_other_constants_mod, ONLY: msn, max_snow_depth, init_snow_rho1l
USE ancil_info,                ONLY: land_pts, nsoilt, soilt_index, soilt_pts
USE jules_surface_types_mod,   ONLY: ntype
USE prognostics,               ONLY: tsnow_surft, t_soil_soilt, snow_surft

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Initialises CABLE soilsnow variables
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

REAL :: snow_tile(land_pts, ntype)
INTEGER :: i,j,t,p


snow_tile = MIN(max_snow_depth, snow_surft)

ssnow%snowd(:) = 0.0
ssnow%ssdnn(:) = init_snow_rho1l
ssnow%isflag(:) = 0

DO j = 1, msn
  ssnow%tggsn(:,j) = PACK(tsnow_surft(:,:,j), l_tile_pts)
END DO

DO t = 1,nsoilt
  DO i = 1,soilt_pts(t)
    p = soilt_index(i,t)
    ssnow%tgg(p,:) = t_soil_soilt(i,t,:)
  END DO
END DO

ssnow%snage(:) = 0.0

! This is normally the snow depth from the previous time step
ssnow%osnowd(:)  = PACK(snow_tile, l_tile_pts)

END SUBROUTINE init_cable_soilsnow

END MODULE init_cable_mod

#endif
