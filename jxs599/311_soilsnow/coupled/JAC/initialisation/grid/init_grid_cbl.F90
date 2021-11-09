MODULE init_grid_mod_cbl

PUBLIC init_grid_cbl
PUBLIC  mp
PUBLIC  L_tile_pts

INTEGER :: mp ! CABLE uses a 1-D vector of active cells

! mask ID-s active cells in (land_pts,ntiles) format, used in mapping to 1-D 
! vector of length mp format used by CABLE 
LOGICAL, ALLOCATABLE :: L_tile_pts(:,:) 

CONTAINS

SUBROUTINE init_grid_cbl() 

USE init_model_grid_arrays_mod_cbl, ONLY: init_model_grid_arrays_cbl
USE vegin_pars_mod_cbl,             ONLY: vegin
USE soilin_pars_mod_cbl,            ONLY: soilin
USE cable_fields_mod,               ONLY: veg     => veg_pars_data_cbl,       &
                                          veg_ptr => veg_pars_cbl,            &
                                          soil    => soil_pars_data_cbl,      &
                                          soil_ptr => soil_pars_cbl
USE ancil_info,                     ONLY: nsurft, land_pts,                   &
                                          surft_pts
USE jules_fields_mod,               ONLY: ainfo

IMPLICIT NONE

!local vars
INTEGER :: i, j

!can we call tilepts from here to get surft_pts? and hence mp
!this is done here B/C and only IF after call init_ic  
! Determine the number of active tiles
mp = SUM(surft_pts) 

! Determine active tiles map
IF ( .NOT. ALLOCATED(l_tile_pts)) ALLOCATE( l_tile_pts(land_pts, nsurft) )

l_tile_pts(:,:) = .FALSE.

DO j = 1, nsurft
  DO i = 1, land_pts
    IF ( ainfo%frac_surft(i,j)  >   0.0 ) THEN
      l_tile_pts(i,j) = .TRUE.
    END IF
  END DO
END DO 

CALL init_model_grid_arrays_cbl( mp, nsurft, land_pts, surft_pts,             &
                         ainfo%frac_surft, l_tile_pts, veg, veg_ptr,          &
                         vegin, soil, soil_ptr, soilin )

RETURN
END SUBROUTINE init_grid_cbl

END MODULE init_grid_mod_cbl
