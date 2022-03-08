MODULE cable_land_sf_explicit_mod

CONTAINS

SUBROUTINE cable_land_sf_explicit( progs, work, pars )

! In general CABLE utilizes a required subset of tbe JULES types, however;
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep
USE params_io_mod_cbl,  ONLY: params_io_type      ! and veg/soil parameters

implicit NONE

!CABLE TYPES containing field data
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs
TYPE(work_vars_type), INTENT(OUT)        :: work
TYPE(params_io_type), INTENT(IN)         :: pars

DO n = 1,ntype
  DO l = 1, land_pts
    tile_frac(l,n) = frac(l,n)
  END DO
END DO

WRITE(6,*) "Currently CABLE explicit@6.3 is not implemented"
RETURN

END SUBROUTINE cable_land_sf_explicit

END MODULE cable_land_sf_explicit_mod


