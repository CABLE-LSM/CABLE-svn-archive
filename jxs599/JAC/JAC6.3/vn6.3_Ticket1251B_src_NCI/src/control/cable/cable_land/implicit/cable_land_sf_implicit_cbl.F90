MODULE cable_land_sf_implicit_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
                  ModuleName='CABLE_LAND_SF_IMPLICIT_MOD'

CONTAINS
!  SUBROUTINE CABLE_LAND_SF_IMPLICIT --------------------------------
!
!  Purpose: 
!
!--------------------------------------------------------------------
!    Arguments :-
SUBROUTINE cable_land_sf_implicit ( progs, work, pars )
! In general CABLE utilizes a required subset of tbe JULES types, however;
USE progs_cbl_vars_mod, ONLY: progs_cbl_vars_type ! CABLE requires extra progs
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep
USE params_io_mod_cbl,  ONLY: params_io_type      ! and veg/soil parameters

IMPLICIT NONE

!CABLE TYPES containing field data
TYPE(progs_cbl_vars_type), INTENT(IN OUT) :: progs
TYPE(work_vars_type), INTENT(OUT)        :: work
TYPE(params_io_type), INTENT(IN)         :: pars

!fields returned @ hydrology (surf_couple_extra) level of JULES
!computed in CABLE @ implicit (surf_couple_implicit) level
work%snow_tile     = 0.0 ! will be UNPACKed from CABLE field
work%lying_snow    = 0.0 ! will be UNPACKed from CABLE field
work%tot_tfall     = 0.0 ! will be UNPACKed from CABLE field
work%surf_roff     = 0.0 ! will be UNPACKed from CABLE field
work%sub_surf_roff = 0.0 ! will be UNPACKed from CABLE field

WRITE(6,*) "Currently CABLE implicit@6.3 is not implemented"

RETURN
END SUBROUTINE cable_land_sf_implicit
END MODULE cable_land_sf_implicit_mod
