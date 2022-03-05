!    SUBROUTINE HYDROL_CBL--------------------------------------------------------

! Description:

MODULE hydrol_mod_cbl
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYDROL_MOD_CBL'

CONTAINS

SUBROUTINE hydrol_cbl ( snow_surft, snow_mass_gb, tot_tfall_gb,                &
                        surf_roff_gb, sub_surf_roff_gb,                        &
                        land_pts, nsurft,                                      &
                        work_cbl ) 

! In general CABLE utilizes a required subset of tbe JULES types, however;
USE work_vars_mod_cbl,  ONLY: work_vars_type      ! and some kept thru timestep

IMPLICIT NONE

!CABLE unpacked to work_cbl%  @ "implicit". pass to JULES fields here
REAL, INTENT(OUT) ::                                                           &
  snow_surft(land_pts,nsurft),                                                 &
  snow_mass_gb(land_pts),                                                      &
  tot_tfall_gb(land_pts),                                                      &
    ! Total throughfall (kg/m2/s).
  sub_surf_roff_gb(land_pts),                                                  &
    ! Sub-surface runoff (kg/m2/s).
  surf_roff_gb(land_pts)
    ! Surface runoff (kg/m2/s).

INTEGER, INTENT(IN) :: land_pts, nsurft

!CABLE TYPES containing field data (IN OUT)
TYPE(work_vars_type), INTENT(IN OUT)       :: work_cbl

CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDROL_CBL'

! End of header---------------------------------------------------------------

!fields returned @ hydrology (surf_couple_extra) level of JULES
!computed in CABLE @ implicit (surf_couple_implicit) level
snow_surft       = work_cbl%snow_tile
snow_mass_gb     = work_cbl%lying_snow
tot_tfall_gb     = work_cbl%tot_tfall
surf_roff_gb     = work_cbl%surf_roff
sub_surf_roff_gb = work_cbl%sub_surf_roff

WRITE(6,*) "Currently CABLE extra@6.3 is not implemented"

RETURN
END SUBROUTINE hydrol_cbl
END MODULE hydrol_mod_cbl
