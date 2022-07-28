MODULE cable_surface_types_mod
! Index identifiers of types in veg%iveg, soil%isoilm
! JAC has file of same name, however there the Surface type (both PFT&nveg) 
! indices are all read from a namelist. SoilType info is complicated, and 
! so had-wired there (as it is here)
IMPLICIT NONE
!#PFT/vegertated types 
INTEGER, PARAMETER :: evergreen_needleleaf = 0
INTEGER, PARAMETER :: evergreen_broadleaf = 0
INTEGER, PARAMETER :: deciduous_needleleaf = 0
INTEGER, PARAMETER :: deciduous_broadleaf = 0
INTEGER, PARAMETER :: shrub_cable = 0
INTEGER, PARAMETER :: c3_grassland = 0
INTEGER, PARAMETER :: c4_grassland = 0
INTEGER, PARAMETER :: tundra = 0
INTEGER, PARAMETER :: c3_cropland = 0
INTEGER, PARAMETER :: c4_cropland = 0
INTEGER, PARAMETER :: wetland = 0
INTEGER, PARAMETER :: empty1 = 0
INTEGER, PARAMETER :: empty2 = 0
!#nveg types
INTEGER, PARAMETER :: barren_cable = 0
INTEGER, PARAMETER :: urban_cable = 0
INTEGER, PARAMETER :: lakes_cable = 16 
INTEGER, PARAMETER :: ICE_cable = 17

!#soil types
INTEGER, PARAMETER :: ICE_SoilType = 9

END MODULE cable_surface_types_mod
