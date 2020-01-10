#if !defined(UM_JULES)
! Module containing variables for coastal tiling

MODULE coastal

IMPLICIT NONE

REAL, DIMENSION(:),   ALLOCATABLE :: fland
  ! Land fraction on land tiles
  ! For offline JULES, this might be better considered as a "0/1 flag to
  !indicate which gridboxes are to be modelled", since we can use it to select
  !land points from a larger set of land points.
  ! As long as JULES is not simulating sea points, FLAND and FLANDG
  !(which is what's input) should be 0 or 1 (not intermediate values).
REAL, DIMENSION(:,:), ALLOCATABLE :: flandg
  ! Land fraction on all tiles.
  !   Divided by 2SQRT(2) on land points only (m)
REAL, DIMENSION(:,:), ALLOCATABLE :: tstar_land_ij
  ! Land mean surface temperature (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: tstar_sea_ij
  ! Open sea surface temperature (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: tstar_sice_ij
  ! Sea-ice surface temperature (K)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: tstar_sice_sicat
  ! Category sea-ice surface temperature (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: tstar_ssi_ij
  ! mean sea surface temperature (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: taux_land_ij
  ! W'ly component of sfc wind stress (N/sq m).
  !   (On U-grid with first and last rows undefined or,
  !    at present, set to missing data)
REAL, DIMENSION(:,:), ALLOCATABLE :: taux_ssi_ij
  ! W'ly component of sfc wind stress (N/sq m).
  !   (On U-grid with first and last rows undefined or,
  !    at present, set to missing data)
REAL, DIMENSION(:,:), ALLOCATABLE :: tauy_land_ij
  ! S'ly component of sfc wind stress (N/sq m).
  !   On V-grid; comments as per TAUX
REAL, DIMENSION(:,:), ALLOCATABLE :: tauy_ssi_ij
  ! S'ly component of sfc wind stress (N/sq m).
  !   On V-grid; comments as per TAUX
REAL, DIMENSION(:,:), ALLOCATABLE :: vshr_land_ij
  ! Magnitude of surface-to-lowest atm level wind shear (m per s)
REAL, DIMENSION(:,:), ALLOCATABLE :: vshr_ssi_ij
  ! Magnitude of surface-to-lowest atm level wind shear (m per s)
REAL, DIMENSION(:,:), ALLOCATABLE :: surf_ht_flux_land_ij
  ! Net downward heat flux at surface over land fraction of gridbox (W/m2)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: surf_ht_flux_sice_sicat
  ! Net downward heat flux at surface over sea-ice fraction of gridbox (W/m2)

REAL, ALLOCATABLE ::                                                          &
  taux_land_star(:,:),                                                        &
  tauy_land_star(:,:),                                                        &
  taux_ssi_star(:,:),                                                         &
  tauy_ssi_star(:,:)
   

END MODULE coastal
#endif
