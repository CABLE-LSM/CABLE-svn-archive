#if !defined(UM_JULES)
! Module containing variables for coastal tiling

MODULE coastal

  IMPLICIT NONE

  REAL, DIMENSION(:),   ALLOCATABLE :: FLAND
!                                        ! Land fraction on land tiles
!        For offline JULES, this might be better considered as a "0/1 flag to indicate which
!        gridboxes are to be modelled", since we can use it to select land points from
!        a larger set of land points.
!        As long as JULES is not simulating sea points, FLAND and FLANDG (which is what's input)
!        should be 0 or 1 (not intermediate values).
  REAL, DIMENSION(:,:), ALLOCATABLE :: FLANDG
!                                    ! Land fraction on all tiles. 
!                                    !   Divided by 2SQRT(2) on land points only (m)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TSTAR_LAND
!                                    ! Land mean surface temperature (K)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TSTAR_SEA
!                                    ! Open sea surface temperature (K)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TSTAR_SICE
!                                    ! Sea-ice surface temperature (K)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: TSTAR_SICE_NCAT
!                                    ! Category sea-ice surface temperature (K)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TSTAR_SSI
!                                    ! mean sea surface temperature (K)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TAUX_LAND
!                                    ! W'ly component of sfc wind stress (N/sq m).
!                                    !   (On U-grid with first and last rows undefined or,
!                                    !    at present, set to missing data)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TAUX_SSI
!                                    ! W'ly component of sfc wind stress (N/sq m).
!                                    !   (On U-grid with first and last rows undefined or,
!                                    !    at present, set to missing data)
  REAL, DIMENSION(:,:), ALLOCATABLE :: TAUY_LAND
!                                    ! S'ly component of sfc wind stress (N/sq m).
!                                    !   On V-grid; comments as per TAUX
  REAL, DIMENSION(:,:), ALLOCATABLE :: TAUY_SSI
!                                    ! S'ly component of sfc wind stress (N/sq m).
!                                    !   On V-grid; comments as per TAUX
  REAL, DIMENSION(:,:), ALLOCATABLE :: VSHR_LAND
!                                    ! Magnitude of surface-to-lowest atm level wind shear (m per s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: VSHR_SSI
!                                    ! Magnitude of surface-to-lowest atm level wind shear (m per s)
  REAL, DIMENSION(:,:), ALLOCATABLE :: SURF_HT_FLUX_LAND
!                                    ! Net downward heat flux at surface over land fraction of gridbox (W/m2)
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: SURF_HT_FLUX_SICE
!                                    ! Net downward heat flux at surface over sea-ice fraction of gridbox (W/m2)

END MODULE coastal
#endif
