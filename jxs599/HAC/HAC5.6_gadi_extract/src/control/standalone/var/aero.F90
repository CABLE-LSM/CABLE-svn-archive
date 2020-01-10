#if !defined(UM_JULES)
! Module containing variables required for aerosol calculations

MODULE aero

IMPLICIT NONE

REAL                              :: co2_mmr = 5.24100e-4
                                  ! CO2 Mass Mixing Ratio
REAL, DIMENSION(:,:), ALLOCATABLE :: co2_3d_ij
                                  ! 3D CO2 field if required
REAL, DIMENSION(:,:), ALLOCATABLE :: rho_cd_modv1_ij
                                  ! Surface air density * drag coef *
                                  ! mod(v1 - v0) before interp
REAL, DIMENSION(:,:), ALLOCATABLE :: rho_aresist_ij
                                  ! RHOSTAR*CD_STD*VSHR for CLASSIC aerosol
                                  ! scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: aresist_ij
                                  ! 1/(CD_STD*VSHR) for CLASSIC aerosol
                                   ! scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: resist_b_ij
                                  ! (1/CH-1/(CD_STD)/VSHR for CLASSIC aerosol
                                  ! scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: rho_aresist_surft
                                  ! RHOSTAR*CD_STD*VSHR on land tiles for
                                  !CLASSIC aerosol scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: aresist_surft
                                  ! 1/(CD_STD*VSHR) on land tiles for CLASSIC
                                  ! aerosol scheme
REAL, DIMENSION(:,:), ALLOCATABLE :: resist_b_surft
                                  ! (1/CH-1/CD_STD)/VSHR on land tiles for
                                  ! CLASSIC aerosol scheme
REAL,DIMENSION(:,:,:),ALLOCATABLE :: r_b_dust_ij
                                  ! surf layer res for mineral dust
REAL, DIMENSION(:,:), ALLOCATABLE :: cd_std_dust_ij
                                  ! Bulk transfer coef. for momentum,
                                  ! excluding orographic effects
REAL, DIMENSION(:,:), ALLOCATABLE :: u_s_std_surft
                                  ! Surface friction velocity (standard value)

END MODULE aero
#endif
