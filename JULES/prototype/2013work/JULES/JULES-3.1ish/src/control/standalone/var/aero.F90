#if !defined(UM_JULES)
! Module containing variables required for aerosol calculations

  MODULE aero

  REAL                              :: CO2_MMR = 5.24100e-4
!                                    ! CO2 Mass Mixing Ratio
  REAL, DIMENSION(:,:), ALLOCATABLE :: CO2_3D
!                                    ! 3D CO2 field if required
  REAL, DIMENSION(:,:), ALLOCATABLE :: RHO_CD_MODV1
!                                    ! Surface air density * drag coef * mod(v1 - v0) before interp
  REAL, DIMENSION(:,:), ALLOCATABLE :: RHO_ARESIST
!                                    ! RHOSTAR*CD_STD*VSHR for Sulphur cycle
  REAL, DIMENSION(:,:), ALLOCATABLE :: ARESIST
!                                    ! 1/(CD_STD*VSHR) for Sulphur cycle
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESIST_B
!                                    ! (1/CH-1/(CD_STD)/VSHR for Sulphur cycle
  REAL, DIMENSION(:,:), ALLOCATABLE :: RHO_ARESIST_TILE
!                                    ! RHOSTAR*CD_STD*VSHR on land tiles
  REAL, DIMENSION(:,:), ALLOCATABLE :: ARESIST_TILE
!                                    ! 1/(CD_STD*VSHR) on land tiles
  REAL, DIMENSION(:,:), ALLOCATABLE :: RESIST_B_TILE
!                                    ! (1/CH-1/CD_STD)/VSHR on land tiles
  REAL,DIMENSION(:,:,:),ALLOCATABLE :: R_B_DUST
!                                    ! surf layer res for dust
  REAL, DIMENSION(:,:), ALLOCATABLE :: CD_STD_DUST
!                                    ! Bulk transfer coef. for momentum, excluding orographic effects
  REAL, DIMENSION(:,:), ALLOCATABLE :: U_S_STD_TILE
!                                    ! Surface friction velocity (standard value)


  NAMELIST /jules_aero/ co2_mmr

  END MODULE aero
#endif
