#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE OASIS4_atm_data_mod
!
! Description: This data module contains items required
!              at various points when coupling the UM atmosphere
!              using OASIS4. 
!
! History
! -------
! Version    Date      Description
! -------- ----------  --------------------------------------
!  6.4      Jan 2007   Original code. R. Hill
!============================================================
      IMPLICIT NONE
        
      ! Variables needed while we await a true conservative
      ! coupler. These allow us to decompose over up to 5
      ! PEs, using a 1xN decomposition and STILL obtain
      ! bit comparable results despite the fact that 
      ! Pre 2007 OASIS4 was not bit comparable (or even
      ! closely comparable) when employing different numbers
      ! of PEs. The work around has been to give 2 dummy
      ! antarctic rows to PEs 0 -> NPROC-2 and the rest of
      ! the domain to PE NPROC-1.
      INTEGER :: O4PE  ! MYPE in O4 context operations
      INTEGER :: O4CNTLPE  ! The OASIS4 pe where all the real
                           ! data is dealt with.
      INTEGER :: O4NPROC   ! Same as NPROC - accessible to 
                           ! special routines.
      INTEGER :: O4_FLD_TYPE ! Set as fld_type_p/u/v as
                             ! appropriate
      INTEGER :: OASIS4_COMP_ID ! Component model ID

      ! Test arrays - ultimately we'll remove these
      REAL, DIMENSION(:,:), ALLOCATABLE :: TEST_VALUES
      REAL, DIMENSION(:,:), ALLOCATABLE :: TEST_WME
      REAL, DIMENSION(:,:), ALLOCATABLE :: TEST_OCN_SST

      ! The UM atmos doesn't have land sea masks available
      ! to us as a matter of course. So we have to go
      ! through the agony of calculating these ourselves
      ! on each grid point type. Here, we define the arrays
      ! which hold the masks. Even then, that's not the whole
      ! story because we may need to cater for coastal tiling.

      ! Mask arrays, local. Whe have two versions of each mask
      ! because PRISM/OASIS4 and the UM atmosphere use opposite
      ! conventions for defining land (the UM has a land point
      ! set to "TRUE", OASIS4 has a land point set to "FALSE").
      ! OASIS4 prefixes indicate masks which are OASIS4 compatible
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS4_TMASK
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS4_UMASK
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: OASIS4_VMASK
        
      ! UM  prefixes indicate masks which are UM compatible
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_TMASK
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_UMASK
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: UM_VMASK

      ! Mask arrays, global.
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GMASKU
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GMASKV

      ! The current prism timestep.
      REAL :: PRISM_TIMESTEP

      ! Indices used to index coupling field arrays.
      ! Outgoing data fields
      INTEGER , PARAMETER :: VIND_TOT_ATM_HFLUX = 1
      INTEGER , PARAMETER :: VIND_TAUX = 2
      INTEGER , PARAMETER :: VIND_TAUY = 3
      INTEGER , PARAMETER :: VIND_WME = 4
      INTEGER , PARAMETER :: VIND_PEN_SOLAR = 5
      INTEGER , PARAMETER :: VIND_PME = 6
      INTEGER , PARAMETER :: VIND_RUNOFF = 7
      INTEGER , PARAMETER :: VIND_SNOW = 8
      INTEGER , PARAMETER :: VIND_SUBLIM = 9
      INTEGER , PARAMETER :: VIND_TOPMELT = 10
      INTEGER , PARAMETER :: VIND_BOTMELT = 11
      INTEGER , PARAMETER :: VIND_CO2 = 12

      ! INCOMING DATA FIELDS
      INTEGER , PARAMETER :: VIND_OCN_SST = 20
      INTEGER , PARAMETER :: VIND_OCN_U = 21
      INTEGER , PARAMETER :: VIND_OCN_V = 22   
      INTEGER , PARAMETER :: VIND_OCN_FREEZE = 23

      INTEGER :: VAR_IND
END MODULE OASIS4_atm_data_mod
#endif
