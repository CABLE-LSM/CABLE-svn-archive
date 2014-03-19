#if !defined(UM_JULES)
! Module containing logical switches for diagnostic outputs

  MODULE diag_swchs

  LOGICAL, PARAMETER :: SQ1P5=.TRUE.             ! Flag for Q1P5M
  LOGICAL, PARAMETER :: ST1P5=.TRUE.             ! Flag for T1P5M
  LOGICAL, PARAMETER :: SU10=.TRUE.              ! Flag for U10M
  LOGICAL, PARAMETER :: SV10=.TRUE.              ! Flag for V10M
  LOGICAL, PARAMETER :: SFME=.TRUE.              ! Flag for FME
  LOGICAL, PARAMETER :: SIMLT=.TRUE.             ! Flag for SICE_MLT_HTF
  LOGICAL, PARAMETER :: SMLT=.TRUE.              ! Flag for SNOMLT_SURF_HTF
  LOGICAL, PARAMETER :: SLH=.TRUE.               ! Flag for LATENT_HEAT
  LOGICAL, PARAMETER :: STF_HF_SNOW_MELT=.TRUE.  ! Flag for snowmelt heat flux
  LOGICAL, PARAMETER :: STF_SUB_SURF_ROFF=.TRUE. ! Flag for sub-surface runoff  
  LOGICAL, PARAMETER :: SZ0HEFF = .TRUE.         ! Stash flag for Z0H_EFF

  END MODULE diag_swchs
#endif
