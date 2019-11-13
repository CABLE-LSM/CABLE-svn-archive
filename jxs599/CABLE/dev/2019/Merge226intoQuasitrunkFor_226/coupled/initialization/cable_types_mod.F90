MODULE cable_types_mod

USE cable_other_constants_mod,  ONLY: nsl, nscs, nvcs, nrb

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------
INTEGER, parameter :: mstype=9      ! # total no of soil types _ this shou;ld
INTEGER, parameter :: msn_cable =3  ! # total no of snow layers

INTEGER :: mp                       ! # total no of patches/tiles

INTEGER, parameter :: npft_max = 17 !cable_veg_params=17 !nsurft

integer, allocatable :: SurfaceTypeID_cbl(:) !veg%iveg
 
!Vegetation parameters
REAL, allocatable  :: VegXfang(:)   !leaf angle PARAMETER (veg%xfang)
REAL, allocatable  :: VegTaul(:,:)  !PARAMETER leaf transmisivity (veg%taul)
REAL, allocatable  :: VegRefl(:,:)  !PARAMETER leaf reflectivity (veg%refl)

LOGICAL, ALLOCATABLE :: l_tile_pts(:,:)

TYPE vegin_type
  REAL ::                                                                  &
       canst1(npft_max),                                                      &
       dleaf(npft_max),                                                       &
       length(npft_max),                                                      &
       width(npft_max),                                                       &
       vcmax(npft_max),                                                       &
       ejmax(npft_max),                                                       &
       hc(npft_max),                                                          &
       xfang(npft_max),                                                       &
       rp20(npft_max),                                                        &
       rpcoef(npft_max),                                                      &
       rs20(npft_max),                                                        &
       wai(npft_max),                                                         &
       rootbeta(npft_max),                                                    &
       shelrb(npft_max),                                                      &
       vegcf(npft_max),                                                       &
       frac4(npft_max),                                                       &
       xalbnir(npft_max),                                                     &
       extkn(npft_max),                                                       &
       tminvj(npft_max),                                                      &
       tmaxvj(npft_max),                                                      &
       vbeta(npft_max),                                                       &
       a1gs(npft_max),                                                        &
       d0gs(npft_max),                                                        &
       alpha(npft_max),                                                       &
       convex(npft_max),                                                      &
       cfrd(npft_max),                                                        &
       gswmin(npft_max),                                                      &
       conkc0(npft_max),                                                      &
       conko0(npft_max),                                                      &
       ekc(npft_max),                                                         &
       eko(npft_max),                                                         &
       g0(npft_max),                                                          &
       g1(npft_max),                                                          &
       zr(npft_max),                                                          &
       clitt(npft_max),                                                       &
       froot(nsl,npft_max),                                                   &
       csoil(nscs,npft_max),                                                  &
       ratecs(nscs,npft_max),                                                 &
       cplant(nvcs,npft_max),                                                 &
       ratecp(nvcs,npft_max),                                                 &
       refl(nrb,npft_max),                                                    &
       taul(nrb,npft_max)
END TYPE vegin_type

TYPE(vegin_type),  SAVE  :: vegin


TYPE soilin_type
  REAL ::                                                                     &
       silt(mstype),                                                        &
       clay(mstype),                                                        &
       sand(mstype),                                                        &
       swilt(mstype),                                                       &
       sfc(mstype),                                                         &
       ssat(mstype),                                                        &
       bch(mstype),                                                         &
       hyds(mstype),                                                        &
       sucs(mstype),                                                        &
       rhosoil(mstype),                                                     &
       css(mstype)
END TYPE soilin_type

TYPE(soilin_type), SAVE :: soilin

END MODULE cable_types_mod
