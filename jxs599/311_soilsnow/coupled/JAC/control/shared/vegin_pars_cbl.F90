MODULE vegin_pars_mod_cbl 

USE grid_constants_mod_cbl, ONLY:                                                 &
          ntype_max, & ! # veg types [13],non-veg=4,ntiles=17
          nrb,       & ! # spectral bANDS VIS/NIR/(LW-not used)
          nsl,       & ! # soil layers
          nscs,      & ! # soil carbon stores
          nvcs         ! # vegetation carbon stores

IMPLICIT NONE

PUBLIC :: vegin_type
PUBLIC :: vegin

!-----------------------------------------------------------------------------
! Description:
!   Defines variable types and variables for CABLE standalone runs.
!   Based on cable_def_types_mod.F90 from the CABLE trunk.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in CABLE SCIENCE
!-----------------------------------------------------------------------------

! Vegetation parameters I/O:
TYPE vegin_type
  REAL ::                                                                     &
       canst1(ntype_max),                                                     &
       dleaf(ntype_max),                                                      &
       length(ntype_max),                                                     &
       width(ntype_max),                                                      &
       vcmax(ntype_max),                                                      &
       ejmax(ntype_max),                                                      &
       hc(ntype_max),                                                         &
       xfang(ntype_max),                                                      &
       rp20(ntype_max),                                                       &
       rpcoef(ntype_max),                                                     &
       rs20(ntype_max),                                                       &
       wai(ntype_max),                                                        &
       rootbeta(ntype_max),                                                   &
       shelrb(ntype_max),                                                     &
       vegcf(ntype_max),                                                      &
       frac4(ntype_max),                                                      &
       xalbnir(ntype_max),                                                    &
       extkn(ntype_max),                                                      &
       tminvj(ntype_max),                                                     &
       tmaxvj(ntype_max),                                                     &
       vbeta(ntype_max),                                                      &
       a1gs(ntype_max),                                                       &
       d0gs(ntype_max),                                                       &
       alpha(ntype_max),                                                      &
       convex(ntype_max),                                                     &
       cfrd(ntype_max),                                                       &
       gswmin(ntype_max),                                                     &
       conkc0(ntype_max),                                                     &
       conko0(ntype_max),                                                     &
       ekc(ntype_max),                                                        &
       eko(ntype_max),                                                        &
       g0(ntype_max),                                                         &
       g1(ntype_max),                                                         &
       zr(ntype_max),                                                         &
       clitt(ntype_max),                                                      &
       froot(nsl,ntype_max),                                                  &
       csoil(nscs,ntype_max),                                                &
       ratecs(nscs,ntype_max),                                               &
       cplant(nvcs,ntype_max),                                               &
       ratecp(nvcs,ntype_max),                                               &
       refl(nrb,ntype_max),                                                   &
       taul(nrb,ntype_max)
END TYPE vegin_type

!Instantiate type here is fine in this case
TYPE(vegin_type) :: vegin !read from namelist

END MODULE vegin_pars_mod_cbl 


