!=== cable fpp directives 
!
!=============================================================
!--- this file should be included in any src file which attempts to use 
!--- any of the below defs or else's. 
!=============================================================
!
!--- *** DEFS TO ENABLE MISC. CHANGES TO CODE ***
! used to use cable_common_mod in other SRs
#define jhan_anychanges
!
!*** triggers
!
#  ifdef jhan_anychanges
!
!=============================================================
!--- defs to enable misc./temp. changes to code
! used to use cable_common_mod in other SRs
! (in some casses) this can prb. be replaced  entirely later w #ifdef CABLE
!when CABLE is defined in UPDEFS thru umui
!
#define cable_common
!
! used to print total_nsteps at end of atm_step.F90
!
#define jhan_atm_step
!
! used to print total_nsteps at end of atmos_physics2.F90
!
#define jhan_atmos_phys2
!
!used to print misc in cable_cbm.F90 
!
!#define jhan_cablecbm 
!
!used to print ktau_gl in hydrol.F90
!
!#define jhan_hydrol
!
!used to test _directives.h in /include/ dir. can be seen 
!
#define jhan_cable_Exum
!
!used to send appropriate section of d1
!
!#define jhan_atmstep_d1
!
!used to artificially increment solar constant cable_radiation.F90
!
!#define jhan_solcon
!
!=============================================================
!
!!--- USE CABLE_DIAG  
!=============================================================
!--- below directives enable generic use of diag 
!--- : 
!---
!--- : cable_diag.F90 - includes module in build
!
#define cable_diag_cable_diag
!
!*** triggers
!
#  ifdef cable_diag_cable_diag
!#     define cable_diag_cbm_ssoil
#     define cable_diag_atm_step
#     define cable_diag_cable_checks
#  endif
!
!--- SOIL LAYERS
!
!#define jhan_soil_layers
#ifdef jhan_soil_layers
!=============================================================
!--- below directives enable debug in the specific file
!--- : 
!---
!--- : cable_cbm.F90
!#define cable_diag_cbm
!
!--- within _cbm.F90 there are _debugs which may be turned on/off 
!
!--- : atmosphere/landsurfacce/hydrol.F90 
!
!#define cable_diag_hydrol
!
!--- : atmosphere/cable/
!
!#define cable_diag_cable_I
!
!--- for soil _diags
!
!#define cable_diag_cable_E
!#define cable_diag_cable_H
!
!--- : control/top_level/atmos_physics2.F90 
!
!#define cable_diag_atmos_phys2
!
!============================================================
!
!
!============================================================
!--- how many layers
!
#  define six_soil_layers
!*** triggers
!
#     ifdef six_soil_layers
!     else uses four layers
!
!---  this was once used rcf_headaddress to change SoilDepths, however i now 
!---  think that!---  the '4' here  might be an index and not a # of levels var. 
!
!#    define six_soil_depths
!#    define six_soil_depths_inancila
!
#     endif
!
!--- which layers are written
!
#  define top_layers 
#  define bottom_layers 
!=============================================================
#endif
!
!
!--- NTILES
!
!#define jhan_ntiles
#ifdef jhan_ntiles
!=============================================================
!--- below directives enable debug in the specific file
!--- : 
!---
!--- for ntiles _diags
#  define cable_diag_cable_Ex_ntiles
!
!---
!
#  define jhan_17_tiles
!---else reverts to 9 tile format
!
!=============================================================
#endif
!
!
!=============================================================
!--- def to enable jhfort make process
!#define jhfort
!=============================================================
!
#endif

