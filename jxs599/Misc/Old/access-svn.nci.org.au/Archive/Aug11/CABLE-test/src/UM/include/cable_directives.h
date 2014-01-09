!=============================================================
!--- this file should be included in any src file which attempts to use 
!--- any of the below defs or else's. include UPDEF here #ifdef CABLE
!=============================================================
#define ONLINE_UM 

!#define CABLE_SOIL_LAYERS 
!#define CABLE_17TILES

#  define FRUNTIME_VARS_FILE '/home/599/jxs599/Public/cable_runtime_vars.txt'            
!--- this also has to be set
!check where this still exists and can it be seen
#  define JHNTILES 17


!--- met%qvair is artificially restricted to be in a certain range 
!--- not currently used
#define FQVAIR_LIMITING "on"

!---this is only here for now as cable_17tiles is scattered thru the code
!#define cable_17tiles

!--- def number soil types (cable_....)
#  define JHNSOIL 10
#  define JHNTILES 17



!--- cable_common vars
#define cable_common

!=============================================================




!=============================================================
!!--- USE CABLE_DIAG  
!=============================================================
!--- below directives enable generic build/use of diag 
!--- : cable_diag.F90 - includes module in build
!
!#define cable_diag
!#ifdef cable_diag
#define inc_cable_diag
!#endif

