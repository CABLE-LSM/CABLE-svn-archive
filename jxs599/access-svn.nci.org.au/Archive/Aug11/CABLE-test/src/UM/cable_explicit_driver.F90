#include "include/cable_directives.h"

subroutine cable_explicit_driver(                      &
#                 include "include/args/um_expl_driver.h"
                &  )

   use cable_um_tech_mod, only : cable_um_runtime_vars
   use cable_common_module, only : cable_runtime, cable_user
   use cable_diag_module, only : cable_stat
   use io_variables, only: filename
   use cable_um_init_mod
   use cbm_module, only : cbm
   use cable_variables 

   implicit none
   
   !___ declare IN args from sf_exch()
#  include "include/decs/um_expl_driver.h"   
     
   !___ declare local vars 
   !------------------------------------- 
   !___ location of namelist file defining runtime vars
   character(len=200), parameter ::   & 
         runtime_vars_file = '/home/599/jxs599/CABLE_standalone/data/cable.nml'

   !___ 1st call here is necessary in UM & coupled 
   logical :: first_cable_call = .true.
   integer :: itimestep 
   !------------------------------------- 

      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cable_explicit_driver')

      !--- timestep width passed as integer to interface_UM_data 
      itimestep = floor(timestep) 
      
      !--- initialize cable_runtime% switches 
      if(first_cable_call) & 
         call init_runtime_switches()

      !--- internal FLAGS def. specific call of CABLE from UM
      !--- from cable_common_module
      cable_runtime%um_explicit = .true.

      !--- user FLAGS, variables etc def. in cable.nml is read on 
      !--- first time step of each run. these variables are read at 
      !--- runtime and for the most part do not require a model rebuild.
      if(first_cable_call) then
         call cable_um_runtime_vars(runtime_vars_file) 
         first_cable_call = .false.
      endif      
     
      !---------------------------------------------------------------------!
      !--- initialize CABLE using UM forcings etc. these args are passed ---!
      !--- down from UM.                                                 ---! 
      !---------------------------------------------------------------------!
      call interface_UM_data( & 
#              include "include/args/um_expl_pack.h"   
             & )  

      !---------------------------------------------------------------------!
      !--- real(timestep) width, CABLE types passed to CABLE "engine" as ---!  
      !--- req'd by Mk3L  --------------------------------------------------!
      !---------------------------------------------------------------------!
      CALL cbm(TIMESTEP, air, bgc, canopy, met, bal, &
                  rad, rough, soil, ssoil, sum_flux, veg )

      !---------------------------------------------------------------------!
      !--- pass land-surface quantities calc'd by CABLE in explicit call ---!
      !--- back to UM. #include contains UM vars which are updated, other---!
      !--- args are CABLE vars which contain these quantities.           ---!
      !---------------------------------------------------------------------!
      call cable_expl_unpack(   &
#              include "include/args/um_expl_unpack.h"   
              & , ssoil%snowd, ssoil%cls, air%rlam, air%rho, canopy%fe,        &
                canopy%fh, canopy%us, canopy%cdtq, canopy%fwet, canopy%wetfac_cs, & 
                canopy%rnet, canopy%zetar, canopy%epot, met%ua, rad%trad,         &
                rad%transd, rough%z0m, rough%zref_tq )

      cable_runtime%um_explicit = .false.

   return
end subroutine cable_explicit_driver


subroutine init_runtime_switches()
   use cable_common_module, only : cable_runtime
      cable_runtime%um = .true.
      cable_runtime%um_explicit = .false.
      cable_runtime%um_implicit = .false.
      cable_runtime%um_radiation = .false.
      cable_runtime%offline = .false.
      cable_runtime%mk3l = .false.
   return      
end subroutine init_runtime_switches



!---------------------------------------------------------------------!
!--- pass land-surface quantities calc'd by CABLE in explicit call ---!
!--- back to UM. #include contains UM vars which are updated, other---!
!--- args are renamed CABLE vars which contain these quantities.   ---!
!---------------------------------------------------------------------!

subroutine cable_expl_unpack( &
#                       include "include/args/um_expl_unpack.h"   
                 &, ssoil_snowd, ssoil_cls, air_rlam, air_rho, canopy_fe, &
                  canopy_fh, canopy_us, canopy_cdtq, canopy_fwet, &
                  canopy_wetfac_cs, canopy_rnet, canopy_zetar, canopy_epot, &
                  met_ua, rad_trad, rad_transd, rough_z0m, rough_zref_tq )

   use define_dimensions, only : mp 
   use physical_constants, only : niter, capp
   use cable_um_tech_mod, only : um1
   use cable_common_module, only : cable_runtime, cable_user, &
                                    ktau_gl, knode_gl 
   use cable_diag_module, only : cable_stat
   implicit none         

   !___ decs of UM variables to recieve unpacked CABLE vars
#  include "include/decs/um_expl_unpack.h"   
   !___ decs of intent(in) CABLE variables to be unpacked
#  include "include/decs/cable_expl_unpack.h"   
        
   !___vars in local calc. of latent heat fluxes
   REAL, dimension(um1%land_pts,um1%ntiles) ::                  &
      FQW_TILE_CAB,  &
      LE_TILE

   !___vars in local calc of Surface friction velocities
   REAL, dimension(um1%land_pts,um1%ntiles) ::                  &
      CD_CAB_TILE,   &  
      CH_CAB_TILE,   &  ! (bulk transfer) coeff. for momentum
      U_S_TILE
   REAL, dimension(mp)  :: &
      CDCAB,CHCAB

   !___local miscelaneous
   real, dimension(mp)  :: &
   THETAST,fraca_cab,rfsfs_cab, RECIPLMOTILE, fe_dlh
   integer :: i,j,k,N,L
   real :: miss = 0.0
   logical :: first_cable_call = .true.

      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('cable_expl_unpack')
      
      !___return fluxes
      FTL_TILE_CAB = unpack(canopy_fh,  um1%l_tile_pts, miss)
      FTL_CAB = SUM(um1%TILE_FRAC * FTL_TILE_CAB,2)
      FQW_TILE_CAB = unpack(canopy_fe,  um1%l_tile_pts, miss)
      LE_TILE_CAB = unpack(canopy_fe,  um1%l_tile_pts, miss)
      LE_CAB = SUM(um1%TILE_FRAC * LE_TILE_CAB,2)
      fe_dlh = canopy_fe/(air_rlam*ssoil_cls)
      FTL_TILE = unpack(canopy_fh,  um1%l_tile_pts, miss)
      FTL_TILE = FTL_TILE / capp
      FQW_TILE = unpack(fe_dlh, um1%l_tile_pts, miss)
      
      !___return temp and roughness
      TSTAR_TILE_CAB = unpack(rad_trad, um1%l_tile_pts, miss)
      TSTAR_CAB = SUM(um1%TILE_FRAC * TSTAR_TILE_CAB,2)
      TSTAR_TILE = unpack(rad_trad,  um1%l_tile_pts, miss)
      Z0M_TILE = unpack(rough_z0m,  um1%l_tile_pts, miss)
      Z0H_TILE = Z0M_TILE
      
      !___return friction velocities/drags/ etc
      U_S_TILE  =  unpack(canopy_us, um1%l_tile_pts, miss)
      U_S_CAB  = SUM(um1%TILE_FRAC *  U_S_TILE,2)
      CDCAB = canopy_us**2/met_ua**2   ! met%ua is always above umin = 0.1m/s
      ! for Cable CD*
      CD_CAB_TILE =  unpack(CDCAB,um1%l_tile_pts, miss)
      CD_CAB= SUM(um1%TILE_FRAC * CD_CAB_TILE,2)
      ! for Cable CH*
      CH_CAB_TILE =  unpack(canopy_cdtq,um1%l_tile_pts, miss)
      CH_CAB= SUM(um1%TILE_FRAC * CH_CAB_TILE,2)

      U_S_STD_TILE=U_S_TILE
      CD_TILE = CD_CAB_TILE
      CH_TILE = CH_CAB_TILE

      U_S = 0.
      DO N=1,um1%ntiles
      DO K=1,um1%TILE_PTS(N)
        L = um1%TILE_INDEX(K,N)
        J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
        I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
        U_S(I,J) = U_S(I,J)+FLAND(L)*um1%TILE_FRAC(L,N)*U_S_TILE(L,N)
      ENDDO
      ENDDO


      !___return miscelaneous 
      fraca_cab = canopy_fwet * (1.-rad_transd)
      where ( ssoil_snowd > 1.0) fraca_cab = 1.0
      rfsfs_cab = min(1.,max(0.01,canopy_wetfac_cs - fraca_cab)/max(0.01,1. - fraca_cab) )
      FRACA = unpack(fraca_cab, um1%l_tile_pts, miss)
      RESFT =  unpack( canopy_wetfac_cs,um1%l_tile_pts, miss)
      RESFS = unpack( rfsfs_cab , um1%l_tile_pts, miss)

      RADNET_TILE = unpack( canopy_rnet , um1%l_tile_pts, miss)
      THETAST = abs(canopy_fh)/(air_rho*capp*canopy_us)
      RECIPLMOTILE =  canopy_zetar(:,niter) / rough_zref_tq
      RECIP_L_MO_TILE = unpack(RECIPLMOTILE,um1%l_tile_pts, miss)
      EPOT_TILE = unpack(canopy_epot,um1%l_tile_pts, miss)
      
      if(first_cable_call) then 
         !open(unit=713941,file='i3.dat', status="replace",action="write" )
         !   write (713941,*) um1%L_tile_pts 
         !close(713941)
         l_tile_pts = um1%l_tile_pts
         first_cable_call = .false.
      endif

   return
end subroutine cable_expl_unpack
    
!============================================================================
!============================================================================
!============================================================================

