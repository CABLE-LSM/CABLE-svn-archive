
module cable_um_init_subrs
   implicit none

   contains
!========================================================================
!========================================================================
!========================================================================                                   
         
   subroutine initialize_maps(latitude,longitude, tile_index_mp)
      use cable_data_module, only : cable, const 
      use cable_um_tech_mod, only : um1
      use define_dimensions, only : mp

      use cable_diag_module, only : cable_diag 
      use cable_common_module, only : ktau_gl, knode_gl, cable_user 
         
      implicit none
      real, intent(in), dimension(um1%row_length,um1%rows) :: &
         latitude, longitude
      integer, intent(in), dimension(um1%land_pts, um1%ntiles) :: &
         tile_index_mp  ! index of tile
          
      logical, save :: first_call = .true.
      
      INTEGER :: j

      real :: dlon
      real, dimension(um1%row_length) :: tlong,acoslong
      real, dimension(um1%row_length, um1%rows) :: new_longitude


           
            allocate( cable%lat(mp), cable%lon(mp), cable%tile(mp), cable%tile_frac(mp) )

            !-------------------------------------   
            !---make indexes for tile, lat, lon
            !-------------------------------------   
           
            !--- get latitude index corresponding to cable points
            call um2cable_rr( (asin(latitude)/const%math%pi180), cable%lat )

            !--- get longitude index corresponding to cable points.
            !--- this is not so straight forward as UM longitude index 
            !--- contains ambiguity. thus define "new_longitude" first
            acoslong =  acos( longitude(:,1) ) /const%math%pi180  
       
            tlong(1) = acoslong(1)
            do j=2, um1%row_length
               if( acoslong(j) < acoslong(j-1) ) then  
                  dlon = acoslong(j) - acoslong(j-1)
                  tlong(j) = tlong(j-1) - dlon   
               else 
                  tlong(j) = acoslong(j)
               endif           
            enddo
            
            do j=1, um1%row_length
               new_longitude(j,:) = tlong(j)
            enddo
            
            call um2cable_rr( new_longitude, cable%lon )
         
         
            !--- get tile index/fraction  corresponding to cable points
            cable%tile = pack(tile_index_mp, um1%l_tile_pts)
            cable%tile_frac = pack(um1%tile_frac, um1%l_tile_pts)

         !--- write all these maps.  cable_user%initialize_mapping can be 
         !--- set in namelist cable.nml
         if ( cable_user%initialize_mapping ) then
            !write indexes for tile, lat, lon
            call cable_diag( 1, 'latitude', um1%rows, 1, ktau_gl,  & 
                  knode_gl, 'latitude', ( asin( latitude(1,:) ) /const%math%pi180 ) ) 
            
            call cable_diag( 1, 'longitude', um1%row_length, 1, ktau_gl,  & 
                  knode_gl, 'longitude', ( new_longitude(:,1) ) ) 
        
            !write indexes for tile, lat, lon
            call cable_diag( 1, 'lat_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'lat', cable%lat )
            call cable_diag( 1, 'lon_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'lon', cable%lon )
            
            !this should be integer-ed. typecast for now
            call cable_diag( 1, 'tile_index', mp, 1, ktau_gl,  & 
                  knode_gl, 'tile', real(cable%tile) )
            
            call cable_diag( 1, 'tile_frac', mp, 1, ktau_gl,  & 
                  knode_gl, 'tile_frac', cable%tile_frac )
            
          endif  
         
      
      return
   end subroutine initialize_maps
  
  
         
   subroutine initialize_soil( bexp, hcon, satcon, sathh, smvcst, smvcwt, &
                  smvccl, albsoil, tsoil_tile, sthu, sthu_tile, dzsoil ) 
      use define_dimensions, only : ms, mstype
      use cable_um_tech_mod, only : um1, soil, veg, ssoil 
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user, &
            soilin, get_type_parameters
      implicit none

      real, intent(in), dimension(um1%land_pts) :: &
         bexp, &
         hcon, &
         satcon, & 
         sathh, &
         smvcst, &
         smvcwt, &
         smvccl, &
         albsoil 
      
      real, intent(in), dimension(um1%land_pts, um1%sm_levels) :: sthu
      
      real, intent(in), dimension(um1%land_pts, um1%ntiles, um1%sm_levels) :: &
         sthu_tile,     &
         tsoil_tile

      real, intent(in), dimension(um1%sm_levels) :: dzsoil

      !___defs 1st call to CABLE in this run
      logical, save :: first_call= .true.
      integer :: i,j,k,L,n, logn=6
      real, allocatable:: tempvar(:)
      logical :: vegparmnew=.true., skip =.true. 

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_soil')
            
         if (first_call) then 

            ssoil%pudsto = 0.0; ssoil%pudsmx = 0.0
         
            !--- soil%isoilm defines soiltype. currently is either 2 (arbitrarily) or 9.
            !--- type 9 -> permanent ice points which are dealt with by CABLE, but not  
            !--- UM-MOSES. spatially explicit soil properties are used by the UM anyway,
            !--- and is only really an issue for soil%css & soil%rhosoil, which are set 
            !--- to either 2 or 9. dealing with this in CASACNP is another issue.
            !--- %isoilm=10 for Lakes
            soil%isoilm  =  2
           
            !--- set CABLE-var soil%albsoil from UM var albsoil (see below ~ um2cable_lp)
            call um2cable_lp( albsoil, albsoil, soil%albsoil(:,1), soil%isoilm, skip )

            ! this is incorrect and need to be changed ??
            ! in UM soil%albsoil = 0.75 for permanet ice points!!!
            where ( soil%albsoil(:,1) > 0.7 ) soil%isoilm = 9



            !--- defined in soil_thick.h in UM
            soil%zse = dzsoil
            soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints
            soil%zshh(ms+1)=0.5*soil%zse(ms)
            soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))


            
            !--- read in soil (and veg) parameters 
            !--- logn, vegparmnew can be set thru cable.nml
            !--- logn=6(write to std out)
            !--- vegparmnew=.true. (read std veg params not CASA file) 
            call  get_type_parameters(logn,vegparmnew)



            !-------------------------------------------------------------------
            !--- UM met forcing vars needed by CABLE which have UM dimensions
            !---(land_pts,ntiles)[_lp], which is no good to cable. These have to be 
            !--- re-packed in a single vector of active tiles. Hence we use 
            !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
            !--- if the land point is/has an active tile. generic format:
            !---     um2cable_lp( UM var, 
            !---                 default value for snow tile, 
            !---                 peranent ice point to be treatedas a snowtile, 
            !---                 CABLE var, 
            !---                 mask )
            !--- where mask tells um2cable_lp whether or not to use default value 
            !--- for snow tile 
            !-------------------------------------------------------------------
            
            ! parameter b in Campbell equation 
            call um2cable_lp( BEXP, soilin%bch, soil%bch, soil%isoilm)
            
            allocate( tempvar(um1%land_pts) )
            tempvar = soilin%sand(9)*0.3  + soilin%clay(9)*0.25 + &
                        soilin%silt(9)*0.265
            call um2cable_lp( HCON, tempvar, soil%cnsd, soil%isoilm)
            deallocate( tempvar )
            ! hydraulic conductivity @ saturation (satcon [mm/s], soilin%hyds[m/s] )
            call um2cable_lp( satcon,soilin%hyds*1000.0, soil%hyds, soil%isoilm)
            call um2cable_lp( sathh, soilin%sucs, soil%sucs, soil%isoilm)
            call um2cable_lp( smvcst, soilin%ssat, soil%ssat, soil%isoilm)
            call um2cable_lp( smvcwt, soilin%swilt, soil%swilt, soil%isoilm)
            call um2cable_lp( smvccl, soilin%sfc, soil%sfc, soil%isoilm)
   
    
            
            !--- (re)set values for CABLE
            !-------------------------------------------------------------------
     
            soil%ibp2    =  NINT(soil%bch)+2
            soil%i2bp3   =  2*NINT(soil%bch)+3
            !--- satcon in UM is in mm/s; Cable needs m/s
            soil%hyds    =  soil%hyds / 1000.
            soil%sucs    =  abs( soil%sucs )
            !jhan:coupled runs 
            soil%hsbh    =  soil%hyds*ABS(soil%sucs)*soil%bch
            
            where (soil%ssat > 0. ) soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
           
            !--- these are temporary 
            soil%rhosoil =  soilin%rhosoil(soil%isoilm)
            soil%css     =  soilin%css(soil%isoilm)
         
            
            first_call= .false.
         endif

      return
   end subroutine initialize_soil
 
!========================================================================
!========================================================================
!========================================================================
          
   subroutine initialize_veg( canht_ft, lai_ft) 
      use cable_um_tech_mod
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user, vegin
      implicit none
      real, intent(in), dimension(um1%land_pts, um1%npft) :: canht_ft, lai_ft 
      !___defs 1st call to CABLE in this run
      logical, save :: first_call= .true.
      integer :: logn=6
      logical :: vegparmnew=.true.

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_veg')
   
         !---clobbers veg height, lai and resets ivegt for CABLE tiles
         call clobber_height_lai( canht_ft, lai_ft )
         
         !--- veg params were read from initialize_soil() 
         if (first_call) &
            call init_veg_pars_fr_vegin() 
         
         first_call= .false.
        
      return
   end subroutine initialize_veg

!========================================================================
!========================================================================
!========================================================================

   subroutine clobber_height_lai( um_htveg, um_lai )
      use cable_um_tech_mod, only : um1, kblum_veg, veg
      implicit none

      real, intent(in), dimension(um1%land_pts, um1%npft) :: &
                                                      um_htveg, um_lai
      integer :: i,j,n
       
      DO N=1,um1%NTILES
         DO J=1,um1%TILE_PTS(N)
            i = um1%TILE_INDEX(j,N)  ! It must be landpt index
            if( um1%TILE_FRAC(i,N) .gt. 0.0 ) then
               if(N < 4 ) then
                  kblum_veg%IVEGT(i,N) = N
                  kblum_veg%LAIFT(i,N) = max(0.01,um_lai(i,N)) 
                  kblum_veg%HTVEG(i,N) = max(1.,um_htveg(i,N)) 
               else if(N > 3 .AND. N < 14 ) then
                  kblum_veg%IVEGT(i,N) = N
                  kblum_veg%LAIFT(i,N) = max(0.01, um_lai(i,N)) 
                  kblum_veg%HTVEG(i,N) = max(0.1, um_htveg(i,N)) 
                else if(N > 13 ) then
                  kblum_veg%IVEGT(i,N) = N
                  kblum_veg%LAIFT(i,N) = 0. 
                  kblum_veg%HTVEG(i,N) = 0.
               endif
            endif
         ENDDO
      ENDDO
     
      veg%iveg   = pack(kblum_veg%ivegt, um1%L_TILE_PTS)
      veg%vlai   = pack(kblum_veg%laift, um1%L_TILE_PTS)
      veg%hc     = pack(kblum_veg%htveg, um1%L_TILE_PTS)

      return
   end subroutine clobber_height_lai

!========================================================================
!========================================================================
!========================================================================

   subroutine init_veg_pars_fr_vegin() 
      use cable_common_module, only : vegin
      use cable_um_tech_mod, only : veg, soil 
      use define_dimensions, only : mp
      implicit none

      integer :: k

         !jhan:UM reads from ancil. & resets thru kblum_veg   
         !jh:veg%hc      = vegin%hc(veg%iveg)
         veg%canst1  = vegin%canst1(veg%iveg)
         veg%ejmax   = 2.*vegin%vcmax(veg%iveg)
         veg%frac4   = vegin%frac4(veg%iveg)
         veg%tminvj  = vegin%tminvj(veg%iveg)
         veg%tmaxvj  = vegin%tmaxvj(veg%iveg)
         veg%vbeta   = vegin%vbeta(veg%iveg)
         veg%rp20    = vegin%rp20(veg%iveg)
         veg%rpcoef  = vegin%rpcoef(veg%iveg)
         veg%shelrb  = vegin%shelrb(veg%iveg)
         veg%vegcf   = vegin%vegcf(veg%iveg)
         veg%extkn   = vegin%extkn(veg%iveg)
         veg%vcmax   = vegin%vcmax(veg%iveg)
         veg%xfang   = vegin%xfang(veg%iveg)
         veg%dleaf   = vegin%dleaf(veg%iveg)
         veg%xalbnir = vegin%xalbnir(veg%iveg)
         
         !---cable uses as soil type
         soil%rs20 = vegin%rs20(veg%iveg)

         do k=1,2
           veg%refl(:,k)   = vegin%reflin(k,veg%iveg)
           veg%taul(:,k)   = vegin%taulin(k,veg%iveg)
         enddo

         !jhan:can't we read this in
         veg%froot(:,1) = 0.05
         veg%froot(:,2) = 0.20
         veg%froot(:,3) = 0.20
         veg%froot(:,4) = 0.20
         veg%froot(:,5) = 0.20
         veg%froot(:,6) = 0.15
      return
   end subroutine init_veg_pars_fr_vegin

!========================================================================
!========================================================================
!========================================================================
        
   subroutine initialize_radiation( sw_down, lw_down, cos_zenith_angle, &
                     surf_down_sw, sin_theta_latitude, ls_rain, ls_snow,   &
                     tl_1, qw_1, vshr_land, pstar, co2_mmr ) 
      use define_dimensions, only : mp
      use physical_constants, only : tfrz
      use other_constants, only : RAD_THRESH
      use cable_um_tech_mod, only : um1, rad, soil, met, & 
                              conv_rain_prevstep, conv_snow_prevstep
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none

      real, intent(inout), dimension(um1%row_length, um1%rows) :: sw_down
      
      real, intent(in), dimension(um1%row_length, um1%rows) :: & 
         lw_down, &
         sin_theta_latitude
      real, intent(inout), dimension(um1%row_length, um1%rows) :: cos_zenith_angle

      real, intent(in), dimension(um1%row_length, um1%rows, 4) :: surf_down_sw 
      
      real, intent(in), dimension(um1%row_length, um1%rows) :: & 
         ls_rain, &
         ls_snow, &   
         tl_1,    &
         qw_1,    &
         vshr_land, & 
         pstar
      
      real, intent(in) :: co2_mmr
                
      !___defs 1st call to CABLE in this run. OK in UM & coupled
      logical, save :: first_call= .true.
     
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_radiation')
         
         if (first_call) then
            rad%albedo_T = soil%albsoil(:,1)
            first_call = .false.
            allocate( conv_rain_prevstep(mp), conv_snow_prevstep(mp) )
            conv_rain_prevstep = 0. 
            conv_snow_prevstep = 0.
         endif   
      
         !--- re-set UM rad. forcings to suit CABLE. also called in explicit call to 
         !--- CABLE from subr cable_um_expl_update() 
         call update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
      
         !--- set met. and rad. forcings to CABLE. also called in radiation call to 
         !--- CABLE from subr cable_rad_() !jhan?
         !--- subr.  um2cable_met_rad_alb() USES CABLE types met%, rad%, soil%
         !--- and kblum% rad. calculated in  update_kblum_radiation() above 
         call um2cable_met_rad( cos_zenith_angle)
         
         !--- UM met forcing vars needed by CABLE which have UM dimensions
         !---(row_length,rows)[_rr], which is no good to cable. These have to be 
         !--- re-packed in a single vector of active tiles. Hence we use 
         !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
         !--- if the land point is/has an active tile
         !--- generic format:
         !--- um2cable_rr( UM var, CABLE var)
         
         !--- CABLE met type forcings, not set by um2cable_met_rad()
         call um2cable_rr( LW_DOWN, met%fld)
         call um2cable_rr( (LS_RAIN*um1%TIMESTEP), met%precip)
         call um2cable_rr( (LS_SNOW*um1%TIMESTEP), met%precip_sn)
         call um2cable_rr( TL_1, met%tk)
         call um2cable_rr( QW_1, met%qv)
         call um2cable_rr( VSHR_LAND, met%ua)
         call um2cable_rr( PSTAR*0.01, met%pmb)
      
         !---re-set some of CABLE's forcing variables
         met%precip   =  (met%precip + conv_rain_prevstep) &
                        + (met%precip_sn +  conv_rain_prevstep)
         met%tvair =     met%tk
         met%tvrad =     met%tk
         met%tc =        met%tk - tfrz
         met%coszen =    max(met%coszen,1e-8) 
  
         !---this is necessary clobrring at present 
         where (met%ua < 0.001 ) met%ua = 0.001
         
         ! rml 24/2/11 Set atmospheric CO2 seen by cable to CO2_MMR (value seen by radiation
         !   scheme).  Option in future to have cable see interactive (3d) CO2 field
         !   Convert CO2 from kg/kg to mol/mol (m_air, 28.97 taken from UKCA include file, c_v_m.h)
         met%ca =        CO2_MMR * 28.97/44.
            
         where (met%coszen < RAD_THRESH ) 
            rad%fbeam(:,1) = real(0) 
            rad%fbeam(:,2) = real(0) 
            rad%fbeam(:,3) = real(0) 
         endwhere

         !--- CABLE radiation type forcings, not set by um2cable_met_rad(
         !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
         call um2cable_rr( um1%longitude,rad%longitude )

      return
   end subroutine initialize_radiation

!========================================================================
!========================================================================
!========================================================================
          
   subroutine initialize_canopy(canopy_tile)
      use cable_um_tech_mod, only : um1, canopy 
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
      real, intent(in),dimension(um1%land_pts, um1%ntiles) :: canopy_tile
      !___defs 1st call to CABLE in this run. OK in UM & coupled
      logical, save :: first_call= .true.
      
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_canopy')
         
         !--- %ga is computed (on LHS only) in define_canopy and 
         !--- then used in soilsnow() in implicit call, then unpacked
         if (first_call) then
            canopy%ga = 0.
            first_call = .false.
         endif   
        !---set canopy storage (already in dim(land_pts,ntiles) ) 
        canopy%cansto = pack(CANOPY_TILE, um1%l_tile_pts)
        canopy%oldcansto=canopy%cansto

      return
   end subroutine initialize_canopy

!========================================================================
!========================================================================
!========================================================================
 
        subroutine initialize_soilsnow( smvcst, tsoil_tile, sthf_tile, smcl_tile, &
                        snow_tile, snow_rho1l, snage_tile, isnow_flg3l, snow_rho3l, &
                        snow_cond, snow_depth3l, snow_mass3l, snow_tmp3l, fland, & 
                        sin_theta_latitude ) 

      use define_dimensions, only : mp
      use physical_constants, only : tfrz
      use cable_um_tech_mod, only : um1, soil, ssoil, met, bal, veg
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
      
      real, intent(in), dimension(um1%land_pts) :: smvcst
      
      real, intent(in), dimension(um1%land_pts, um1%ntiles, um1%sm_levels) :: &
         sthf_tile, &   !
         smcl_tile, &   !
         tsoil_tile     !

      integer, intent(in), dimension(um1%land_pts, um1%ntiles) :: isnow_flg3l 

      real, intent(inout), dimension(um1%land_pts, um1%ntiles) :: snow_tile

      real, intent(in), dimension(um1%land_pts, um1%ntiles) :: &
         snow_rho1l, &  !
         snage_tile     !

      real, intent(inout), dimension(um1%land_pts, um1%ntiles,3) :: snow_cond

      real, intent(in), dimension(um1%land_pts, um1%ntiles,3) :: & 
         snow_rho3l, &     !
         snow_depth3l, &   !
         snow_mass3l,  &   !
         snow_tmp3l        !
      
      real, intent(in), dimension(um1%land_pts) :: fland 
      
      real, intent(in), dimension(um1%row_length, um1%rows) :: sin_theta_latitude

      
      
      integer :: i,j,k,L,n
      real  :: zsetot, max_snow_depth=50000.
      real, allocatable:: fwork(:,:,:), sfact(:), fvar(:), rtemp(:)
      logical :: skip =.true. 
      logical :: first_call = .true.

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_soilsnow')
         
         ssoil%wbtot1 = 0
         ssoil%wbtot2 = 0
!jhan:check this with Eva
         snow_tile = min(max_snow_depth, snow_tile)

         ssoil%snowd  = pack(SNOW_TILE,um1%l_tile_pts)
         ssoil%ssdnn  = pack(SNOW_RHO1L,um1%l_tile_pts)  
         ssoil%isflag = pack(ISNOW_FLG3L,um1%l_tile_pts)  
         !jhan:replace hard_wiring   
         do J=1,3  
            ssoil%sdepth(:,J)= pack(SNOW_DEPTH3L(:,:,J),um1%l_tile_pts)
            ssoil%smass(:,J) = pack(SNOW_MASS3L(:,:,J),um1%l_tile_pts)  
            ssoil%ssdn(:,J)  = pack(SNOW_RHO3L(:,:,J),um1%l_tile_pts)  
            ssoil%tggsn(:,J) = pack(SNOW_TMP3L(:,:,J),um1%l_tile_pts)  
            ssoil%sconds(:,J)= pack(SNOW_COND(:,:,J),um1%l_tile_pts)  
            where ( veg%iveg == 16 )
               ssoil%wbtot1 = ssoil%wbtot1 + REAL(ssoil%wb(:,J)) * 1000.0 * soil%zse(J)
               !ssoil%wb(:,J) = 0.95*soil%ssat
               !jhan:coupled!temp fix for lakes
               ssoil%wb(:,J) = soil%sfc
               ssoil%wbtot2 = ssoil%wbtot2 + REAL(ssoil%wb(:,J)) * 1000.0 * soil%zse(J)
            endwhere
         enddo 
          
         do J=1,um1%sm_levels
            ssoil%tgg(:,J) = pack(TSOIL_TILE(:,:,J),um1%l_tile_pts)
         enddo 
         
         ssoil%snage = pack(SNAGE_TILE, um1%l_tile_pts)
         ssoil%wb_lake = max( ssoil%wbtot2 - ssoil%wbtot1, 0.)

         if( first_call) then 
            ssoil%wbtot = 0.
            ssoil%wb_lake = 0.0
            ssoil%tggav = 0.
            ssoil%rtsoil = 50.
            ssoil%t_snwlr = 0.05

            !jhan: snow depth from prev timestep
            ssoil%osnowd  = pack(SNOW_TILE,um1%l_tile_pts)  

            !jhan:things which should not need to be redone if restart is matured
            !:!ssoil%snowd  = min(1100., ssoil%snowd)  
            !:!ssoil%osnowd = min(1100., ssoil%osnowd)  
    
            !!:!WHERE( soil%isoilm == 9 )
            !!:!   ssoil%tggsn(:,3) = min(ssoil%tggsn(:,3),ssoil%tgg(:,1)) ! 
            !!:!   ssoil%tggsn(:,2) = min(ssoil%tggsn(:,2),ssoil%tgg(:,1)) ! 
            !!:!   ssoil%tggsn(:,1) = min(ssoil%tggsn(:,1),tfrz )
            !!:!ENDWHERE
        
            zsetot = sum(soil%zse)
            DO k = 1, um1%sm_levels
               ssoil%tggav = ssoil%tggav  + soil%zse(k)*ssoil%tgg(:,k)/zsetot
            END DO
        
            !--- not updated 
            allocate( sfact( mp ) )
            sfact = 0.68
            WHERE (soil%albsoil(:,1) <= 0.14) 
               sfact = 0.5
            ELSEWHERE (soil%albsoil(:,1) > 0.14 .and. soil%albsoil(:,1) <= 0.20)
              sfact = 0.62
            END WHERE
            ssoil%albsoilsn(:,2) = 2. * soil%albsoil(:,1) / (1. + sfact)
            ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
            deallocate( sfact )
   
            allocate( fvar(um1%land_pts ) )
            DO N=1,um1%NTILES
               DO K=1,um1%TILE_PTS(N)
                  L = um1%TILE_INDEX(K,N)
                  fvar(L) = real(L)
               ENDDO
            ENDDO
            call um2cable_lp( fland, fland, ssoil%fland, soil%isoilm, skip )
            call um2cable_lp( fvar, fvar, ssoil%ifland, soil%isoilm, skip )
            deallocate( fvar )
            
            !--- updated via smcl,sthf etc 
            allocate( fwork(um1%land_pts,um1%ntiles,2*um1%sm_levels) )
            DO N=1,um1%NTILES                                                   
              DO K=1,um1%TILE_PTS(N)                                           
              I = um1%TILE_INDEX(K,N)                                      
                DO J = 1,um1%SM_LEVELS
                  fwork(I,N,J) = SMCL_TILE(I,N,J)/soil%zse(j) / um1%RHO_WATER
                  fwork(I,N,J+um1%SM_LEVELS) = STHF_TILE(I,N,J)*SMVCST(I)
                ENDDO ! J
              ENDDO
            ENDDO
      
            DO J = 1,um1%SM_LEVELS
              ssoil%wb(:,J)  = pack(fwork(:,:,J),um1%l_tile_pts)
              ssoil%wbice(:,J) = pack(fwork(:,:,J+um1%SM_LEVELS),um1%l_tile_pts)
              ssoil%wbice(:,J) = max(0.,ssoil%wbice(:,J))
               where ( veg%iveg == 16 ) &
                  ssoil%wb(:,J) = 0.95*soil%ssat
            ENDDO
            
            deallocate( fwork )
            
            ssoil%owetfac = MAX(0., MIN(1.0, &
                (ssoil%wb(:,1) - soil%swilt) / (soil%sfc - soil%swilt)))
            !Temporay fixer for accounting for reduction of soil evaporation due to freezing
            where ( ssoil%wbice(:,1) > 0. ) &
               ! Prevents divide by zero at glaciated points where wb and wbice=0.
               ! rml 15/2/11 change ssoil%wetfac to ssoil%owetfac, since ssoil%wetfac not defined yet
               ssoil%owetfac = ssoil%owetfac * ( 1.0 - ssoil%wbice(:,1)/ssoil%wb(:,1) )**2
         
            !jhan: do we want to do this before %owetfac is set 
            DO J = 1, um1%sm_levels
            WHERE( soil%isoilm == 9 )
                ssoil%wb(:,J) = 0.95*soil%ssat
                ssoil%wbice(:,J) = 0.8*soil%ssat
            ENDWHERE
              ssoil%wbtot  = ssoil%wbtot + ssoil%wb(:,j) * soil%zse(j) * 1000.0
            ENDDO
        
            bal%wbtot0 = ssoil%wbtot
   
            !---set antartic flag using  sin_theta_latitude(row_length,rows)
            allocate( fwork(1,um1%land_pts,um1%ntiles) )
             DO N=1,um1%NTILES                     
               DO K=1,um1%TILE_PTS(N)
                  L = um1%TILE_INDEX(K,N)
                  J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
                  I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
                  if( sin_theta_latitude(I,J) .lt. -0.91 ) fwork(1,L,N) = 1.0
               ENDDO
            ENDDO
            ssoil%iantrct = pack(fwork(1,:,:),um1%L_TILE_PTS)
           
            deallocate( fwork )
   
            first_call = .false.
         endif        
      return
   end subroutine initialize_soilsnow
 
!========================================================================
!========================================================================
!========================================================================
          
   subroutine initialize_roughness( z1_tq, z1_uv, htveg )  
      use cable_um_tech_mod, only : um1, rough, veg
      use cable_common_module, only : ktau_gl
      use define_dimensions, only : mp
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
      real, intent(in), dimension(um1%row_length, um1%rows) ::  z1_tq, z1_uv
      real, intent(inout), dimension(um1%land_pts, um1%ntiles) :: htveg
      integer :: i,j,k,L,n
      real, allocatable, dimension(:,:) :: jhruff, jhwork

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_roughness')

         !--- CABLE roughness type forcings
         call um2cable_rr( Z1_TQ, rough%za_tq)
         call um2cable_rr( Z1_UV, rough%za_uv)

         allocate(jhwork (um1%land_pts,um1%ntiles) ) 
         allocate(jhruff (um1%land_pts,um1%ntiles) ) 

         !Veg height changes seasonally in MOSES hence no updates here due to snow
         jhwork = 0.
         DO N=1,um1%NTILES
           DO K=1,um1%TILE_PTS(N)
             I = um1%TILE_INDEX(K,N)
             jhWORK(I,N) = MAX(.01,HTVEG(I,N))
           ENDDO
         ENDDO

         jHRUFF= 0.01 
         DO l=1,um1%land_pts
           DO n=1,um1%ntiles     
             IF( jHRUFF(L,N) .lt. jhwork(l,n)) jHRUFF(L,:) =  jhwork(l,n)
           ENDDO
         ENDDO
         
         !jhan: what is the redundancy here veg%hc -> htveg?
         rough%hruff= MAX(0.01,veg%hc)
         rough%hruff_grmx = pack(jHRUFF, um1%l_tile_pts) 

         deallocate(jhruff ) 
         deallocate(jhwork) 
   
      return
   end subroutine initialize_roughness

!========================================================================
!========================================================================
!========================================================================
      subroutine update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
         use cable_um_tech_mod!, only : um1, um_rad, kblum_rad
         implicit none
         real, intent(inout), dimension(um1%row_length, um1%rows) :: sw_down
         real, intent(in), dimension(um1%row_length, um1%rows) :: cos_zenith_angle
         real, intent(in), dimension(um1%row_length, um1%rows, 4) :: surf_down_sw 

            !jhan: do you really want to be changing sw_down 
            !NOTE: in which case passing SW_DOWN is unecessary           
            SW_DOWN = ( surf_down_sw(:,:,1)      &
                              + surf_down_sw(:,:,2)     &
                              + surf_down_sw(:,:,3)     &
                              + surf_down_sw(:,:,4) )   &
                              * cos_zenith_angle(:,:)

            kblum_rad%SW_DOWN_DIR = ( surf_down_sw(:,:,1)      &
                              + surf_down_sw(:,:,3) )          &
                              * cos_zenith_angle(:,:)

            kblum_rad%SW_DOWN_DIF = ( surf_down_sw(:,:,2)      & 
                                    + surf_down_sw(:,:,4) )    &
                                    *cos_zenith_angle(:,:)

            kblum_rad%SW_DOWN_VIS = (surf_down_sw(:,:,1)      & 
                                    + surf_down_sw(:,:,2) )    &
                                    * cos_zenith_angle(:,:)

            kblum_rad%SW_DOWN_NIR = ( surf_down_sw(:,:,3)        &
                                    + surf_down_sw(:,:,4) ) &
                                    *cos_zenith_angle(:,:)
            ! fbeam for VIS
            kblum_rad%FBEAM(:,:,1) = surf_down_sw(:,:,1)         &
                                    * cos_zenith_angle(:,:)      &
                                       / max( 0.1, kblum_rad%SW_DOWN_VIS )
            ! fbeam for NIR
            kblum_rad%FBEAM(:,:,2) = surf_down_sw(:,:,3)         &
                                    * cos_zenith_angle(:,:)  &
                                    / max( 0.1, kblum_rad%SW_DOWN_NIR )
            !---fbeam for all solar 
            kblum_rad%FBEAM(:,:,3) = kblum_rad%SW_DOWN_DIR / &
                                    max( 0.1, SW_DOWN )
             
         return
      end subroutine update_kblum_radiation

!========================================================================
!========================================================================
!========================================================================

   subroutine  um2cable_met_rad( cos_zenith_angle)
      use cable_um_tech_mod, only :um1, kblum_rad, rad, met
      implicit none
      !___ from UM, cosine zenith angle and soil albedo
      real, intent(inout) :: cos_zenith_angle(um1%row_length, um1%rows)

         !--- CABLE met type forcings
         call um2cable_rr( cos_zenith_angle, met%coszen)
         call um2cable_rr( kblum_rad%SW_DOWN_VIS, met%fsd(:,1))
         call um2cable_rr( kblum_rad%SW_DOWN_NIR, met%fsd(:,2))
         
         !--- CABLE radiation type forcings
         !--- kblum_rad% vars are computed in subroutine update_kblum_radiation 
         call um2cable_rr( kblum_rad%FBEAM(:,:,1), rad%fbeam(:,1))
         call um2cable_rr( kblum_rad%FBEAM(:,:,2), rad%fbeam(:,2))
         call um2cable_rr( kblum_rad%FBEAM(:,:,3), rad%fbeam(:,3))
      return
   end subroutine  um2cable_met_rad

!========================================================================
!========================================================================
!========================================================================

   !--- UM met forcing vars needed by CABLE commonly have UM dimensions
   !---(row_length,rows), which is no good to CABLE. These have to be 
   !--- re-packed in a single vector of active tiles. Hence we use 
   !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
   !--- if the land point is/has an active tile
   subroutine um2cable_rr(umvar,cablevar)
         use define_dimensions, only : mp
         use cable_um_tech_mod, only :um1
      implicit none
    
      real, intent(in), dimension(um1%row_length, um1%rows) :: umvar   
      real, intent(inout), dimension(mp) :: cablevar
      real, dimension(um1%land_pts,um1%ntiles) :: fvar   
      integer :: n,k,l,j,i

         fvar = 0.0
         DO N=1,um1%NTILES                     
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               J=(um1%LAND_INDEX(L)-1)/um1%row_length + 1
               I = um1%LAND_INDEX(L) - (J-1)*um1%row_length
               fvar(L,N) = umvar(I,J)
            ENDDO
         ENDDO
         cablevar =  pack(fvar,um1%l_tile_pts)
      return
   end subroutine um2cable_rr


!========================================================================
!========================================================================
!========================================================================

   !--- UM met forcing vars needed by CABLE which have UM dimensions
   !---(land_points)[_lp], which is no good to cable. These have to be 
   !--- re-packed in a single vector of active tiles. Hence we use 
   !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
   !--- if the land point is/has an active tile
   subroutine um2cable_lp(umvar, defaultin, cablevar, soiltype, skip )
         use define_dimensions, only : mp
         use cable_um_tech_mod, only :um1
      implicit none
      real, intent(in), dimension(um1%land_pts) :: umvar
      real, intent(in), dimension(10) :: defaultin    
      real, intent(inout), dimension(mp) :: cablevar
      integer, intent(inout), dimension(mp) :: soiltype
      real, dimension(:,:), allocatable:: fvar   
      logical, optional :: skip
      integer :: n,k,l,i

         
         allocate( fvar(um1%land_pts,um1%ntiles) )
         fvar = 0.0

         DO N=1,um1%NTILES
            DO K=1,um1%TILE_PTS(N)
               L = um1%TILE_INDEX(K,N)
               fvar(L,N) = umvar(L)
               if(.not. present(skip) ) then
                  if( N == um1%ntiles ) then
                     fvar(L,N) =  defaultin(9)
                  endif
               endif
            ENDDO
         ENDDO
     
         cablevar     =  pack(fvar,um1%l_tile_pts)
     
         if(.not. present(skip) ) then
            do i=1,mp
               if(soiltype(i)==9) cablevar(i) =  defaultin(9)         
            enddo        
         endif
      
         deallocate(fvar)

         return
   end subroutine um2cable_lp
 
!========================================================================
!========================================================================
!========================================================================

subroutine init_bgc_vars() 
   use define_types
   use define_dimensions, only : ncs, ncp
   use cable_um_tech_mod, only : bgc, veg
   use cable_common_module, only : vegin
   implicit none
   integer :: k

   ! note that ratecp and ratecs are the same for all veg at the moment. (BP)
    do k=1,ncp
       bgc%cplant(:,k) = vegin%cplant(k,veg%iveg)
       bgc%ratecp(k) = vegin%ratecp(k,1)
    enddo
    do k=1,ncs
      bgc%csoil(:,k) = vegin%csoil(k,veg%iveg)
      bgc%ratecs(k) = vegin%ratecs(k,1)
    enddo
return
end subroutine init_bgc_vars

!========================================================================
!========================================================================
!========================================================================

   subroutine init_sumflux_zero() 
      use cable_um_tech_mod, only : sum_flux
      implicit none
         sum_flux%sumpn = 0.; sum_flux%sumrp = 0.; sum_flux%sumrpw = 0.
         sum_flux%sumrpr = 0.; sum_flux%sumrs = 0.; sum_flux%sumrd = 0.
         sum_flux%dsumpn = 0.; sum_flux%dsumrp = 0.; sum_flux%dsumrs = 0.
         sum_flux%dsumrd = 0.; sum_flux%sumxrp = 0.;  sum_flux%sumxrs = 0.
      return
   end subroutine init_sumflux_zero 

!========================================================================
!========================================================================
!========================================================================

   subroutine alloc_cable_types()
      use define_dimensions, only : mp
      use define_types, only : alloc_cbm_var
      use cable_um_tech_mod, only : air, canopy, met, bal, rad, rough, &
                              soil, ssoil, sum_flux, veg, bgc
      implicit none 


         call alloc_cbm_var(air, mp)
         call alloc_cbm_var(canopy, mp)
         call alloc_cbm_var(met, mp)
         call alloc_cbm_var(bal, mp)
         call alloc_cbm_var(rad, mp)
         call alloc_cbm_var(rough, mp)
         call alloc_cbm_var(soil, mp)
         call alloc_cbm_var(ssoil, mp)
         call alloc_cbm_var(sum_flux, mp)
         call alloc_cbm_var(veg, mp)
         call alloc_cbm_var(bgc, mp)

      return
   end subroutine alloc_cable_types

!========================================================================
!========================================================================
!========================================================================


end module cable_um_init_subrs




