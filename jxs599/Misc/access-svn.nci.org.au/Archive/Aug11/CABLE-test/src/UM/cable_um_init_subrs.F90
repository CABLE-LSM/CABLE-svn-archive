#include "include/cable_directives.h"

module cable_um_init_subrs
   implicit none

   contains
!========================================================================
!========================================================================
!========================================================================                                   
         
         

   subroutine initialize_soil( & 
#                    include "include/args/um_soil.h"
                  &  )
      use define_dimensions, only : ms
      use cable_variables, only : soil, veg, ssoil 
      use cable_um_tech_mod
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
#     include "include/decs/um_soil.h"
      !___defs 1st call to CABLE in this run
      logical, save :: first_call= .true.
      integer :: i,j,k,L,n
      real, allocatable:: tempvar(:)
      logical :: skip =.true. 

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_soil')
      
         if (first_call) then 
            !--- this is text which hardwiressome soil params.
#           include "include/cable_soil_params.h"
            
            ssoil%pudsto = 0.0
            ssoil%pudsmx = 0.0
         
            !--- defined in soil_thick.h in UM
            soil%zse = dzsoil
            soil%zshh(1)=0.5*soil%zse(1) ! distance between consecutive layer midpoints
            soil%zshh(ms+1)=0.5*soil%zse(ms)
            soil%zshh(2:ms) = 0.5 * (soil%zse(1:ms-1) + soil%zse(2:ms))
            
            !-------------------------------------------------------------------
            !--- UM met forcing vars needed by CABLE which have UM dimensions
            !---(land_pts,ntiles)[_lp], which is no good to cable. These have to be 
            !--- re-packed in a single vector of active tiles. Hence we use 
            !--- conditional "mask" l_tile_pts(land_pts,ntiles) which is .true.
            !--- if the land point is/has an active tile
            !--- generic format:
            !--- um2cable_lp( UM var, default value for snow tile, CABLE var, mask )
            !--- where mask tells um2cable_lp whether or not to use default value 
            !--- for snow tile 
            !-------------------------------------------------------------------
            
            call um2cable_lp( BEXP, soilin%bch, soil%bch )
            
            allocate( tempvar(um1%land_pts) )
            tempvar = soilin%sand(9)*0.3  + soilin%clay(9)*0.25 + &
                        soilin%silt(9)*0.265
            call um2cable_lp( HCON, tempvar, soil%cnsd)
            deallocate( tempvar )
            
            call um2cable_lp( satcon,SATCON*1000.0, soil%hyds)
            call um2cable_lp( sathh, soilin%sucs, soil%sucs)
            call um2cable_lp( smvcst, soilin%ssat, soil%ssat)
            call um2cable_lp( smvcwt, soilin%swilt, soil%swilt)
            call um2cable_lp( smvccl, soilin%sfc, soil%sfc)
            call um2cable_lp( albsoil, albsoil, soil%albsoil(:,1), skip )
            !-------------------------------------------------------------------
            
            !--- (re)set values for CABLE
            !-------------------------------------------------------------------
            soil%ibp2    =  NINT(soil%bch)+2
            soil%i2bp3   =  2*NINT(soil%bch)+3
            !--- satcon in UM is in mm/s; Cable needs m/s
            soil%hyds    =  soil%hyds / 1000.
            soil%sucs    =  abs( soil%sucs )
            !jhan:coupled runs 
            soil%sucs    =   soil%sucs + 0.1 
            soil%hsbh    =  soil%hyds*ABS(soil%sucs)*soil%bch
            
            where (soil%ssat > 0. ) soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
           
            !jhan***
            !   soil_type is only a SCM variable in UM, but in Cable is needed to
            !   initialise  some parameters. For permanent ice use Cable's soil type 9
            soil%isoilm  =  2
   
            ! this is incorrect and need to be changed
            where ( soil%albsoil(:,1) > 0.7 )   ! in UM soil%albsoil = 0.75 for permanet ice points!!!
               soil%isoilm = 9
            endwhere
            !--- these are temporary 
            soil%rhosoil =  soilin%rhosoil(soil%isoilm)
            soil%css     =  soilin%css(soil%isoilm)
         
            !--- for permanent ice use Cable's soil properties ********for ver 6.3 old soil data *****
            WHERE ( soil%isoilm == 9 )
               soil%bch     =  soilin%bch(soil%isoilm)    ! Cable's
               soil%ibp2    =  NINT(soil%bch)+2
               soil%i2bp3   =  2*NINT(soil%bch)+3
               soil%rhosoil =  soilin%rhosoil(soil%isoilm)
               soil%cnsd    =  soilin%sand(soil%isoilm)*0.3 + & ! Cable's
                               soilin%clay(soil%isoilm)*0.25 + &
                               soilin%silt(soil%isoilm)*0.265
               soil%css     =  soilin%css(soil%isoilm)
               soil%hyds    =  soilin%hyds(soil%isoilm)  ! Cable's
               soil%sucs    =  soilin%sucs(soil%isoilm)  ! Cable's
               soil%sucs    =  abs( soil%sucs )
               soil%sucs    =  max(0.106,soil%sucs )
               soil%hsbh    =  soil%hyds*ABS(soil%sucs)*soil%bch 
               soil%ssat    =  soilin%ssat(soil%isoilm)   ! Cable's
               soil%swilt   =  soilin%swilt(soil%isoilm)  ! Cable's
               soil%sfc     =  soilin%sfc(soil%isoilm)    ! Cable's
               soil%pwb_min =  (soil%swilt / soil%ssat )**soil%ibp2
            END WHERE
            
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
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
      real, intent(in), dimension(um1%land_pts, um1%npft) :: canht_ft, lai_ft 
      !___defs 1st call to CABLE in this run
      logical, save :: first_call= .true.

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_veg')
   
         if (first_call) & 
            call read_veg_pars(um1%ntiles) 
         
         !---clobbers veg height, lai and resets ivegt for CABLE tiles
         call clobber_height_lai( canht_ft, lai_ft )

         if (first_call) &
            call init_veg_pars_fr_vegin() 
         
         first_call= .false.
        
      return
   end subroutine initialize_veg

!========================================================================
!========================================================================
!========================================================================

   subroutine initialize_radiation( &
#                       include "include/args/um_rad.h"
                        , &
#                       include "include/args/um_met.h"
                      & )    
      use cable_variables, only : rad, soil, met 
      use physical_constants, only : tfrz
      use other_constants, only : RAD_THRESH
      use cable_um_tech_mod
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
#        include "include/decs/um_rad.h"
#        include "include/decs/um_met.h"
                
      !___defs 1st call to CABLE in this run. OK in UM & coupled
      logical, save :: first_call= .true.
     
         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_radiation')
         
         if (first_call) then
            rad%albedo_T = soil%albsoil(:,1)
            first_call = .false.
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
         call um2cable_rr( LS_RAIN*um1%TIMESTEP, met%precip)
         call um2cable_rr( LS_SNOW*um1%TIMESTEP, met%precip_sn)
         call um2cable_rr( TL_1, met%tk)
         call um2cable_rr( QW_1, met%qv)
         call um2cable_rr( VSHR_LAND, met%ua)
         call um2cable_rr( PSTAR*0.01, met%pmb)
      
         !---re-set some of CABLE's forcing variables
         met%precip   =  met%precip + met%precip_sn
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
      use cable_variables, only : canopy 
      use cable_um_tech_mod
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
 
   subroutine initialize_soilsnow( &  
#                    include "include/args/um_soil.h"
                     , &
#                    include "include/args/um_snow.h"                  
                   &, sin_theta_latitude ) 

      use cable_variables, only :soil, ssoil, met, bal, veg
      use define_dimensions, only : mp
      use physical_constants, only : tfrz
      use cable_um_tech_mod, only : um1
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none

#     include "include/decs/um_soil.h"
#     include "include/decs/um_snow.h"
      real, intent(in), dimension(um1%row_length, um1%rows) :: sin_theta_latitude

      integer :: i,j,k,L,n
      real  :: zsetot
      real, allocatable:: fwork(:,:,:), sfact(:), fvar(:)
      logical :: skip =.true. 
      logical :: first_call = .true.

         if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
            call cable_stat('initialize_soilsnow')
         
         ssoil%wbtot1 = 0
         ssoil%wbtot2 = 0
         snow_tile = min(50000., snow_tile)

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

            !jhan:things which should not need to be resone if restart is matured
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
            call um2cable_lp( fland, fland, ssoil%fland, skip )
            call um2cable_lp( fvar, fvar, ssoil%ifland, skip )
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
          
   subroutine initialize_roughness( & 
#               include "include/args/um_met.h"
                & , htveg )  
      use cable_variables, only : rough, veg
      use cable_um_tech_mod
      use cable_common_module, only : ktau_gl
      use define_dimensions, only : mp
      use cable_diag_module, only : cable_stat
      use cable_common_module, only : cable_runtime, cable_user
      implicit none
         
#     include "include/decs/um_met.h"
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
         
         rough%hruff= MAX(0.01,veg%hc)
         rough%hruff_grmx = pack(jHRUFF, um1%l_tile_pts) 

         deallocate(jhruff ) 
         deallocate(jhwork) 
   
      return
   end subroutine initialize_roughness

!========================================================================
!========================================================================
!========================================================================

   subroutine read_veg_pars(ntiles) 
      use cable_common_module
      use io_variables, only : filename 
      use define_dimensions
      use cable_um_tech_mod
      implicit none
   
      integer, intent(in) :: ntiles
      integer :: ioerror,a
      character(len=999) :: comments
      character(len=99) :: vegtypetmp, vegnametmp 
      integer :: jveg 
      real :: notused
      character*200 :: paramfile

          paramfile=filename%VEG

          OPEN(40,FILE=trim(paramfile), &
                  STATUS='old',ACTION='READ',IOSTAT=ioerror)
            write(6,*)'CABLE_log:Opening veg params file: ',trim(paramfile)

         IF(ioerror/=0) then
            print *, 'CABLE_log:Cannot open vegetation type definitions.'
            stop
         endif
         ! assume using IGBP vegetation types
         READ(40,*) comments
         READ(40,*) mvtype
         DO a = 1,mvtype 
            READ(40,*) jveg, vegtypetmp, vegnametmp
            IF (jveg .GT. mvtype) then 
              print *,'CABLE_log:jveg out of range in paramgeter file'
              stop
            endif
            !jh:no veg_desc here
            !veg_desc(jveg) = vegnametmp
           
            ! BP adding the next conditional statement to move the 4 non-veg
            ! types to type number 14-17
            !rml test for 15 or 17 types in parameter file
            IF (jveg .GT. 11.and. mvtype.eq.15) jveg = jveg + 2
            !jh:no veg_desc here
            !veg_desc(jveg) = vegnametmp
            ! vegtype(jveg) = vegtypetmp     ! not yet used in offline mode
            READ(40,*) vegin%hc(jveg), vegin%xfang(jveg), notused,  &
                     &   vegin%dleaf(jveg), vegin%frac4(jveg)

            !kdcorbin, 09/10 - added refl and taul variables to read
            READ(40,*) vegin%reflin(1,jveg),vegin%reflin(2,jveg),notused, &
                       notused,notused,notused
            READ(40,*) vegin%taulin(1,jveg),vegin%taulin(2,jveg),notused, &
                       notused,notused,notused
            READ(40,*) notused, notused, notused, vegin%xalbnir(jveg)
            !jhan:why is this reset
            vegin%xalbnir(jveg) = 1.0
            
            READ(40,*) notused, notused, vegin%canst1(jveg), &
                     & vegin%shelrb(jveg), vegin%vegcf(jveg), vegin%extkn(jveg)
            READ(40,*) vegin%vcmax(jveg), vegin%rp20(jveg), &
                       vegin%rpcoef(jveg), &
                       vegin%rs20(jveg)
            READ(40,*) vegin%tminvj(jveg), vegin%tmaxvj(jveg), &
                       vegin%vbeta(jveg), &
            !jh: no & vegin%rootbeta(jveg)
                       & notused
            READ(40,*) vegin%cplant(1:3,jveg), vegin%csoil(1:2,jveg)
            ! rates not currently set to vary with veg type
            READ(40,*) vegin%ratecp(1:3,jveg), vegin%ratecs(1:2,jveg)
         END DO
         !jhan:Evauses this
         !jh:vegin%vcmax(2) = 1.2*vegin%vcmax(2)
         !jh:vegin%vcmax(1) = 0.5*vegin%vcmax(1)

            write(6,*)'CABLE_log:Closing veg params file: ',trim(paramfile)
         CLOSE(40)

         !jhan:dow e have    a fix for this
         ! EAK temporary fix
         vegin%tminvj = 0.0
         vegin%tmaxvj = 15.0

         ! kdcorbin, 08/10 - commenting fill loop as file now has 17 types
         ! BP changing the total number from 15 to 17
         ! jhan:check this because i think BP issed below push subr.
         ! BP added the next loop to fill in values for types 12 and 13
         !    in case some WHERE statements uses the vegin variables.
         ! rml add test for 15 or 17 types in file
         if (mvtype.eq.15) then
           mvtype = mvtype + 2
           DO jveg = 12, 13
             vegin%hc(jveg) = vegin%hc(5)
             vegin%xfang(jveg) = vegin%xfang(5)
             vegin%dleaf(jveg) = vegin%dleaf(5)
             vegin%frac4(jveg) = vegin%frac4(5)
             vegin%xalbnir(jveg) = vegin%xalbnir(5)
             vegin%extkn(jveg) = vegin%extkn(5)
             vegin%canst1(jveg) = vegin%canst1(5)
             vegin%shelrb(jveg) = vegin%shelrb(5)
             vegin%vegcf(jveg) = vegin%vegcf(5)  
             vegin%vcmax(jveg) = vegin%vcmax(5)
             vegin%rp20(jveg) = vegin%rp20(5)
             vegin%rpcoef(jveg) = vegin%rpcoef(5)
             vegin%rs20(jveg) = vegin%rs20(5)
             vegin%tminvj(jveg) = vegin%tminvj(5)
             vegin%tmaxvj(jveg) = vegin%tmaxvj(5)
             vegin%vbeta(jveg) = vegin%vbeta(5)
             vegin%cplant(1:3,jveg) = vegin%cplant(1:3,5)
             vegin%csoil(1:2,jveg) = vegin%csoil(1:2,5)
             vegin%ratecp(1:3,jveg) = vegin%ratecp(1:3,5)
             vegin%ratecs(1:2,jveg) = vegin%ratecs(1:2,5)
           END DO
         endif

         return

   end subroutine read_veg_pars 

!========================================================================
!========================================================================
!========================================================================
   
   subroutine clobber_height_lai( um_htveg, um_lai )
      use cable_variables, only : veg
      use cable_um_tech_mod, only : um1, kblum_veg
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
      use cable_um_tech_mod
      use cable_variables, only : veg, soil 
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
        
      subroutine update_kblum_radiation( sw_down, cos_zenith_angle, surf_down_sw )
         use cable_um_tech_mod!, only : um1, um_rad, kblum_rad
         implicit none
         real, intent(inout), dimension(um1%row_length, um1%rows) :: sw_down
         real, intent(in), dimension(um1%row_length, um1%rows) :: cos_zenith_angle
         real, intent(in), dimension(um1%row_length, um1%rows, 4) :: surf_down_sw 

            !jhan: do you really want to be changing sw_down            
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
      use cable_um_tech_mod, only : kblum_rad
      use cable_variables, only : met, rad
      use cable_um_tech_mod, only :um1
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
   subroutine um2cable_lp(umvar, defaultin, cablevar, skip )
         use define_dimensions, only : mp
         use cable_um_tech_mod, only :um1
      implicit none
      real, intent(in), dimension(um1%land_pts) :: umvar, defaultin    
      real, intent(inout), dimension(mp) :: cablevar
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
         deallocate(fvar)

         return
   end subroutine um2cable_lp
 
!========================================================================
!========================================================================
!========================================================================

subroutine init_bgc_vars() 
   use define_types
   use define_dimensions, only : ncs, ncp
   use cable_variables, only : bgc, veg
   use cable_um_tech_mod
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
      use cable_variables, only : sum_flux
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
      use define_types
      use cable_variables
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




