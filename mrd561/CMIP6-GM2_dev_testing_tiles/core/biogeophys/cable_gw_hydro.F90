!==============================================================================
! This source code is part of the 
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CABLE Academic User Licence Agreement 
! (the "Licence").
! You may not use this file except in compliance with the Licence.
! A copy of the Licence and registration form can be obtained from 
! http://www.accessimulator.org.au/cable
! You need to register and read the Licence agreement before use.
! Please contact cable_help@nf.nci.org.au for any questions on 
! registration and the Licence.
!
! Unless required by applicable law or agreed to in writing, 
! software distributed under the Licence is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the Licence for the specific language governing permissions and 
! limitations under the Licence.
! ==============================================================================
!
! Purpose: Soil hydrology, including conceptual representations of subgrid scale
! soil moisture, runoff, and flux processes
!
!
! Contact: Mark Decker
!          m.decker@unsw.edu.au
!          mdatmosci@gmail.com
!History:  A long time ago, in a galaxy far far away the following advancements
!  where developed utilizing CLM
!   1) GWsoilfreeze: Decker and Zeng 2006
!       soil ice-moisture-temperature function derived from in situ data
!        Empirically derived. Theoretical formulations (assumptions include
!        liquid and ice vapor pressure equilibrium, ignore contact anlges, pore
!        sizes, liquid, and ice, hysteresis, etc.
!        The empirical method matches nearly sat and very dry sites (out of
!        sample) better than other methods
!    2) smoistgw, aquifer_recharge Decker and Zeng 2010, Zeng and Decker 2009
!       finite volume based, non iterative soln to 1d richards eqn applicable
!       in hetero veriablly sat soils
!    3) ovrlandflux, saturated_fraction, subsurface_drainage Decker 2015
!       use assumed pdf to find sat/unsat fractions, subgrid topo to 
!       parameterize subgrid runoff/soil moisture processes
!    4) Paramters (GWdz, drain_dens, GWhyds_vec) Decker et al. 2018 (hopefully)
!
!   All snow processes use the default soilsnow subroutines
!  Soil thermal fluxes use GWstempv in soilsnow
!
!  OPTIONS:
!   gw_params%ssgw_ice_switch: use liquid content to diagnose the soil matric
!   potential, and alternative conductivity function 
!      physically more realistic, as freezing (thawing) oil has been deomonstrated to
!      behave similarly to drying (wetting) soils
!      seems ok now but not rigorously tested
!
!   gw_params%subsurface_sat_drainage:
!      TRUE: remove base flow (subsurface lateral flow) from layers that
!      are below the water table depth.  rate exp function of wtd
!      the efolding wtd depth depends on user parameter andA
!       drainage density
!      exponential function can be derived using topmodel approach and exp
!      distribution of topo index.

!      FALSE: Full qhlev for liq > sfc
!      else multiply by liq like true
!      but qhlev is wbliq2b+3*hyds/(sum)
!     
!
! ==============================================================================

MODULE cable_gw_hydro_module
   
   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp           

   USE cable_common_module, ONLY : gw_params,cable_user,&
                                   cable_runtime,&
                                   max_glacier_snowd,ktau_gl

   USE cable_soil_snow_module, ONLY : trimb, snow_processes_soil_thermal

   USE cable_data_module, only: issnow_type,point2constants


   IMPLICIT NONE

   PRIVATE

   TYPE(issnow_type), SAVE :: C

   !mrd561 GW params
   REAL(r_2), SAVE :: smp_cor = 8.0
   REAL(r_2), PARAMETER :: sucmin       = -1.0e8      ! minimum soil pressure head [mm]
   REAL(r_2), PARAMETER :: volwatmin    = 0.9          !min soil water coeff (volwatmin*watr)
   REAL(r_2), PARAMETER :: wtd_uncert   = 10.0         ! uncertaintiy in wtd calcultations [mm]
   REAL(r_2), PARAMETER :: wtd_max      = 75000.0  ! maximum wtd [mm]
   REAL(r_2), PARAMETER :: wtd_min      = 0.1      ! minimum wtd [mm]
   REAL(r_2), PARAMETER :: close_to_one = 0.9999999999
   REAL(r_2), PARAMETER :: m_to_mm = 1000.0
   REAL(r_2), PARAMETER :: mm_to_m = 0.001
   REAL(r_2), PARAMETER :: den_rat = 0.921
   REAL(r_2), PARAMETER :: zero = 0.0

  INTEGER, PARAMETER :: wtd_iter_max = 20! maximum number of iterations to find the water table depth                    

  ! ! This module contains the following subroutines that
  !are called from other modules
   PUBLIC :: soil_snow_gw,calc_srf_wet_fraction,sli_hydrology,& 
             pore_space_relative_humidity,set_unsed_gw_vars

CONTAINS


! -----------------------------------------------------------------------------

SUBROUTINE GWsoilfreeze(dels, soil, ssnow)
  !NOTE: this is only included because gw_model uses parameters XXX_vec
  !these are r_2.  this breaks bitwise compatibility with trunk
  !if acceptable this routine does the same thing but with r_2 soil params
  ! if max_ice_frac always set to frozen_limit and tgg_tmp is always C%TFRZ

   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
   REAL     , DIMENSION(mp,ms) :: tgg_old,tgg_new,tgg_tmp !tgg_old is previous point of when crosses freezing
   REAL(r_2), DIMENSION(mp)           :: sicefreeze
   REAL(r_2), DIMENSION(mp)           :: sicemelt
   REAL(r_2), DIMENSION(mp,ms)        :: wbice_delta,avail_por,delta_ice_vol
   REAL(r_2), DIMENSION(mp)           :: ice_mass,liq_mass,tot_mass
   INTEGER :: i,j,k
   REAL(r_2) :: func,funcderv,Aconst,Dconst,t_zero,t_one,dtmp
   REAL, DIMENSION(mp,ms) :: gammzz_snow
   REAL(r_2),DIMENSION(mp,ms) :: xx,max_ice_frac,iceF,den_css  !Decker and Zeng 20.9
   REAL(r_2) :: delta_wbliq,tmp_var

   call point2constants( C )
 
   max_ice_frac(:,:) = zero
   delta_ice_vol(:,:) = zero
   tgg_old(:,:) = ssnow%otgg(:,:)
   tgg_new(:,:) = ssnow%tgg(:,:)
   tgg_tmp(:,:) = tgg_old(:,:)

   gammzz_snow(:,:) = zero
   k=1
   do i=1,mp
      if (ssnow%isflag(i) .eq. zero.and. soil%isoilm(i) .ne. 9) then
           gammzz_snow(i,k) = real(C%cgsnow,r_2) * real(ssnow%snowd(i),r_2)
      end if
   end do

   do k=1,ms
   do i=1,mp

      ssnow%wmice(i,k)  = ssnow%wbice(i,k)*soil%zse_vec(i,k)*real(C%density_ice,r_2)
      ssnow%wmliq(i,k)  = ssnow%wbliq(i,k)*soil%zse_vec(i,k)*real(C%density_liq,r_2)
      ssnow%wmtot(i,k)  = ssnow%wmice(i,k) + ssnow%wmliq(i,k)

      if  ((ssnow%tgg(i,k) .lt. C%TFRZ)  .and. &
           (ssnow%tgg(i,k) .lt. ssnow%otgg(i,k))) then

            ssnow%otgg(i,k) = min(ssnow%otgg(i,k),C%TFRZ)

            iceF(i,k) = max(0.3_r_2,min(0.7_r_2,ssnow%wbliq(i,k)/max(ssnow%wb(i,k),1.0-8)))

            tgg_tmp(i,k) = (1._r_2 - iceF(i,k))*ssnow%otgg(i,k) + &
                           iceF(i,k)*ssnow%tgg(i,k)

            Aconst = 2.0_r_2*( (max(0.2,ssnow%wb(i,k)/soil%ssat_vec(i,k))**2.0_r_2))
            Dconst = exp(1.  - max(0.2,ssnow%wb(i,k)/soil%ssat_vec(i,k)))

            if (tgg_tmp(i,k) .lt. C%TFRZ) then
            max_ice_frac(i,k) = max(zero,ssnow%wb(i,k)-soil%Watr(i,k))*&
                                 (1._r_2 - exp(Aconst*(tgg_tmp(i,k)-C%TFRZ)))/Dconst
            else
                max_ice_frac(i,k) = zero
            end if

            delta_ice_vol(i,k) = max_ice_frac(i,k) - ssnow%wbice(i,k)

            !check amount of water we have
            delta_ice_vol(i,k) = min((ssnow%wbliq(i,k)-soil%watr(i,k))/den_rat, max(zero, delta_ice_vol(i,k) ) )

            delta_ice_vol(i,k) = min(delta_ice_vol(i,k), &
                                     max(zero,(ssnow%otgg(i,k)-ssnow%tgg(i,k))*ssnow%gammzz(i,k)/C%HLF) )

      elseif ((ssnow%tgg(i,k) .gt. C%TFRZ) .and. &
              (ssnow%tgg(i,k) .gt. ssnow%otgg(i,k)) .and. &
              ssnow%wbice(i,k) .gt. zero ) then

            ssnow%otgg(i,k) = min(ssnow%otgg(i,k),C%TFRZ)

            tgg_tmp(i,k) = C%TFRZ

            Aconst = 2.0_r_2*( (max(0.2,ssnow%wb(i,k)/soil%ssat_vec(i,k))**2.0_r_2))
            Dconst = exp(1.  - max(0.2,ssnow%wb(i,k)/soil%ssat_vec(i,k)))

            if (tgg_tmp(i,k) .lt. C%TFRZ) then
            max_ice_frac(i,k) = max(zero,ssnow%wb(i,k)-soil%Watr(i,k))*&
                                 (1._r_2 - exp(Aconst*(tgg_tmp(i,k)-C%TFRZ)))/Dconst
            else
                max_ice_frac(i,k) = zero
            end if

            delta_ice_vol(i,k) = max(zero, ssnow%wbice(i,k) - max_ice_frac(i,k))

            delta_ice_vol(i,k) = min(delta_ice_vol(i,k), &
                                  max(zero, (ssnow%tgg(i,k)-ssnow%otgg(i,k)) * ssnow%gammzz(i,k) / C%HLF)  )
      
      endif
   end do
   end do

   DO k = 1, ms
   DO i=1,mp

     
      if (ssnow%tgg(i,k) .lt. C%TFRZ .and. ssnow%tgg(i,k).lt.ssnow%otgg(i,k) .and. &
            delta_ice_vol(i,k) .gt. zero) then

         sicefreeze(i) = delta_ice_vol(i,k)*soil%zse_vec(i,k)*real(C%density_ice,r_2)

         ssnow%wbice(i,k) = ssnow%wbice(i,k) +&
                                  sicefreeze(i)/soil%zse_vec(i,k)/real(C%density_ice,r_2)

         delta_wbliq = max(zero,min(ssnow%wbliq(i,k),delta_ice_vol(i,k)*real(C%density_ice/C%density_liq,r_2)))

         ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),  &
                              (1.0- soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k) &
                            + (ssnow%wbliq(i,k) - delta_wbliq) * REAL(C%cswat*C%density_liq,r_2)   &
                            + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
                            )*soil%zse_vec(i,k)  + gammzz_snow(i,k)


         ssnow%tgg(i,k) = ssnow%tgg(i,k) + real(sicefreeze(i) )&
                             * C%hlf / real(ssnow%gammzz(i,k) )

         ssnow%wmice(i,k) = ssnow%wbice(i,k)*den_rat*soil%zse_vec(i,k)*m_to_mm
         ssnow%wmliq(i,k) = ssnow%wmtot(i,k) - ssnow%wmice(i,k)

         ssnow%wbliq(i,k) = ssnow%wmliq(i,k) / (soil%zse_vec(i,k)*m_to_mm)
         ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)

     elseif (ssnow%tgg(i,k) .gt. C%TFRZ .and. ssnow%tgg(i,k) .gt. ssnow%otgg(i,k) .and.&
               delta_ice_vol(i,k) .gt. zero) then

         sicemelt(i) =delta_ice_vol(i,k) * soil%zse_vec(i,k)*real(C%density_ice,r_2) 
        
         ssnow%wbice(i,k) = ssnow%wbice(i,k)  -  delta_ice_vol(i,k) 

         delta_wbliq = delta_ice_vol(i,k)*real(C%density_ice/C%density_liq,r_2)

         ssnow%gammzz(i,k) =max(soil%heat_cap_lower_limit(i,k), &
                           (1.0 - soil%ssat_vec(i,k))*soil%css_vec(i,k)*soil%rhosoil_vec(i,k)&
              + (ssnow%wbliq(i,k)+delta_wbliq) * REAL(C%cswat*C%density_liq,r_2)   &
              + ssnow%wbice(i,k) * REAL(C%csice*C%density_ice,r_2)&
                             )*soil%zse_vec(i,k) + gammzz_snow(i,k)           

         ssnow%tgg(i,k) = ssnow%tgg(i,k) - real(sicemelt(i) )&
                             * C%hlf / real(ssnow%gammzz(i,k) )

         ssnow%wmice(i,k) = ssnow%wbice(i,k)*den_rat*soil%zse_vec(i,k)*m_to_mm
         ssnow%wmliq(i,k) = ssnow%wmtot(i,k) - ssnow%wmice(i,k)

         ssnow%wbliq(i,k) = ssnow%wmliq(i,k) / (soil%zse_vec(i,k)*m_to_mm)
         ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)

      END IF

    
   END DO
   END DO

END SUBROUTINE GWsoilfreeze
! -----------------------------------------------------------------------------
!
!! -----------------------------------------------------------------------------
!
SUBROUTINE remove_transGW(dels, soil, ssnow, canopy, veg)
  !NOTE: this is only included because gw_model uses parameters XXX_vec
  !these are r_2.  this breaks bitwise compatibility with trunk
  !if acceptable this routine does the same thing but with r_2 soil params
 
   ! Removes transpiration water from soil.
   REAL, INTENT(IN)                    :: dels ! integration time step (s)
   TYPE(canopy_type), INTENT(INOUT)         :: canopy
   TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow
   TYPE(soil_parameter_type), INTENT(INOUT)    :: soil
   TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
   REAL(r_2), DIMENSION(mp,0:ms+1) :: diff 
   REAL(r_2), DIMENSION(mp)      :: xx,xxd,evap_cur
   REAL(r_2), DIMENSION(mp,ms) :: zse_mp_mm
   INTEGER :: k,i

   do k=1,ms
      do i=1,mp
         zse_mp_mm(i,k)  = real(soil%zse_vec(i,k)*C%density_liq,r_2)
      end do
   end do

   IF (cable_user%FWSOIL_switch.ne.'Haverd2013') THEN
 
      xx(:) = zero
      xxd(:) = zero
      diff(:,:) = zero
   
      DO k = 1,ms  
   
         DO i=1,mp
   
            if (canopy%fevc(i) .gt. zero) then
   
               xx(i) = canopy%fevc(i) * dels / C%hl * veg%froot(i,k) + diff(i,k-1)
               diff(i,k) = max(zero,ssnow%wbliq(i,k)-soil%swilt_vec(i,k)) &
                          * zse_mp_mm(i,k)
               xxd(i) = xx(i) - diff(i,k)
   
               if (xxd(i) .gt. zero) then
                  ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - diff(i,k)/zse_mp_mm(i,k)
                  diff(i,k) = xxd(i)
               else
                  ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - xx(i)/zse_mp_mm(i,k)
                  diff(i,k) = zero
               end if
   
   
             end if  !fvec > 0
   
         END DO  !mp
      END DO     !ms

   ELSE

     WHERE (canopy%fevc .lt. zero)
        canopy%fevw = canopy%fevw+canopy%fevc
        canopy%fevc = zero
     END WHERE
     DO k = 1,ms 
        ssnow%wbliq(:,k) = ssnow%wbliq(:,k) - ssnow%evapfbl(:,k)/zse_mp_mm(:,k)
     ENDDO

  ENDIF

  do k=1,ms
     do i=1,mp
        ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*zse_mp_mm(i,k)!mass
        ssnow%wmtot(i,k) = ssnow%wmliq(i,k) + ssnow%wmice(i,k)  !mass
        ssnow%wb(i,k)    = ssnow%wbliq(i,k) + ssnow%wbice(i,k)  !volume
     end do
  end do


END SUBROUTINE remove_transGW
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
!!!!!!!!!!!!!!MD GW code from here on!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ssnow, soil,veg, canopy,sli_call )
  !compute the surface runoff
  !divided into runoff from sat and unsat
  !unsat runoff occurs due to ice and from exceeding max infiltration rate
  !Decker 2015
  USE cable_common_module, ONLY : gw_params,cable_user

  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy
    LOGICAL, INTENT(IN)       :: sli_call

    !local variables
    INTEGER                   :: nglacier ! zerooriginal, 1 off, 2 new Eva
    INTEGER                   :: k, i, j
    REAL, DIMENSION(mp)       :: rnof5
    REAL, DIMENSION(mp)       :: sgamm
    REAL, DIMENSION(mp)       :: smasstot
    REAL, DIMENSION(mp,0:3)   :: smelt1                   !snow melt
    REAL(r_2), DIMENSION(mp)  :: icef,efpor               !tmp vars, fraction of ice in gridcell
    REAL(r_2)                 :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)  :: satfrac_liqice,S       !saturated fraction of cell, wtd in m
    REAL(r_2)                 :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2)                 :: fice
    REAL(r_2)                 :: slopeSTDmm

   if (sli_call) then
      do i=1,mp
         if (canopy%through(i) .ge. canopy%through_sn(i)) then
           ssnow%fwtop(i)  = max((canopy%through(i)-canopy%through_sn(i))/dels , zero)             ! liq precip rate (m s-1)
         else
           ssnow%fwtop(i) = max(canopy%through(i), zero)
         end if
      end do
   end if
   !amount of ice in surface layer
   do i = 1,mp
      efpor(i) = max(0.01, soil%ssat_vec(i,1) - ssnow%wbice(i,1))
      icemass  = ssnow%wmice(i,1)
      liqmass  = ssnow%wmliq(i,1)
      totmass  = max(ssnow%wmtot(i,1),0.01)
      icef(i)     = max(zero,min(1._r_2,gw_params%IceBeta*icemass / totmass))
   end do

   !sat fraction assuming topo controlled subgrid soil moisture distribution
   !called from cable_canopy for srf wet fraction alrady
   !call saturated_fraction(ssnow,soil,veg)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
      fice  = min(max(fice,zero),1._r_2)
      satfrac_liqice(i)   = max(zero,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
   end do

   do i=1,mp
      tmpa = ssnow%wbliq(i,1) / efpor(i)
      tmpb = max( (tmpa-satfrac_liqice(i))/max(0.01,(1._r_2-satfrac_liqice(i))), zero)
      tmpa = -2._r_2*soil%bch_vec(i,1)*soil%sucs_vec(i,1)/soil%zse_vec(i,1)/m_to_mm
      qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hyds_vec(i,1)*exp(-gw_params%hkrz*(0.5*soil%zse_vec(i,1)-gw_params%zdepth))

      ssnow%rnof1(i) = satfrac_liqice(i) * ssnow%fwtop(i) + &
                         (1._r_2-satfrac_liqice(i))*max((ssnow%fwtop(i)-qinmax) , zero)

      ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

   end do  !mp

  !add back to the lakes to keep saturated instead of drying
  do i=1,mp
     if (veg%iveg(i) .eq. 16) then
        ssnow%fwtop(i) = ssnow%fwtop(i) + ssnow%rnof1(i)
        ssnow%rnof1(i) = zero
     end if
  end do
           
   !---  glacier formation
   rnof5= zero

   if (sli_call .or. cable_runtime%UM .or. cable_user%gw_model) then
      nglacier = 0
   else
     nglacier = 2
   end if

   IF (nglacier == 2) THEN
      smelt1 = 0.0
      WHERE( ssnow%snowd > max_glacier_snowd )

         rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

         !---- change local tg to account for energy - clearly not best method
         WHERE( ssnow%isflag == zero)
            smasstot = zero
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * C%hlf                    &
                             / REAL( ssnow%gammzz(:,1) )
            ssnow%snowd = ssnow%snowd - rnof5
         ELSEWHERE
            smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
         END WHERE

      END WHERE

      DO k = 1, 3
         
         WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > zero)
            sgamm = ssnow%ssdn(:,k) * C%cgsnow * ssnow%sdepth(:,k)
            smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                          0.2 * ssnow%smass(:,k) )
            ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
            ssnow%snowd = ssnow%snowd - smelt1(:,k)
         END WHERE
      
      END DO
   
      WHERE( ssnow%isflag > zero) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)

      ssnow%rnof1 = ssnow%rnof1 + rnof5/dels   !include this runoff in suface runoff term
   
   END IF

  END SUBROUTINE ovrlndflx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  SUBROUTINE iterative_wtd (ssnow, soil, veg, include_aquifer)
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
  LOGICAL, INTENT(IN)                       :: include_aquifer  !use GWwb or only wb to find wtd?
 
  !Local vars 
  REAL(r_2), DIMENSION(mp,ms)   :: tmp_def
  REAL(r_2), DIMENSION(mp)      :: temp
  REAL(r_2), DIMENSION(mp)      :: def,defc,total_depth_column
  REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
  REAL(r_2), DIMENSION(mp)      :: lam,Nsucs_vec  !inverse of C&H B,Nsucs_vec
  INTEGER :: k,i,wttd,jlp

  !make code cleaner define these here 
  lam(:)        = 1._r_2/soil%bch_vec(:,ms)                                !1 over C&H B
  Nsucs_vec(:)  = soil%sucs_vec(:,ms)                                !psi_saturated mm

  if (include_aquifer) then  !do we include the aquifer in the calculation of wtd?

     do i=1,mp
        total_depth_column(i) = soil%GWdz(i)*m_to_mm
        def(i) = max(zero,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*m_to_mm
     end do   

  else
     def(:) = zero
     total_depth_column(:) = zero
  end if

  !total depth of soil column
  do k=1,ms
     do i=1,mp
         total_depth_column(i) = total_depth_column(i) + soil%zse_vec(i,k)*m_to_mm
     end do
  end do
  
  !comute the total mass away from full saturation
  do k=1,ms
     do i=1,mp

       def(i) = def(i) +                                                           &
                max(zero,(soil%ssat_vec(i,k)-(ssnow%wbliq(i,k)+ssnow%wbice(i,k)))*soil%zse_vec(i,k)*m_to_mm)
      end do  !mp
  end do  !ms

  !find the deficit if the water table is at the bottom of the soil column
  do i=1,mp
     defc(i) = (soil%ssat_vec(i,ms))*(total_depth_column(i)+Nsucs_vec(i)/(1._r_2-lam(i))*            &
             (1._r_2-((Nsucs_vec(i)+total_depth_column(i))/Nsucs_vec(i))**(1._r_2-lam(i)))) 
     defc(i) = max(0.1_r_2,defc(i)) 
  end do

  !initial guess at wtd
  ssnow%wtd(:) = total_depth_column(:)*def(:)/defc(:)

 !use newtons method to solve for wtd, note this assumes homogenous column but
 !that is ok 
  do i=1,mp
    if ((soil%isoilm(i) .ne. 9) .and. (veg%iveg(i) .ne. 16)) then      

      if (defc(i) > def(i)) then                 !iterate tfor wtd

        jlp=0

        mainloop: DO

          tempa   = 1._r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-lam(i))
          derv    = (soil%ssat_vec(i,ms))*(tempa-tempb) + &
                                       soil%ssat_vec(i,ms)

          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa   = 1._r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-lam(i))
          deffunc = (soil%ssat_vec(i,ms))*(ssnow%wtd(i) +&
                           Nsucs_vec(i)/(1-lam(i))* &
                     (tempa-tempb)) - def(i)
          calc    = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

            ssnow%wtd(i) = calc
            EXIT mainloop

          ELSEIF (jlp .ge. wtd_iter_max) THEN

            EXIT mainloop

          ELSE

            jlp=jlp+1
            ssnow%wtd(i) = calc

          END IF

        END DO mainloop  !defc .gt. def

      elseif (defc(i) .lt. def(i)) then

        jlp=0

        mainloop2: DO

          tmpc     = Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i)
          tempa    = (abs(tmpc/Nsucs_vec(i)))**(-lam(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-lam(i))
          derv     = (soil%ssat_vec(i,ms))*(tempa-tempb)
          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa    = (abs((Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i))/Nsucs_vec(i)))**(1._r_2-lam(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-lam(i))
          deffunc  = (soil%ssat_vec(i,ms))*(total_depth_column(i) +&
                     Nsucs_vec(i)/(1._r_2-lam(i))*(tempa-tempb))-def(i)
          calc     = ssnow%wtd(i) - deffunc/derv

          IF ((abs(calc-ssnow%wtd(i))) .le. wtd_uncert) THEN

            ssnow%wtd(i) = calc
            EXIT mainloop2

          ELSEIF (jlp==wtd_iter_max) THEN

            EXIT mainloop2

          ELSE

            jlp=jlp+1
            ssnow%wtd(i) = calc

          END IF

        END DO mainloop2  !defc .lt. def

      else  !water table depth is exactly on bottom boundary

        ssnow%wtd(i) = total_depth_column(i)

      endif

    endif  !check veg and soils

  end do   !mp loop

  !limit wtd to be within a psecified range
  do i=1,mp
     if (veg%iveg(i) .ge. 16) ssnow%wtd(i) = wtd_min
     ssnow%wtd(i) = min(wtd_max,max(wtd_min,ssnow%wtd(i) ) )
  end do


  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,veg,canopy)
  ! solves the modified richards equation (Zeng and Decker 20.9) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)     :: veg
    TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables
    
    !Local variables.  
    REAL(r_2), DIMENSION(mp,ms+1)       :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: rt

    INTEGER                             :: k,kk,i
    REAL(r_2), DIMENSION(mp,ms)         :: eff_por,old_wb  !effective porosity (mm3/mm3),wb(mm3/mm3),mass (mm) of eff_por
    REAL(r_2), DIMENSION(mp,ms)         :: msliq,msice             !mass of the soil liquid and ice water    
    REAL(r_2), DIMENSION(mp)            :: den
    REAL(r_2), DIMENSION(mp)            :: dne
    REAL(r_2), DIMENSION(mp)            :: num
    REAL(r_2), DIMENSION(mp)            :: qin
    REAL(r_2), DIMENSION(mp)            :: qout
    REAL(r_2), DIMENSION(mp)            :: dqidw0
    REAL(r_2), DIMENSION(mp)            :: dqidw1
    REAL(r_2), DIMENSION(mp)            :: dqodw0
    REAL(r_2), DIMENSION(mp)            :: dqodw1,dqodw2
    REAL(r_2), DIMENSION(mp)            :: s1,s2,tmpi,temp0voleq1,tempi
    REAL(r_2), DIMENSION(mp,ms)            :: dzmm
    REAL(r_2), DIMENSION(mp,0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(mp,ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi,dmass,dvol,cmass
    REAL(r_2), DIMENSION(mp,ms+1)       :: del_wb
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    !make code cleaner define these here
    dzmm(:,:)     = m_to_mm * soil%zse_vec(:,:)
    GWdzmm(:)     = m_to_mm * soil%GWdz(:)

    zmm(:,1:ms-1) = 0.5*(dzmm(:,1:ms-1)+dzmm(:,2:ms))
    zmm(:,ms)     = 0.5*(dzmm(:,ms) + GWdzmm(:) )

    zaq(:)        = sum(dzmm(:,:),dim=2) + 0.5*GWdzmm(:)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    do k=1,ms
       do i=1,mp
          old_wb(i,k) = ssnow%wb(i,k)
          rt(i,k) = zero
          at(i,k) = zero
          bt(i,k) = zero
          ct(i,k) = zero
       end do
    end do
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    CALL subsurface_drainage(ssnow,soil,veg)

    k = 1     !top soil layer
    do i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = (zmm(i,k+1)-zmm(i,k))
       dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) =  zero
       bt(i,k) =  dzmm(i,k)/dels + dqodw1(i)
       ct(i,k) =  dqodw2(i)      
    end do
    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
          den(i)     = (zmm(i,k) - zmm(i,k-1))
          dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = (zmm(i,k+1)-zmm(i,k))
          dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i)
          at(i,k) = -dqidw0(i)
          bt(i,k) =  dzmm(i,k)/dels - dqidw1(i) + dqodw1(i)
          ct(i,k) =  dqodw2(i)
       end do
    end do
       
    k = ms   !Bottom soil layer
    do i=1,mp
       den(i)     = (zmm(i,k) - zmm(i,k-1))
       dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(i,k)
       dne(i)     = (ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = zero
       dqodw1(i)  = zero
       dqodw2(i)  = zero
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) = -dqidw0(i)
       bt(i,k) =  dzmm(i,k)/dels - dqidw1(i) + dqodw1(i)
       ct(i,k) =  dqodw2(i) 
    end do
       
    CALL aquifer_recharge(dels,ssnow,soil,veg)

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - ssnow%qhlev(i,k)*dels/dzmm(i,k)   !volutermic liquid
       end do
    end do

    do i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/dzmm(i,ms)
    end do
    do i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-ssnow%qhlev(i,ms+1))*dels/GWdzmm(i)
    end do

    !determine the available pore space
    !volumetric
    do k=1,ms
       do i=1,mp
          eff_por(i,k)  = max(zero, soil%ssat_vec(i,k) - ssnow%wbice(i,k) )
       end do
    end do

    do i=1,mp
       xsi = zero
       if (ssnow%GWwb(i) .gt. soil%GWssat_vec(i)) then
          dmass =max(zero,  (ssnow%GWwb(i) - soil%GWssat_vec(i))*GWdzmm(i))
          ssnow%GWwb(i) = soil%GWssat_vec(i)

          cmass = min(dmass,max(zero,(eff_por(i,ms)-ssnow%wbliq(i,ms))*dzmm(i,ms)))
          ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) + cmass/dzmm(i,ms)

          dmass = max(zero, dmass - cmass)
          ssnow%Qrecharge(i) = ssnow%Qrecharge(i) - cmass/dels
          xsi = dmass
       end if


       do k=ms,1,-1
          if (ssnow%wbliq(i,k) .gt. eff_por(i,k)) then
             !oversat due to ice expansion, move down
             dmass = (ssnow%wbliq(i,k) - eff_por(i,k))*dzmm(i,k)
             ssnow%wbliq(i,k) = eff_por(i,k)

             !first add as much as possible to lower layer
             if (k .lt. ms) then
                cmass = min(dmass, max(zero, (eff_por(i,k+1)-ssnow%wbliq(i,k+1))*dzmm(i,k+1) ))
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) + cmass/dzmm(i,k+1)
             else
                cmass = min(dmass, max(zero, (soil%GWssat_vec(i)-ssnow%GWwb(i))*GWdzmm(i) ))
                ssnow%GWwb(i) = ssnow%GWwb(i) + cmass/GWdzmm(i)
             end if

             dmass = max(zero, dmass - cmass)

             !if left over, add to layer above
             if (k .gt. 1 .and. dmass .gt. zero) then
                cmass = min(dmass, max(zero, (eff_por(i,k-1)-ssnow%wbliq(i,k-1))*dzmm(i,k-1) ))
                ssnow%wbliq(i,k-1) = ssnow%wbliq(i,k-1) + cmass/dzmm(i,k-1)
                dmass = max(zero, dmass - cmass)
             end if
             xsi = xsi + dmass
          end if
       end do
       !if any left over (xsi > zero) add bottom up
       do k=ms,1,-1
          if ( xsi .gt. zero ) then
             cmass = min(xsi, max(zero, (eff_por(i,k)-ssnow%wbliq(i,k))*dzmm(i,k) ))
             ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + cmass/dzmm(i,k)
             xsi = max(zero, xsi - cmass)
          end if
       end do
       !any left over xsi must be added to qhz.
       !one day lakes/reservoirs will help..
       if (xsi .gt. zero) then
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
          xsi = zero
       end if
     
       do k = 1,ms
          xsi = zero             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. volwatmin*soil%watr(i,k)) then
             xsi = (soil%Watr(i,k)*volwatmin - ssnow%wbliq(i,k))*dzmm(i,k)  !in mm
             ssnow%wbliq(i,k) = soil%watr(i,k)*volwatmin
             if (k .lt. ms) then
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(i,k+1)
             else
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
             end if
          end if
       end do  !ms loop
 
       if ( (ssnow%GWwb(i) .lt. volwatmin*soil%GWwatr(i)) ) then
          xsi = (soil%GWwatr(i)*volwatmin - ssnow%GWwb(i)) / GWdzmm(i)  !mm
          ssnow%GWwb(i) = soil%GWwatr(i)
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       end if

   end do

   do k=1,ms
      do i=1,mp
         
       !update mass variables
         ssnow%wmliq(i,k)      = ssnow%wbliq(i,k) * &
                                         soil%zse_vec(i,k)*C%density_liq
         ssnow%wmice(i,k)      = ssnow%wbice(i,k) * &
                                         soil%zse_vec(i,k)*C%density_ice
         ssnow%wb(i,k)         = ssnow%wbliq(i,k) + ssnow%wbice(i,k)
         ssnow%wmtot(i,k)      = ssnow%wmliq(i,k) + ssnow%wmice(i,k)
      end do
   end do

   do i=1,mp
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var 
   end do  !mp loop
       

 END SUBROUTINE smoistgw

SUBROUTINE soil_snow_gw(dels, soil, ssnow, canopy, met, bal, veg)
   !Main hydro driver subroutine
   REAL                     , INTENT(IN)     :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(canopy_type)        , INTENT(INOUT)  :: canopy
   TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg
   TYPE(met_type)           , INTENT(INOUT)  :: met ! all met forcing
   TYPE (balances_type)     , INTENT(INOUT)  :: bal

   INTEGER             :: k,i
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: GWwb_ic
   REAL, DIMENSION(mp) :: tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx
   REAL(r_2), DIMENSION(mp,ms) :: gammzz_snow
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =zero
   REAL(r_2) :: wb_lake_T, rnof2_T
   LOGICAL :: use_sli
   LOGICAL, SAVE :: first_gw_hydro_call = .true.

   use_sli = .false. 

   ktau = ktau +1 
   
   zsetot = sum(soil%zse) 
   ssnow%tggav = zero

   CALL point2constants( C )

   DO k = 1, ms

      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot

   END DO


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) ssnow%t_snwlr = 0.01

   ssnow%otgg = ssnow%tgg
   do i=1,mp
      ssnow%fwtop1(i) = zero
      ssnow%fwtop2(i) = zero
      ssnow%fwtop3(i) = zero
      ssnow%runoff(i) = zero ! initialise total runoff
      ssnow%rnof1(i) = zero ! initialise surface runoff
      ssnow%rnof2(i) = zero ! initialise deep drainage
      ssnow%smelt(i) = zero ! initialise snowmelt
      ssnow%dtmlt(i,:) = zero 
      ssnow%osnowd(i) = ssnow%snowd(i)
      ! Scaling  runoff to kg/m^2/s (mm/s) to match rest of the model
      ssnow%sinfil(i) = zero   
      ssnow%qhz(i) = zero
   end do

   IF (cable_user%soil_thermal_fix) THEN
      soil%heat_cap_lower_limit(:,:) = zero  !allow /zeroto show bugs
   ELSE
      soil%heat_cap_lower_limit(:,:) = soil%css_vec(:,:) * soil%rhosoil_vec(:,:)
   END IF

   IF( (.NOT.cable_user%cable_runtime_coupled ) .and. (first_gw_hydro_call)) THEN
   
         IF (cable_runtime%um) canopy%dgdtg = zero ! RML added um condition
                                                  ! after discussion with BP
            
            WHERE( ssnow%tgg(:,:) < C%TFRZ)  
               ssnow%wbice(:,:) = (ssnow%wb(:,:)-soil%watr(:,:))*&
                                  (1._r_2 - exp(2.0*((ssnow%wb(:,:)/soil%ssat_vec(:,:))**2.0*&
                                  (ssnow%tgg(:,:)-C%TFRZ))))/exp(1._r_2 - ssnow%wb(:,:)/soil%ssat_vec(:,:))
            ENDWHERE

         WHERE ( soil%isoilm .eq. 9 .and. ssnow%snowd .le. 0.1*max_glacier_snowd) 

            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0

         END WHERE
   END IF

   gammzz_snow(:,:) = zero
   do i=1,mp
      gammzz_snow(i,1) = real((1. - ssnow%isflag(i)) * C%cgsnow * ssnow%snowd(i),r_2)
   end do

   IF( first_gw_hydro_call ) THEN
     do k=1,ms
       do i=1,mp
            ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),&
                               (1._r_2 - soil%ssat_vec(i,k))*&
                               soil%css_vec(i,k) * soil%rhosoil_vec(i,k)  &
                & + ssnow%wbliq(i,k) * real(C%cswat*C%density_liq,r_2)  &
                & + ssnow%wbice(i,k) * real(C%csice*C%density_liq*0.9,r_2) )* &
                 soil%zse_vec(i,k) +   gammzz_snow(i,k)
    
       end do
     end do
   ENDIF  ! if(.NOT.cable_runtime_coupled) and first_gw_hydro_call

   !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
   !find the mass of ice and liq from the prognostic volumetric values
   do k=1,ms
      do i=1,mp
         ssnow%wbliq(i,k) = ssnow%wb(i,k) - ssnow%wbice(i,k)                     !liquid volume
         ssnow%wmice(i,k) = ssnow%wbice(i,k)*real(C%density_ice*soil%zse(k),r_2) !ice mass
         ssnow%wmliq(i,k) = ssnow%wbliq(i,k)*real(C%density_liq*soil%zse(k),r_2) !liquid mass
         ssnow%wmtot(i,k) = ssnow%wmice(i,k) + ssnow%wmliq(i,k)                  !liq+ice mass

         ssnow%wblf(i,k)   = max(ssnow%wbliq(i,k)/soil%ssat_vec(i,k),0.01)
         ssnow%wbfice(i,k) = max(ssnow%wbice(i,k)/soil%ssat_vec(i,k),zero)

      end do
   end do



   do i=1,mp
      !initial water in the soil column
      wbtot_ic(i)  = sum(ssnow%wbliq(i,:)*C%density_liq*soil%zse(:),1) + &
                     sum(ssnow%wbice(i,:)*C%density_ice*soil%zse(:),1) + &
                     ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq
                  
      GWwb_ic(i) = ssnow%GWwb(i)

   end do

   gammzz_snow(:,:) = zero
   do i=1,mp
      gammzz_snow(i,1) = real((1. - ssnow%isflag(i)) * C%cgsnow * ssnow%snowd(i),r_2)
   end do

   do k=1,ms
     do i=1,mp
          ssnow%gammzz(i,k) = max(soil%heat_cap_lower_limit(i,k),&
                             (1._r_2 - soil%ssat_vec(i,k))*&
                             soil%css_vec(i,k) * soil%rhosoil_vec(i,k)  &
              & + ssnow%wbliq(i,k) * real(C%cswat*C%density_liq,r_2)  &
              & + ssnow%wbice(i,k) * real(C%csice*C%density_ice,r_2) )* &
               soil%zse_vec(i,k) +   gammzz_snow(i,k)
   
     end do
   end do
   !improve hiding, call single soilsnow subroutine to do all the
   !snow processes and thermal soil calculations

   CALL snow_processes_soil_thermal(dels,ssnow,soil,veg,canopy,met,bal)

   !leave here for now, could move into soilsnow as well
   CALL remove_transGW(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  GWsoilfreeze(dels, soil, ssnow)

   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]   

   CALL iterative_wtd (ssnow, soil, veg, .true. )  

   CALL ovrlndflx (dels, ssnow, soil, veg, canopy,use_sli )         !surface runoff, incorporate ssnow%pudsto?

   ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/C%hl               !remove soil evap from throughfall

   CALL smoistgw (dels,ktau,ssnow,soil,veg,canopy)               !vertical soil moisture movement. 
  
   ! correction required for energy balance in online simulations
   IF( cable_runtime%um ) THEN 

      !cls package - rewritten for flexibility
      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      !canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%dfe_ddq * ssnow%ddq_dtg)
      canopy%fes_cor = ssnow%dtmlt(:,1)*ssnow%dfe_dtg

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor

      !REV_CORR associated changes to other energy balance terms
      !NB canopy%fns changed not rad%flws as the correction term needs to
      !pass through the canopy in entirety, not be partially absorbed
      IF (cable_user%L_REV_CORR) THEN 
         canopy%fns_cor = ssnow%dtmlt(:,1)*ssnow%dfn_dtg
         canopy%ga_cor = ssnow%dtmlt(:,1)*canopy%dgdtg

         canopy%fns = canopy%fns + canopy%fns_cor
         canopy%ga = canopy%ga + canopy%ga_cor

         canopy%fess = canopy%fess + canopy%fes_cor
      ENDIF
   ENDIF

   do i=1,mp
      ssnow%pudsto(i) = zero  !no puddle
      ssnow%smelt(i)  = ssnow%smelt(i)/dels    !change units to mm/s.  cable_driver then reverts back to mm
      ssnow%runoff(i) = (ssnow%rnof1(i) + ssnow%rnof2(i))!*dels  !cable_driver converts from mm/s to mm
                                                        !rnof1 and rnof2 are already in mm/s
      ! Set weighted soil/snow surface temperature
      ssnow%tss(i) =  (1-ssnow%isflag(i))*ssnow%tgg(i,1) + ssnow%isflag(i)*ssnow%tggsn(i,1)

      !total water mass at the end of the soilsnow_GW routine
      ssnow%wbtot(i)  = sum(ssnow%wbliq(i,:)*C%density_liq*soil%zse(:),dim=1) + &
                     sum(ssnow%wbice(i,:)*C%density_ice*soil%zse(:),dim=1) + &
                     ssnow%GWwb(i)*soil%GWdz(i)*C%density_liq
      !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
      del_wbtot(i)   = dels * (ssnow%sinfil(i) - ssnow%rnof2(i) - canopy%fevc(i) / C%hl)
      !set below to keep track of water imbalance within the GW module explicitly.  also must change cable_checks
      !ssnow%wbtot(i) = ssnow%wbtot(i)-wbtot_ic(i)

   end do

   first_gw_hydro_call=.false.

END SUBROUTINE soil_snow_gw

SUBROUTINE calc_equilibrium_water_content(ssnow,soil)
  !find layer mean soil moisture and potential at equilibrium with wtd
  !Decker and Zeng 2009 and Zeng and Decker 2009
  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)   :: zaq    ,& !node depth of the aquifer
                                  GWdzmm ,& !aquifer thickness at each tile
                                  GWzimm  !aquifer layer interface depth

    REAL(r_2), dimension(mp,ms) :: dzmm,& !layer thickness for single tile
                                   zmm    !node depths in mm

    REAL(r_2), dimension(mp,0:ms):: zimm !layer interface depth in mm  

    REAL(r_2), dimension(mp,ms) :: eff_ssat,& !ssat - watr
                                   lam        !1.0/bch
    REAL(r_2), dimension(mp)    :: eff_GWssat,& !ssat - watr aquifer
                                   GWlam  !1.0/bch aquifer

    REAL(r_2)   :: intgrl_one,intgrl_zed,&  !values of f(wb) integrated over
                                            !layer that is at hydrostat w/ respect to wtd
                   sat_depth,unsat_depth,&  !depth [mm] of sat and unsat parts of layer
                   rel_pot_eq,           &  !hydrostat potential at layer
                   voleq1,wbrat             !equilibrium wbliq partway, relative wb

    INTEGER :: k,i



    !make code cleaner define these here
    dzmm(:,:) = m_to_mm * soil%zse_vec(:,:)
    zimm(0,:) = zero
    do k=1,ms
       zimm(:,k) = zimm(:,k-1) + dzmm(:,k)
    end do 
    zmm(:,1:ms)  = zimm(:,0:ms-1) + 0.5_r_2*dzmm(:,1:ms)

    GWdzmm(:) = m_to_mm*soil%GWdz(:)
    GWzimm(:) = sum(dzmm(:,:),dim=2) + GWdzmm(:)
    zaq(:)    = GWzimm(:) - 0.5*GWdzmm(:)

    !inverse B param always called lam
    !actually hydrologists say b is inverse lam
    lam(:,:) =  1._r_2 / soil%bch_vec(:,:)
    GWlam(:) = 1._r_2/soil%GWbch_vec(:)
    !temp var
    eff_ssat(:,:) = soil%ssat_vec(:,:) - soil%watr(:,:)
    eff_GWssat(:) = soil%GWssat_vec(:) - soil%GWwatr(:)
    !equilibrium water content
    do k=1,ms
       do i=1,mp

          if ((ssnow%wtd(i) .le. zimm(i,k-1))) then         !fully saturated

             ssnow%wbeq(i,k) = soil%ssat_vec(i,k)
             
          elseif ((ssnow%wtd(i) .le. zimm(i,k)) .and. &
                  (ssnow%wtd(i) .gt. zimm(i,k-1))) then

             intgrl_one    = 1._r_2
             rel_pot_eq = (soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/soil%sucs_vec(i,k)
             intgrl_zed    = rel_pot_eq**(1._r_2-lam(i,k))


             sat_depth   = zimm(i,k)-ssnow%wtd(i)  !unsat bottom part
             unsat_depth = dzmm(i,k) - sat_depth !sat top part

             !layer part saturated so sat/unsat treated seperately
             !voleq1 is wbeq from integrating over unsat portion
             !total is unsat_depth*voleq1 + sat_depth*eff_ssat + watr
             voleq1 = -soil%sucs_vec(i,k)*eff_ssat(i,k)/(1._r_2-lam(i,k))/&
                       (unsat_depth)*(intgrl_one - intgrl_zed) 

             ssnow%wbeq(i,k) = (voleq1*unsat_depth +eff_ssat(i,k)*sat_depth)/&
                                   dzmm(i,k) + soil%watr(i,k)
          else

             rel_pot_eq = (soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k))/soil%sucs_vec(i,k)
             intgrl_one = rel_pot_eq**(1._r_2-lam(i,k))

             rel_pot_eq = (soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(i,k-1))/soil%sucs_vec(i,k)
             intgrl_zed = rel_pot_eq**(1._r_2-lam(i,k))   

             ssnow%wbeq(i,k) = -soil%sucs_vec(i,k)*eff_ssat(i,k)/(1._r_2 - lam(i,k))/&
                               dzmm(i,k)*(intgrl_one - intgrl_zed)+soil%watr(i,k)
          end if
          !limited by possible range (never actually needed) 
          ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),soil%watr(i,k)),soil%ssat_vec(i,k))

          !limit to have some moisture to prevent zq -> -inf
          wbrat = max((ssnow%wbeq(i,k) - soil%watr(i,k))/eff_ssat(i,k), 0.01_r_2)

          ssnow%zq(i,k) = max(-soil%sucs_vec(i,k)*wbrat**(-soil%bch_vec(i,k)),sucmin)
       end do  !mp
    end do  !ms

    ! 
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(i,ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          rel_pot_eq = (soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/soil%GWsucs_vec(i)
          intgrl_one = rel_pot_eq**(1._r_2 - GWlam(i))

          rel_pot_eq = (soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/soil%GWsucs_vec(i)
          intgrl_zed = rel_pot_eq**(1._r_2 - GWlam(i))

          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*eff_GWssat(i)/(1._r_2-GWlam(i))/GWdzmm(i)*&
                           (intgrl_one - intgrl_zed)  + soil%GWwatr(i)   

       else           

          rel_pot_eq = 1._r_2 !by definition
          intgrl_zed = 1._r_2

          rel_pot_eq = (soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(i,ms))/soil%GWsucs_vec(i)
          intgrl_zed = rel_pot_eq**(1._r_2 - GWlam(i))

          sat_depth   = GWzimm(i) - ssnow%wtd(i)
          unsat_depth = GWdzmm(i) - sat_depth

          !layer part saturated so sat/unsat treated seperately
          !voleq1 is wbeq from integrating over unsat portion
          !total is unsat_depth*voleq1 + sat_depth*eff_ssat + watr
          voleq1 = -soil%GWsucs_vec(i)*eff_GWssat(i)/(1._r_2-GWlam(i))/unsat_depth*&
                    (intgrl_one - intgrl_zed)

          ssnow%GWwbeq(i) =        (voleq1*unsat_depth + &
                             eff_GWssat(i)*sat_depth   )/GWdzmm(i)+soil%GWwatr(i)

       end if

       ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))

       wbrat = max(0.01_r_2, (ssnow%GWwbeq(i)-soil%GWwatr(i))/eff_GWssat(i))

       ssnow%GWzq(i) = max(-soil%GWsucs_vec(i)*wbrat**(-soil%GWbch_vec(i)), sucmin)
       
    end do

END SUBROUTINE calc_equilibrium_water_content

SUBROUTINE calc_srf_wet_fraction(ssnow,soil,met,veg)
  !compute limitation of soil evaporation (wetfac)
  !or_evap: resistances used so wetfac=1
  !gw_model: nonlinear function (cos) from Sakaguchi et al. 2010
  !else: default linear function
  IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE (met_type), INTENT(IN)       :: met
    TYPE (veg_parameter_type), INTENT(IN)    :: veg

    !local variables
    REAL(r_2), DIMENSION(mp) :: icef,satfrac_liqice,S
    REAL(r_2)                :: fice,xx
    REAL(r_2)                :: dzmm_one,liqmass,icemass,totmass
    INTEGER                  :: i,j,k
    REAL(r_2)                :: wb_unsat,wb_lin,funcval
    REAL(r_2)                :: derv,slopeSTDmm,func_step
    REAL(r_2)                :: wb_evap_threshold

    CALL point2constants( C )


    IF (cable_user%or_evap) THEN

       call saturated_fraction(ssnow,soil,veg)

       ssnow%wetfac(:) = 1.0

      do i=1,mp
         IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9
   
         IF ( veg%iveg(i) == 16 .and. met%tk(i) >= C%TFRZ + 5. )   &
                 ssnow%wetfac(i) = 1.0! lakes: hard-wired number to be removed
   
         IF( veg%iveg(i) == 16 .and. met%tk(i) < C%TFRZ + 5. )   &
                 ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
      end do

    ELSEIF (cable_user%gw_model) THEN

       call saturated_fraction(ssnow,soil,veg)

   
       do i = 1,mp
          dzmm_one  = m_to_mm * soil%zse_vec(i,1)
          icemass  = ssnow%wmice(i,1)
          liqmass  = ssnow%wmliq(i,1)
          totmass  = max(ssnow%wmtot(i,1),real(1e-2,r_2))
          icef(i)     = max(zero,min(1._r_2, gw_params%IceBeta*icemass / totmass))
      end do
   
   
      !srf frozen fraction.  should be based on topography
      do i = 1,mp
         fice = (exp(-gw_params%IceAlpha*(1._r_2-icef(i)))-&
                 exp(-gw_params%IceAlpha))/&
                 (1._r_2-exp(-gw_params%IceAlpha))
         fice = min(1._r_2,max(zero,fice))
   
         satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)
   
         wb_unsat = (ssnow%wbliq(i,1) -&
                     ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
         wb_unsat = min(soil%ssat_vec(i,1),max(0.01,wb_unsat))
   
         wb_evap_threshold = min( max( &
                             gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), &
                             soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )
   
         !Sakguchi and Zeng 20.9
         if (wb_unsat .ge. wb_evap_threshold) then
            xx = 1.
         else
            xx = 0.25 * (1._r_2 - cos(3.14159_r_2*wb_unsat/(wb_evap_threshold)))**2.0
         end if
   
         ssnow%wetfac(i) = max(zero,min(1.0,satfrac_liqice(i) +&
                                        (1. - satfrac_liqice(i))*xx ) )

      end do

   ELSE  !Default formulation

       !call saturated_fraction(ssnow,soil,veg)
       ssnow%satfrac(:) = 1.0e-8
       ssnow%rh_srf(:)  = 1.0

       ssnow%wetfac = MAX( 1.e-6, MIN( 1.0,&
            ( REAL (ssnow%wb(:,1) ) - soil%swilt/ 2.0)                  &
            / ( soil%sfc - soil%swilt/2.0) ) )
   
       DO i=1,mp
   
          IF( ssnow%wbice(i,1) > zero )&
               ssnow%wetfac(i) = ssnow%wetfac(i) * &
                                real(MAX( 0.5_r_2, 1._r_2 - MIN( 0.2_r_2, &
                                 ( ssnow%wbice(i,1) / ssnow%wb(i,1) )**2 ) ) )
   
          IF( ssnow%snowd(i) > 0.1) ssnow%wetfac(i) = 0.9
   
          IF ( veg%iveg(i) == 16 .and. met%tk(i) >= C%tfrz + 5. )   &
               ssnow%wetfac(i) = 1.0! lakes: hard-wired number to be removed
   
          IF( veg%iveg(i) == 16 .and. met%tk(i) < C%tfrz + 5. )   &
               ssnow%wetfac(i) = 0.7 ! lakes: hard-wired number to be removed
   
       ENDDO
       ! owetfac introduced to reduce sharp changes in dry regions,
       ! especially in offline runs in which there may be discrepancies b/n
       ! timing of precip and temperature change (EAK apr20.9)
       ssnow%wetfac = 0.5*(ssnow%wetfac + ssnow%owetfac)

   ENDIF  !or_evap, gw_model, or default wetfac parameterization

END SUBROUTINE calc_srf_wet_fraction

SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
   !diagnose the hydraulic conductivity between layers (hk) and 
   !soil matric potential (smp) and the derivatives
   !(dhkdw) (dsmpdw)
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil 
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg

   INTEGER :: i,k,kk

   REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
                               s2, &
                               s_mid, &
                               liq_ratio, &
                               Dliq_ratio_Dz

   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers
   REAL(r_2), dimension(mp,ms+1) ::wb_temp 
   REAL(r_2), DIMENSION(mp,ms+1) :: hk_ice_factor
   REAL(r_2), DIMENSION(mp,ms) :: tot_pore_tmp

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

    do k=1,ms
       do i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01))
          ssnow%fracice(i,k) = (exp(-gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))&
                               -exp(-gw_params%IceAlpha))/(1._r_2-exp(-gw_params%IceAlpha))
       end do
    end do



    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), zero)

    if (gw_params%ssgw_ice_switch) then
       wb_temp(:,1:ms) =  ssnow%wbliq(:,:)

       liq_ratio(:) = 1.0
       k = ms 
       do i=1,mp
          liq_ratio(i) =min(1.,max(zero,wb_temp(i,k)/max(ssnow%wb(i,k),1e-8) ) )
       end do
       !aquifer ice 
       wb_temp(:,ms+1) = liq_ratio(:) * ssnow%GWwb(:)

       do k=1,ms
          kk = min(k+1,ms)
          do i=1,mp
             if (soil%isoilm(i) .eq. 9) then
                hk_ice_factor(i,k) = 10.0**(-gw_params%ice_impedence)
             else
                hk_ice_factor(i,k) = 10.0**(-gw_params%ice_impedence* &
                                         ( 0.5*(ssnow%wbice(i,k)/max(1.0-8,ssnow%wb(i,k)) + &
                                                ssnow%wbice(i,kk)/max(1.0-8,ssnow%wb(i,kk))) ) &
                                         )
             end if
          end do
       end do
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    else
       wb_temp(:,1:ms) = ssnow%wb(:,:)
       wb_temp(:,ms+1) = ssnow%GWwb(:)
       do k=1,ms
          kk = min(k+1,ms)
          do i=1,mp
             hk_ice_factor(i,k) = (1.0- 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))
          end do
       end do
       hk_ice_factor(:,ms+1) = hk_ice_factor(:,ms)
    end if


    !potential from soil water rention function
    !defined as layer average
    do k=1,ms 
       do i=1,mp
          s_mid(i) = (wb_temp(i,k)-soil%watr(i,k))/&  
              (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.01),1._r_2)

          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

          ssnow%smp(i,k) = min(ssnow%smp(i,k),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/s_mid(i)
       end do   
    end do

    !Aquifer potential
    do i=1,mp
       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/&
                    (soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.01),1._r_2)
       s2(i)    = soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%GWbch_vec(i)+2._r_2)

       ssnow%GWhk(i)     =s_mid(i)*s2(i) * hk_ice_factor(i,ms+1)
       ssnow%GWdhkdw(i)  =  (2._r_2*soil%GWbch_vec(i)+3._r_2)*&
                           s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i)) *&
                           hk_ice_factor(i,ms+1)

       s_mid(i) = (wb_temp(i,ms+1)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.01),1._r_2)

       ssnow%GWsmp(i)    = -soil%GWsucs_vec(i)*s_mid(i)**(-soil%GWbch_vec(i))
       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/s_mid(i)
    end do

    !hydraulic conductivity
    !Interfacial so uses layer i and i+1
    do k=1,ms
       kk=min(ms,k+1)
       do i=1,mp

          if (k .lt. ms) then
          tot_pore_tmp(i,k) = 0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k)+&
                                      soil%ssat_vec(i,kk)-soil%watr(i,kk))
          s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                           (wb_temp(i,kk)-soil%watr(i,kk))) / &
                           tot_pore_tmp(i,k)
          else
          tot_pore_tmp(i,k) = 0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k)+&
                                      soil%GWssat_vec(i)-soil%GWwatr(i))
          s1(i) = 0.5_r_2*((wb_temp(i,k)-soil%watr(i,k)) + &
                           (wb_temp(i,kk)-soil%GWwatr(i))) / &
                           tot_pore_tmp(i,k)
          end if
          s1(i) = min(max(s1(i),0.01),1._r_2)
          s2(i) = hk_ice_factor(i,k)*soil%hyds_vec(i,k)*&
                      s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)

          ssnow%hk(i,k)    =  s1(i)*s2(i)
          ssnow%dhkdw(i,k) = (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)/&
                            tot_pore_tmp(i,k)
       end do
    end do

    !conductivity between soil column and aquifer
    k = ms 
    kk = ms+1
       do i=1,mp

          s1(i) = 0.5_r_2*(max(wb_temp(i,ms )-soil%watr(i,ms),zero) + &
                           max(wb_temp(i,ms+1)-soil%GWwatr(i),zero)) / &
                           tot_pore_tmp(i,ms)

          s1(i) = min(max(s1(i),0.01),1._r_2)
          s2(i) = hk_ice_factor(i,ms)*soil%hyds_vec(i,ms)*&
                     s1(i)**(2._r_2*soil%bch_vec(i,ms)+2._r_2)

          ssnow%hk(i,ms)    = s1(i)*s2(i)
          ssnow%dhkdw(i,ms) = (2._r_2*soil%bch_vec(i,ms)+3._r_2)*s2(i)/tot_pore_tmp(i,ms)
       end do
 

END SUBROUTINE calc_soil_hydraulic_props


  SUBROUTINE aquifer_recharge(dt,ssnow,soil,veg)
  ! compute the flux between the soil column and the aquifer using
  ! darcys law

  IMPLICIT NONE
    real, intent(in) :: dt 
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg

    !local vars
    REAL(r_2), DIMENSION(mp,ms) :: dzmm,zmm  !layer thickness, depths in mm
    REAL(r_2), DIMENSION(mp)    :: GWz       !aquifer depth in mm
    integer :: i,k

    dzmm(:,:) = m_to_mm*soil%zse_vec(:,:)
    zmm(:,1) = 0.5*dzmm(:,1)
    do k=2,ms
       zmm(:,k) = zmm(:,k-1) + 0.5*(dzmm(:,k)+dzmm(:,k-1))
    end do

    GWz(:) = sum(dzmm(:,:),dim=1) + m_to_mm*0.5*soil%GWdz(:)
    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
    do i=1,mp
       if ((ssnow%wtd(i) .le. sum(dzmm(i,:),dim=1)) .or. &
           (veg%iveg(i) .ge. 16) .or. &
           (soil%isoilm(i) .eq. 9))  then

          ssnow%Qrecharge(i) = zero
       else
          ssnow%Qrecharge(i) = -0.5*(ssnow%hk(i,ms)*ssnow%GWhk(i))*&
                               ((ssnow%GWsmp(i)-ssnow%smp(i,ms)) -&
                                (ssnow%GWzq(i)-ssnow%zq(i,ms))) / &
                                (min(ssnow%wtd(i),GWz(i)) - &
                                 (zmm(i,ms)-0.5*dzmm(i,ms)))
       end if
    end do


  END SUBROUTINE aquifer_recharge

  SUBROUTINE subsurface_drainage(ssnow,soil,veg)
  ! compute the base flow (subsurface horizontal drainage)
  ! Decker 2015
  ! based on topmodel exp decline (subsurface_sat_drainage = true)
  ! subsurface_sat_drainage = false add lin/nonlin dep on sm amount
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg

    !local variables
    REAL(r_2), dimension(mp)       :: sm_tot,&!sum var
                                      drain_dens_fac,&  !1/drain density (km) 
                                      GWwb_liq,&!aquifer liquid, derived using last layer ice/liq
                                      col_liq,max_col_liq,min_col_liq
    REAL(r_2), dimension(mp,ms)    :: ice_factor,&!ice impact on hk (deps on ssgw_ice_switch
                                      cum_dzmm  !cummulative depth
    REAL(r_2), dimension(mp,ms+1)  :: zmm_node  !node depths in mm
                                       !subsurface drainage
    INTEGER, dimension(mp)         :: k_drain  !layer number with wtd

    !tmp vars
    REAL(r_2) :: rel_ice_bot,dice_dz,dz_tot
    integer :: i,k,kk

    cum_dzmm(:,1) = m_to_mm*soil%zse_vec(:,1)
    do k=2,ms
       cum_dzmm(:,k) = cum_dzmm(:,k-1) + m_to_mm*soil%zse_vec(:,k)
    end do

    do i=1,mp
        !Find GWwb_liq using average dice/dz from lowest 2 layers and
        !interpolate
        dz_tot = soil%zse_vec(i,ms) + soil%zse_vec(i,ms-1)
        rel_ice_bot = ssnow%wbice(i,ms)/soil%ssat_vec(i,ms)
        dice_dz = (soil%zse_vec(i,ms  )/dz_tot*ssnow%wbice(i,ms)/soil%ssat_vec(i,ms)-&
                   soil%zse_vec(i,ms-1)/dz_tot*ssnow%wbice(i,ms-1)/soil%ssat_vec(i,ms-1))/&
                                                (0.5*dz_tot)
        GWwb_liq(i) = ssnow%GWwb(i) - soil%GWssat_vec(i) * &
                      (rel_ice_bot + dice_dz * 0.5*(soil%GWdz(i)+soil%zse_vec(i,ms)))

       drain_dens_fac(i) = min(5.0, max(0.1 , 0.001/soil%drain_dens(i)/2.0))

    end do

   do i=1,mp
      !Note: future revision will have interaction with river here. nned to
      !work on router and add river type cells
      ssnow%qhz(i)  = min(max(soil%slope(i),0.00001),0.9)*&
                   gw_params%MaxHorzDrainRate* &
                    exp(-mm_to_m*ssnow%wtd(i)/max(0.1,(drain_dens_fac(i)+&
                          gw_params%EfoldHorzDrainRate)))  

      !drain from sat layers
      !only used for subsurface_sat_drainage=true
      k_drain(i) = ms+1
      do k=ms,2,-1
         if (ssnow%wtd(i) .le. cum_dzmm(i,k)) then
            k_drain(i) = k + 1
         end if
      end do
      k_drain(i) = max(k_drain(i),3)
   end do


   if (gw_params%ssgw_ice_switch) then
      do k=1,ms
         do i=1,mp
            ice_factor(i,k) = (1.0 - ssnow%wbice(i,k))**3.0
         end do
      end do
   else
      do k=1,ms
         do i=1,mp
            ice_factor(i,k) = (1._r_2-ssnow%fracice(i,k))
         end do
      end do

   end if

   do i=1,mp

       if (gw_params%subsurface_sat_drainage) then
           !only layers below wtd, lin func of wbliq
           do k=k_drain(i),ms
              ssnow%qhlev(i,k) = ssnow%qhz(i)*&
                                max(zero,ssnow%wbliq(i,k)-soil%watr(i,k))
           end do

           ssnow%qhlev(i,ms+1) = max((ssnow%GWwb(i) - soil%GWwatr(i)),zero)*&
                                     ice_factor(i,ms)*ssnow%qhz(i)

       else  !second option
           !all layers
           !weighted by thickness
           !above sfc no reduction, below nonlinear
           sm_tot(i) = zero
           sm_tot(i) = soil%GWdz(i)+sum(soil%zse_vec(i,:),dim=1)
           
           do k=1,ms
              if (ssnow%wbliq(i,k) .ge. soil%sfc_vec(i,k)) then
              ssnow%qhlev(i,k) = ssnow%qhz(i)*soil%zse_vec(i,k)/sm_tot(i)

              else
              ssnow%qhlev(i,k) = ssnow%qhz(i)*soil%zse_vec(i,k)/sm_tot(i) * max(zero, &
                                ( (ssnow%wbliq(i,k)-soil%watr(i,k))/&
                                  (soil%sfc_vec(i,k)-soil%watr(i,k)))**gw_params%sub_exp )
              end if
           end do
           ssnow%qhlev(i,ms+1) = ssnow%qhz(i)*soil%GWdz(i)/sm_tot(i)*&
                                   max(zero,( (GWwb_liq(i) - soil%GWwatr(i))/&
                                   (soil%GWssat_vec(i)-soil%GWwatr(i)))**gw_params%GW_sub_exp)
       end if

    end do

    do i=1,mp
       !Keep "lakes" saturated forcing qhz = zero  runoff only from lakes
       !overflowing
       if (soil%isoilm(i) .eq. 9 .or. veg%iveg(i) .ge. 16) then
          ssnow%qhlev(i,:) = zero
       end if
       ssnow%qhz(i) = sum(ssnow%qhlev(i,:),dim=1)

    end do  


  END SUBROUTINE subsurface_drainage


  SUBROUTINE saturated_fraction(ssnow,soil,veg)
  ! compute the fraction of the tile that is saturated (satfrac)
  ! Decker 2015
  ! uses assumed pdf of subgrid soil moisture
  ! used for srf runoff and soil evap
  ! gw_params%level_for_satfrac determines the number of layers to use for
  ! diagnosing the mean soil moisture
  ! default of 6, but 2 enables large sat frac during extreme precip events

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters

    REAL(r_2), DIMENSION(mp) :: S
    REAL(r_2) :: slopeSTDmm
    INTEGER :: i,k

     S(:) = zero
     do k=1,gw_params%level_for_satfrac
       S(:) = S(:) + max(0.1,min(1.0, &
              (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/&
               max(0.01,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )*soil%zse_vec(:,k)
     end do
     S(:) = S(:)/sum(soil%zse(1:gw_params%level_for_satfrac),dim=1)
     !srf frozen fraction.  should be based on topography
      do i = 1,mp
         !Saturated fraction
          if (gw_params%MaxSatFraction .gt. 1e-7 .and. veg%iveg(i) .lt. 16) then 
             slopeSTDmm = sqrt(min(max(&
                           gw_params%MaxSatFraction*soil%slope_std(i),&
                           1e-5),1.0)) ! ensure some variability
             ssnow%satfrac(i)    = max(zero,min(0.99_r_2,&
                 !note UM wants stdzero, and erf is not included then
                                   1._r_2 - my_erf( slopeSTDmm / sqrt(2.0*S(i)) ) ) )  
          elseif (veg%iveg(i) .lt. 16) then
             ssnow%satfrac(i) = zero
          else
             ssnow%satfrac(i) = 0.975
          end if
      end do


  END SUBROUTINE saturated_fraction

  SUBROUTINE pore_space_relative_humidity(ssnow,soil,veg)
  ! compute the relative humidity in the soil pore space
  ! used for gw_model and or_evap (Decker 2015, Decker et al. 2017)

    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT)    :: soil  ! soil parameters
    TYPE(veg_parameter_type), INTENT(INOUT)      :: veg

    !local variables
    REAL(r_2), DIMENSION(mp) :: unsat_wb,unsat_smp
    INTEGER :: i

    CALL point2constants( C )


    do i=1,mp
       if (veg%iveg(i) .lt. 16 .and. soil%isoilm(i) .ne. 9 .and. &
           ssnow%snowd(i) .le. 1e-8 ) then

          unsat_wb(i) =  (ssnow%wb(i,1) - soil%ssat_vec(i,1)*&
                      min(0.95,max(zero,ssnow%satfrac(i))))/(1.0- min(0.95,max(zero,ssnow%satfrac(i))))

          unsat_wb(i) = max(soil%watr(i,1)+1e-2, min(soil%ssat_vec(i,1), unsat_wb(i) ) )

          unsat_smp(i) = sign(soil%sucs_vec(i,1),-1.0) * min(1.0, &
                         (max(0.01,(unsat_wb(i)-soil%watr(i,1))/(soil%ssat_vec(i,1)-&
                         soil%watr(i,1)) ) )** (-soil%bch_vec(i,1)) )


          unsat_smp(i) = max(-1.0e7,unsat_smp(i) )/mm_to_m

          ssnow%rh_srf(i) = max(zero,min(1., &
                         exp(C%grav*unsat_smp(i)/(ssnow%tgg(i,1)*461.4)) ) )
         
       else

          ssnow%rh_srf(i) = 1.0

       end if
    end do  


  END SUBROUTINE pore_space_relative_humidity

  SUBROUTINE sli_hydrology(dels,ssnow,soil,veg,canopy)
  ! driving routine when hydrology routines used with sli
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(INOUT)    :: veg  ! veg parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy

    LOGICAL, SAVE :: sli_call = .true.

    REAL(r_2), DIMENSION(ms) :: dzmm
    REAL(r_2), DIMENSION(mp) :: zmm
    REAL(r_2), DIMENSION(mp) :: zaq

    CALL point2constants( C )

    call iterative_wtd (ssnow, soil, veg, cable_user%test_new_gw)
 
    CALL calc_soil_hydraulic_props(ssnow,soil,veg)
 
    call  ovrlndflx (dels, ssnow, soil, veg,canopy,sli_call )
 
    CALL subsurface_drainage(ssnow,soil,veg)
 
    call aquifer_recharge(dels,ssnow,soil,veg)




  END SUBROUTINE sli_hydrology


  SUBROUTINE set_unsed_gw_vars(ssnow,soil,canopy)
  ! sets variables used in gw_model when gw_model not used
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT)    :: soil ! soil parameters
    TYPE (canopy_type), INTENT(INOUT)           :: canopy

       ssnow%qhlev = zero
       ssnow%Qrecharge = zero
       ssnow%fwtop = zero
       ssnow%wtd = zero
       ssnow%satfrac = 1.0
       ssnow%qhz = zero
       ssnow%wbliq = ssnow%wb - ssnow%wbice
       canopy%sublayer_dz = zero
       ssnow%rtevap_sat = zero
       ssnow%rtevap_unsat = zero

       ssnow%GWwb = 0.9*soil%ssat


  END SUBROUTINE set_unsed_gw_vars

  real(r_2) function my_erf(x)
  ! an errror function to be used when code must conform to 2003 standard
  !  2003 standard does not have intrinsic erf

  real(r_2), intent(in) :: x 
  real(r_2)             :: tmp_val, tmp

   tmp = 1.0 / ( 1.0 + 0.5 * abs(x) )

   tmp_val =  tmp * exp(-abs(x) * abs(x) - 1.26551223 + tmp *     &
             ( 1.00002368 + tmp * ( 0.37409196 + tmp *          &
         ( 0.09678418 + tmp * (-0.18628806 + tmp *              &
                     ( 0.27886807 + tmp * (-1.13520398 + tmp *          &
         ( 1.48851587 + tmp * (-0.82215223 + tmp * 0.17087277 )))))))))

  if ( x.lt.0.0 ) tmp_val = 2.0 - tmp_val

  my_erf = 1.0 - tmp_val


  if ( x.lt.zero ) tmp_val = 2.0- tmp_val

  my_erf = 1.0- tmp_val

  end function my_erf

END MODULE cable_gw_hydro_module
