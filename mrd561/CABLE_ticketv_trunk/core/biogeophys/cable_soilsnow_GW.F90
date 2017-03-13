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
! Purpose: All routines for calculating soil temperature and moisture
!          and snow calculations
!
! Contact: Eva.Kowalczyk@csiro.au
!
! History: v2.0 Tighter water budget
!          v2.0 Hydraulic redistribution subroutine (with namelist switch). 
!               NB Currently hard-wired to veg types 2 and 7 
!                  (usually evergreen broadleaf and c4 grass)
!          v2.0 ssoil variable renamed ssnow
!          Mark Decker - used ssnow as base for ssgw.  Could be part of same module
!
! ==============================================================================

MODULE cable_soil_snow_gw_module
   
   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp           

   USE cable_data_module, ONLY : issnow_type, point2constants

   USE cable_common_module, ONLY : gw_params,cable_user

   USE cable_soil_snow_module, ONLY: remove_trans,snowdensity, snow_melting, &
                                     snowcheck, snowl_adjust,snow_accum, stempv,&
                                     trimb,soilfreeze


   IMPLICIT NONE

   PRIVATE

!   TYPE ( issnow_type ) :: C 
   
   REAL, PARAMETER ::                                                          &
      cgsnow = 2090.0,     & ! specific heat capacity for snow
      csice = 2.100e3,     & ! specific heat capacity for ice
      cswat = 4.218e3,     & ! specific heat capacity for water
      rhowat = 1000.0,     & ! density of water
      snmin = 1.,          & ! for 3-layer;
      max_ssdn = 750.0,    & !
      max_sconds = 2.51!,   & !
   !   frozen_limit = 0.85    ! EAK Feb2011 (could be 0.95) 
      
   !mrd561 GW params
   !Should read some in from namelist
   REAL(r_2), PARAMETER :: sucmin  = -1.0e8, &! minimum soil pressure head [mm]
                      volwatmin    = 1e-4,        &!min soil water [mm]      
                      wtd_uncert   = 0.1,         &! uncertaintiy in wtd calcultations [mm]
                      wtd_max      = 1000000.0,    &! maximum wtd [mm]
                      wtd_min      = 100.0,       &! minimum wtd [mm]
                      dri          = 1.0           !ratio of density of ice to density of liquid [unitless]

   INTEGER, PARAMETER :: wtd_iter_max = 20 ! maximum number of iterations to find the water table depth                    
   
   REAL :: cp    ! specific heat capacity for air
   
   !jhan:make parameter
   REAL :: max_glacier_snowd
 
   ! This module contains the following subroutines:
   PUBLIC soil_snow_gw,calc_srf_wet_fraction,calc_soil_hydraulic_props ! must be available outside this module
   !PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust ,snow_accum, stempv,trimb
   PRIVATE calc_equilibrium_water_content
   !PRIVATE GWsoilfreeze!, remove_trans,
   PRIVATE iterative_wtd,simple_wtd
   PRIVATE smoistgw, ovrlndflx

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  !-------------------------------------------------------------------------
  SUBROUTINE ovrlndflx (dels, ktau, ssnow, soil,veg, md_prin )
  USE cable_common_module

  IMPLICIT NONE
    REAL, INTENT(IN)                         :: dels ! integration time step (s)
    INTEGER, INTENT(IN)                      :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters
    TYPE(veg_parameter_type) , INTENT(IN)    :: veg  ! veg parameters
    LOGICAL, INTENT(IN)                      :: md_prin
    INTEGER, PARAMETER                       :: ntest = 0 ! for snow diag prints
    INTEGER, PARAMETER                       :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER                                  :: k, i, j
    REAL, DIMENSION(mp)                :: rnof5
    REAL, DIMENSION(mp)                :: sgamm
    REAL, DIMENSION(mp)                :: smasstot
    REAL, DIMENSION(mp,0:3)            :: smelt1                   !snow melt
    REAL(r_2), DIMENSION(mp)           :: icef,efpor               !tmp vars, fraction of ice in gridcell
    REAL(r_2)                          :: tmpa,tmpb,qinmax         !tmp vars, maximum infiltration [mm/s]
    REAL(r_2), DIMENSION(mp)           :: satfrac_liqice,S       !saturated fraction of cell, wtd in m
    REAL(r_2)                          :: liqmass,icemass,totmass  !liquid mass,ice mass, total mass [mm]
    REAL(r_2), parameter               :: pi=3.1415926535898
    REAL(r_2)                          :: fice
    REAL(r_2)                          :: dzmm,slopeSTDmm
    logical                            :: prinall = .false.  !for debugging

    TYPE ( issnow_type ) :: C 


    CALL point2constants( C )
    
   !For now assume there is no puddle
   dzmm = 1000._r_2 * soil%zse(1)

   do i = 1,mp
      efpor(i) = max(0.001_r_2, soil%ssat_vec(i,1) - ssnow%wbice(i,1))
      icemass  = ssnow%wbice(i,1) * dzmm * dri
      liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm
      totmass  = max(liqmass+icemass,real(1e-2,r_2))
      !icef(i)     = max(0._r_2,min(1._r_2,1.25_r_2*icemass / totmass))
      icef(i)     = max(0._r_2,min(1._r_2,gw_params%IceBeta*icemass / totmass))
   end do
   S(:) = 0._r_2
   do k=1,ms
     S(:) = S(:) + max(0.01,min(1.0, (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )*soil%zse(k)
   end do
   S(:) = S(:)/sum(soil%zse(1:ms),dim=1)
   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      !fice = (exp(-3._r_2*(1._r_2-icef(i)))-exp(-3._r_2))/(1._r_2-exp(-3._r_2))
      fice = (exp(gw_params%IceAlpha*(1._r_2-icef(i)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
      fice  = min(max(fice,0._r_2),1._r_2)
      !Saturated fraction
       if (gw_params%MaxSatFraction .gt. 1e-7) then 
          slopeSTDmm = sqrt(min(max(gw_params%MaxSatFraction*soil%slope_std(i),1e-5),10000._r_2)) ! ensure some variability
          ssnow%satfrac(i)    = max(1e-6,min(0.95,1._r_2 - erf( slopeSTDmm / sqrt(2.0* S(i)) ) ) )  
       else
          ssnow%satfrac(i) = 0. 
       end if
       satfrac_liqice(i)   = max(1e-6,min(0.95,fice + (1._r_2-fice)*ssnow%satfrac(i) ) )
   end do

   do i=1,mp
      tmpa = ssnow%wbliq(i,1) / efpor(i)
      tmpb = max( (tmpa-satfrac_liqice(i))/max(0.01_r_2,(1._r_2-satfrac_liqice(i))), 0._r_2)
      tmpa = -2._r_2*soil%bch_vec(i,1)*soil%sucs_vec(i,1)/dzmm
      qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*soil%hyds_vec(i,1)*exp(-gw_params%hkrz*(0.5*dzmm/1000.0_r_2-gw_params%zdepth))

      ssnow%rnof1(i) = satfrac_liqice(i) * ssnow%fwtop(i) + &
                         (1._r_2-satfrac_liqice(i))*max((ssnow%fwtop(i)-qinmax) , 0._r_2)

      ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

   end do  !mp

  !add back to the lakes to keep saturated instead of drying
  where (veg%iveg .eq. 16)
     ssnow%fwtop(:) = ssnow%fwtop(:) + ssnow%rnof1(:)
     ssnow%rnof1(:) = 0._r_2
  endwhere
           
   !---  glacier formation
   rnof5= 0.

   IF (nglacier == 2) THEN
      smelt1=0.
      WHERE( ssnow%snowd > max_glacier_snowd )

         rnof5 = MIN( 0.1, ssnow%snowd - max_glacier_snowd )

         !---- change local tg to account for energy - clearly not best method
         WHERE( ssnow%isflag == 0 )
            smasstot = 0.0
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - rnof5 * C%HLF                    &
                             / REAL( ssnow%gammzz(:,1) )
            ssnow%snowd = ssnow%snowd - rnof5
         ELSEWHERE
            smasstot = ssnow%smass(:,1) + ssnow%smass(:,2) + ssnow%smass(:,3)
         END WHERE

      END WHERE

      DO k = 1, 3
         
         WHERE( ssnow%snowd > max_glacier_snowd  .AND.  ssnow%isflag > 0 )
            sgamm = ssnow%ssdn(:,k) * cgsnow * ssnow%sdepth(:,k)
            smelt1(:,k) = MIN( rnof5 * ssnow%smass(:,k) / smasstot,            &
                          0.2 * ssnow%smass(:,k) )
            ssnow%smass(:,k) = ssnow%smass(:,k) - smelt1(:,k)
            ssnow%snowd = ssnow%snowd - smelt1(:,k)
         END WHERE
      
      END DO
   
      WHERE( ssnow%isflag > 0 ) rnof5 = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)

      ssnow%rnof1 = ssnow%rnof1 + rnof5/dels   !include this runoff in suface runoff term
   
   END IF

  END SUBROUTINE ovrlndflx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE simple_wtd(ssnow, soil, veg, ktau, md_prin)
  IMPLICIT NONE
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
  INTEGER, INTENT(IN)                       :: ktau  ! integration step number
  LOGICAL, INTENT(IN)                       :: md_prin  !print info?

  REAL(r_2), DIMENSION(mp)            :: fz, wmean,ztot
  REAL(r_2), DIMENSION(mp,ms)         :: stot
  INTEGER                             :: k,i

  do i=1,mp
     wmean(i) = 0._r_2
     fz(i)    = 5._r_2
     ztot(i)  = 0._r_2
     stot(i,:) = (ssnow%wb(i,:)-soil%watr(i,:)) / (soil%ssat_vec(i,:)-soil%watr(i,:))
  end do
  do k  = 1, ms
     do i=1,mp
        wmean(i) = wmean(i) + stot(i,k)*soil%zse(k)*1000._r_2
        ztot(i)  = ztot(i) + soil%zse(k)*1000._r_2
     end do
  end do

  do i=1,mp
     wmean(i) = wmean(i) + ssnow%GWwb(i)/soil%GWssat_vec(i) * soil%GWdz(i)*1000._r_2
     ztot(i)  = ztot(i) + soil%GWdz(i)*1000._r_2

     ssnow%wtd(i) = min(200000._r_2, fz(i) * (ztot(i) - wmean(i)))
  end do

  END SUBROUTINE simple_wtd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


  !----------------------------------------------------------------------
  ! SUBROUTINE iterative_wtd
  !
  ! Iteratively calcs the water table depth by equating the mass of water in the
  ! soil column to the mass of a hydrostatic column inegrated from the surface to the 
  ! water table depth
  !  
  SUBROUTINE iterative_wtd (ssnow, soil, veg, ktau, md_prin)
  IMPLICIT NONE
  TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
  TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
  TYPE (veg_parameter_type), INTENT(IN)     :: veg
  INTEGER, INTENT(IN)                       :: ktau  ! integration step number
  LOGICAL, INTENT(IN)                       :: md_prin  !print info?

 
  !Local vars 
  REAL(r_2), DIMENSION(mp,ms)   :: dzmm_mp,tmp_def
  REAL(r_2), DIMENSION(0:ms)    :: zimm
  REAL(r_2), DIMENSION(ms)      :: zmm
  REAL(r_2), DIMENSION(mp)      :: GWzimm,temp
  REAL(r_2), DIMENSION(mp)      :: def,defc,total_depth_column

  REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
  REAL(r_2), DIMENSION(mp)      :: invB,Nsucs_vec  !inverse of C&H B,Nsucs_vec
  INTEGER :: k,i,wttd,jlp

  !make code cleaner define these here 
  invB     = 1._r_2/soil%bch_vec(:,ms)                                !1 over C&H B
  Nsucs_vec  = soil%sucs_vec(:,ms)                                !psi_saturated mm
  dzmm_mp  = real(spread((soil%zse(:)) * 1000.0,1,mp),r_2)    !layer thickness mm
  zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm

  do k=1,ms
    zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2
  end do
  total_depth_column(:) = zimm(ms) + soil%GWdz(:)*1000._r_2
  
  !find the deficit if the water table is at the bottom of the soil column
  do i=1,mp
     defc(i) = (soil%ssat_vec(i,ms))*(total_depth_column(i)+Nsucs_vec(i)/(1._r_2-invB(i))*            &
             (1._r_2-((Nsucs_vec(i)+total_depth_column(i))/Nsucs_vec(i))**(1._r_2-invB(i)))) 
     defc(i) = max(0.1_r_2,defc(i)) 
  end do

  def(:) = 0._r_2
  do k=1,ms
     do i=1,mp

       def(i) = def(i) +                                                           &
                max(0._r_2,(soil%ssat_vec(i,k)-(ssnow%wbliq(i,k)+dri*ssnow%wbice(i,k)))*dzmm_mp(i,k))
      end do  !mp
  end do  !ms

  do i=1,mp
    def(i) = def(i) + max(0._r_2,soil%GWssat_vec(i)-ssnow%GWwb(i))*soil%GWdz(i)*1000._r_2
  end do   

 ssnow%wtd(:) = total_depth_column(:)*def(:)/defc(:)

 !use newtons method to solve for wtd, note this assumes homogenous column but
 !that is ok 
  do i=1,mp
    !mrd561 debug remove iveg = 16 test here
    if ((soil%isoilm(i) .ne. 9) .and. (veg%iveg(i) .ne. 16)) then      

      if (defc(i) > def(i)) then                 !iterate tfor wtd

        jlp=0

        mainloop: DO

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
          derv    = (soil%ssat_vec(i,ms))*(tempa-tempb) + &
                                       soil%ssat_vec(i,ms)

          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa   = 1.0_r_2
          tempb   = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
          deffunc = (soil%ssat_vec(i,ms))*(ssnow%wtd(i) +&
                           Nsucs_vec(i)/(1-invB(i))* &
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
          tempa    = (abs(tmpc/Nsucs_vec(i)))**(-invB(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(-invB(i))
          derv     = (soil%ssat_vec(i,ms))*(tempa-tempb)
          if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

          tempa    = (abs((Nsucs_vec(i)+ssnow%wtd(i)-total_depth_column(i))/Nsucs_vec(i)))**(1._r_2-invB(i))
          tempb    = (1._r_2+ssnow%wtd(i)/Nsucs_vec(i))**(1._r_2-invB(i))
          deffunc  = (soil%ssat_vec(i,ms))*(total_depth_column(i) +&
                     Nsucs_vec(i)/(1._r_2-invB(i))*(tempa-tempb))-def(i)
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
  where (ssnow%wtd(:) .gt. wtd_max) ssnow%wtd(:) = wtd_max
  where (ssnow%wtd(:) .lt. wtd_min) ssnow%wtd(:) = wtd_min


  END SUBROUTINE iterative_wtd

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoistgw (fwtop,dt,ktau,ssnow,soil,prin)
  ! solves the modified richards equation (Zeng and Decker 2009) to find
  ! vertical mocement of soil water.  Bottom boundary condition is determined
  ! using a single layer groundwater module
  !
  SUBROUTINE smoistgw (dels,ktau,ssnow,soil,veg,canopy,md_prin)
  USE cable_common_module

  IMPLICIT NONE
  
    REAL, INTENT(IN)                          :: dels  ! time step size (s)
    INTEGER, INTENT(IN)                       :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
    TYPE(canopy_type), INTENT(INOUT)          :: canopy ! vegetation variables
    LOGICAL, INTENT(IN)                       :: md_prin

    TYPE ( issnow_type ) :: C 
    
    !Local variables.  
    REAL(r_2), DIMENSION(mp,ms+1)       :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms+1)       :: rt

    INTEGER                             :: k,kk,i
    REAL(r_2), DIMENSION(mp,ms)         :: eff_por,old_wb,mss_por  !effective porosity (mm3/mm3),wb(mm3/mm3),mass (mm) of eff_por
    REAL(r_2), DIMENSION(mp,ms)         :: msliq,msice             !mass of the soil liquid and ice water    
    REAL(r_2), DIMENSION(mp,ms)         :: wbrat                   !ratio of volumetric water - watr / saturated - watr
    REAL(r_2), DIMENSION(mp)            :: den
    REAL(r_2), DIMENSION(mp)            :: dne
    REAL(r_2), DIMENSION(mp)            :: num
    REAL(r_2), DIMENSION(mp)            :: qin
    REAL(r_2), DIMENSION(mp)            :: qout
    REAL(r_2), DIMENSION(mp)            :: dqidw0
    REAL(r_2), DIMENSION(mp)            :: dqidw1
    REAL(r_2), DIMENSION(mp)            :: dqodw0
    REAL(r_2), DIMENSION(mp)            :: dqodw1,dqodw2
    REAL(r_2), DIMENSION(mp)            :: s1,s2,tmpi,temp0,voleq1,tempi
    REAL(r_2), DIMENSION(ms)            :: dzmm
    REAL(r_2), DIMENSION(mp,ms)         :: dzmm_mp
    REAL(r_2), DIMENSION(0:ms+1)        :: zimm
    REAL(r_2), DIMENSION(ms)            :: zmm
    REAL(r_2), DIMENSION(mp)            :: GWzimm,xs,zaq,s_mid,GWdzmm
    REAL(r_2), DIMENSION(mp)            :: xs1,GWmsliq!xsi    !mass (mm) of liquid over/under saturation, mass of aquifer water
    REAL(r_2)                           :: xsi
    REAL(r_2), DIMENSION(mp,ms+1)       :: qhlev,del_wb
    REAL(r_2), DIMENSION(mp)            :: sm_tot,drainmod  !total column soil water available for drainage
    INTEGER, DIMENSION(mp)              :: idlev
    logical                             :: prinall = .false.   !another debug flag
    character (len=30)                  :: fmt  !format to output some debug info
    !MD DEBUG VARS
    INTEGER :: imp,ims,k_drain

    drainmod(:) = 1._r_2  !parameter to modify qhrz params by basin or veg type
    !where(veg%iveg .eq. 2) drainmod(:) = 0.1_r_2*drainmod(:)
    !drainmod(:) = (1._r_2 + soil%FOrg(:,1))*drainmod(:)


    fmt='(A6,6(1X,F8.6))'       !not needed.  was used to nicely format debug output
    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    dzmm_mp = spread(dzmm,1,mp)
    zimm(0) = 0._r_2

    zimm(1:ms) = zimm(0:(ms-1)) + dzmm(1:ms)
    zmm(1:ms)  = zimm(0:(ms-1)) + 0.5_r_2*dzmm(1:ms)

    GWdzmm(:) = real(soil%GWdz(:),r_2)*1000._r_2
    GWzimm(:) = zimm(ms)+GWdzmm(:)
    zaq(:)    = zimm(ms) + 0.5_r_2*GWdzmm(:)
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
    old_wb(:,:) = ssnow%wb(:,:)
    
    !equilibrium water content
    CALL calc_equilibrium_water_content(ssnow,soil)

    CALL calc_soil_hydraulic_props(ssnow,soil,veg)

    do i=1,mp

       !Note: future revision will have interaction with river here. nned to
       !work on router and add river type cells
       ssnow%qhz(i)  = min(max(tan(soil%slope(i)),0.001),0.1) * drainmod(i)*gw_params%MaxHorzDrainRate* &
                    exp(-ssnow%wtd(i)/(1000._r_2*(gw_params%EfoldHorzDrainRate)))

 
       !identify first no frozen layer.  drinage from that layer and below
       !drain from sat layers
       k_drain = ms+1
       do k=ms,3,-1
           !below what was in paper
          !if ( ssnow%wbliq(i,k+1) .le. gw_params%frozen_frac*ssnow%wb(i,k+1)  .and.&
          !     ssnow%wbliq(i,k  ) .gt. gw_params%frozen_frac*ssnow%wb(i,k  ) ) then
          if (ssnow%wtd(i) .le. sum(dzmm(1:k),dim=1)) then
             k_drain = k + 1
          end if
       end do
       k_drain = max(k_drain,4)

!       k_drain = ms + 1
!       do k=ms,2,-1
!          if ((all(ssnow%wbliq(i,k:min(k_drain,ms)) .ge. 0.95*soil%ssat_vec(i,k:min(k_drain,ms)))) .and. &
!                 (ssnow%wbliq(i,k-1) .lt. 0.95*soil%ssat_vec(i,k-1)))) then
!             k_drain = k
!          end if
!       end do

       qhlev(i,:) = 0._r_2
       sm_tot(i) = 0._r_2
       if (k_drain .le. ms) then
          do k=k_drain,ms
             sm_tot(i) = sm_tot(i) + max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)!*dzmm(k)
          end do

          sm_tot(i) = sm_tot(i) + &
                      max(ssnow%GWwb(i)-soil%watr(i,ms),0._r_2)*max(1._r_2-ssnow%fracice(i,ms),0._r_2)

          sm_tot(i) = max(sm_tot(i),0.01_r_2)

         do k=k_drain,ms
             qhlev(i,k) = ssnow%qhz(i)*max(ssnow%wbliq(i,k)-soil%watr(i,k),0._r_2)/sm_tot(i)!*dzmm(k)/sm_tot(i)
         end do

         qhlev(i,ms+1) = ssnow%qhz(i)*max(ssnow%GWwb(i)-soil%watr(i,ms),0._r_2)*max(1._r_2-ssnow%fracice(i,ms),0._r_2)/sm_tot(i) 

       else
          qhlev(i,ms+1) = max(1._r_2-ssnow%fracice(i,ms),0._r_2)*ssnow%qhz(i)!*max(ssnow%GWwb(i)-soil%watr(i,ms),0._r_2)/sm_tot(i)
       end if

       !incase every layer is frozen very dry
       ssnow%qhz(i) = qhlev(i,ms+1)
       do k=k_drain,ms
          ssnow%qhz(i) = ssnow%qhz(i) +qhlev(i,k)
       end do

       !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
       !overflowing
       if (soil%isoilm(i) .eq. 9 .or. veg%iveg(i) .ge. 16) then
          ssnow%qhz(i) = 0._r_2
          qhlev(i,:) = 0._r_2
       end if

    end do  


    rt(:,:) = 0._r_2; at(:,:) = 0._r_2     !ensure input to tridiag is valid
    bt(:,:) = 1._r_2; ct(:,:) = 0._r_2

    k = 1     !top soil layer
    do i=1,mp
       qin(i)     = ssnow%sinfil(i)
       den(i)     = (zmm(k+1)-zmm(k))
       dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
       num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
       dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
       dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
       rt(i,k) =  qin(i) - qout(i)! - qhlev(i,k) - ssnow%rex(i,k)
       at(i,k) =  0._r_2
       bt(i,k) =  dzmm(k)/dels + dqodw1(i)
       ct(i,k) =  dqodw2(i)      
    end do
    do k = 2, ms - 1     !middle soil layers
       do i=1,mp
          den(i)     = (zmm(k) - zmm(k-1))
          dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
          num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
          qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
          dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
          den(i)     = (zmm(k+1)-zmm(k))
          dne(i)     = (ssnow%zq(i,k+1)-ssnow%zq(i,k))
          num(i)     = (ssnow%smp(i,k+1)-ssnow%smp(i,k)) - dne(i)
          qout(i)    = -ssnow%hk(i,k)*num(i)/den(i)
          dqodw1(i)  = -(-ssnow%hk(i,k)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k))/den(i)
          dqodw2(i)  = -( ssnow%hk(i,k)*ssnow%dsmpdw(i,k+1) + num(i)*ssnow%dhkdw(i,k))/den(i)
          rt(i,k) =  qin(i) - qout(i)! - qhlev(i,k) - ssnow%rex(i,k)
          at(i,k) = -dqidw0(i)
          bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
          ct(i,k) =  dqodw2(i)
       end do
    end do
       
    k = ms   !Bottom soil layer
    do i=1,mp
       den(i)     = (zmm(k) - zmm(k-1))
       dne(i)     = (ssnow%zq(i,k)-ssnow%zq(i,k-1))
       num(i)     = (ssnow%smp(i,k)-ssnow%smp(i,k-1)) - dne(i)
       qin(i)     = -ssnow%hk(i,k-1)*num(i)/den(i)
       dqidw0(i)  = -(-ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k-1) + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       dqidw1(i)  = -( ssnow%hk(i,k-1)*ssnow%dsmpdw(i,k)   + num(i)*ssnow%dhkdw(i,k-1))/den(i)
       den(i)     = zaq(i) - zmm(k)
       dne(i)     = (ssnow%GWzq(i)-ssnow%zq(i,k))
       num(i)     =  (ssnow%GWsmp(i)-ssnow%smp(i,k)) - dne(i)
       qout(i)    = 0._r_2
       dqodw1(i)  = 0._r_2
       dqodw2(i)  = 0._r_2
       rt(i,k) =  qin(i) - qout(i)
       at(i,k) = -dqidw0(i)
       bt(i,k) =  dzmm(k)/dels - dqidw1(i) + dqodw1(i)
       ct(i,k) =  dqodw2(i) 
    end do
       
    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
!    do i=1,mp
!       if ((ssnow%wtd(i) .le. sum(dzmm,dim=1)) .or. (veg%iveg(i) .ge. 16) .or. (soil%isoilm(i) .eq. 9))  then
!          ssnow%Qrecharge(i) = 0._r_2
!       else
!          ssnow%Qrecharge(i) = -ssnow%hk(i,ms)*((ssnow%GWsmp(i)-ssnow%smp(i,ms)) - (ssnow%GWzq(i)-ssnow%zq(i,ms)))/(zaq(i) - zmm(ms))
!       end if
!    end do

    CALL aquifer_recharge(ssnow,soil,veg,zaq,zmm,dzmm)

    CALL trimb(at,bt,ct,rt,ms)                       !use the defulat cable tridiag solution

    do k=1,ms
       do i=1,mp
          ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + rt(i,k) - qhlev(i,k)*dels/dzmm(k)   !volutermic liquid
       end do
    end do

    do i=1,mp
       ssnow%wbliq(i,ms) = ssnow%wbliq(i,ms) - ssnow%Qrecharge(i)*dels/dzmm(ms)
    end do
    do i=1,mp
       ssnow%GWwb(i) = ssnow%GWwb(i)  +  (ssnow%Qrecharge(i)-qhlev(i,ms+1))*dels/GWdzmm(i)
    end do

    !determine the available pore space
    !volumetric
    eff_por(:,:)  = soil%ssat_vec(:,:) - ssnow%wbice(:,:)
    do i=1,mp
       xsi = 0._r_2

       if (ssnow%GWwb(i) .gt. soil%GWssat_vec(i)) then
          xsi = (ssnow%GWwb(i) - soil%GWssat_vec(i))*GWdzmm(i)
          ssnow%GWwb(i) = soil%GWssat_vec(i)
       end if

       do k=1,ms
          if (ssnow%wbliq(i,k) .gt. eff_por(i,k)) then
             xsi = xsi + (ssnow%wbliq(i,k) - eff_por(i,k))*dzmm(k)
             ssnow%wbliq(i,k) = eff_por(i,k)
           end if
       end do
     
       do k = ms,1,-1  !loop from bottom to top adding extra water to each layer
          if (xsi .gt. 0._r_2) then
             if (xsi .lt. (eff_por(i,k)-ssnow%wbliq(i,k))*dzmm(k)) then
                ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xsi/dzmm(k)
                xsi = 0._r_2
             else
                xsi = xsi - (eff_por(i,k) - ssnow%wbliq(i,k))*dzmm(k)
                ssnow%wbliq(i,k) = eff_por(i,k)
             end if
          end if
       end do  !ms loop

       if (xsi .gt. 0._r_2) then
          ssnow%qhz(i) = ssnow%qhz(i) + xsi/dels
          xsi = 0._r_2
       end if

       do k = 1,ms
          xsi = 0._r_2             !should be a single float (array not needed)
          if (ssnow%wbliq(i,k) .lt. volwatmin) then
             xsi = (volwatmin - ssnow%wbliq(i,k))*dzmm(k)  !in mm
             ssnow%wbliq(i,k) = volwatmin
             if (k .lt. ms) then
                ssnow%wbliq(i,k+1) = ssnow%wbliq(i,k+1) - xsi/dzmm(k+1)
             else
                ssnow%GWwb(i) = ssnow%GWwb(i) - xsi / GWdzmm(i)
             end if
          end if
       end do  !ms loop
 
       if ( (ssnow%GWwb(i) .lt. volwatmin) .and. (soil%isoilm(i) .ne. 9) ) then
          xsi = (volwatmin - ssnow%GWwb(i)) / GWdzmm(i)  !mm
          ssnow%GWwb(i) = volwatmin
          ssnow%qhz(i) = ssnow%qhz(i) - xsi / dels
       end if
         
       !update mass variables
       ssnow%wmliq(i,:)      = ssnow%wbliq(i,:) * dzmm_mp(i,:)
       ssnow%wmice(i,:)      = ssnow%wbice(i,:) * dzmm_mp(i,:) * dri
       ssnow%wb(i,:)         = ssnow%wbliq(i,:) + ssnow%wbice(i,:)
       ssnow%wmtot(i,:)      = ssnow%wmliq(i,:) + ssnow%wmice(i,:)
       ssnow%rnof2(i)        = ssnow%qhz(i)               !rnof2 is default cable deep runoff var 

   end do  !mp loop
       

 END SUBROUTINE smoistgw



! Inputs:
!	 dt_in - time step in sec
!	 ktau_in - time step no.
!	 ga	 - ground heat flux W/m^2
!	 dgdtg	 -
!	 condxpr - total precip reaching the ground (liquid and solid)
!	 scondxpr - precip (solid only)
!	 fev   - transpiration (W/m2)
!	 fes   - soil evaporation (W/m2)
!	 isoil - soil type
!	 ivegt - vegetation type
! Output
!	 ssnow
SUBROUTINE soil_snow_gw(dels, soil, ssnow, canopy, met, bal, veg)
   USE cable_IO_vars_module, ONLY: wlogn

   USE cable_common_module
   REAL                     , INTENT(IN)     :: dels ! integration time step (s)
   TYPE(soil_parameter_type), INTENT(INOUT)  :: soil
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(canopy_type)        , INTENT(INOUT)  :: canopy
   TYPE(veg_parameter_type) , INTENT(INOUT)  :: veg
   TYPE(met_type)           , INTENT(INOUT)  :: met ! all met forcing
   TYPE (balances_type)     , INTENT(INOUT)  :: bal

   INTEGER             :: k,i
   REAL, DIMENSION(mp) :: snowmlt
   REAL, DIMENSION(mp) :: totwet
   REAL, DIMENSION(mp) :: weting,GWwb_ic,wberr
   REAL, DIMENSION(mp) :: tgg_old, tggsn_old,wbtot_ic,del_wbtot
   REAL(r_2), DIMENSION(mp) :: xx,deltat,sinfil1,sinfil2,sinfil3 
   REAL                :: zsetot
   INTEGER, SAVE :: ktau =0 
   LOGICAL :: prin,md_prin
   REAL(r_2) :: wb_lake_T, rnof2_T, ratio

   TYPE ( issnow_type ) :: C
   
   prin = .FALSE.
   md_prin = .false.
   
   CALL point2constants( C ) 
   cp = C%CAPP
    
   ktau = ktau +1 

   !jhan - make switchable 
   ! appropriate for ACCESS1.0
   !max_glacier_snowd = 50000.0
   ! appropriate for ACCESS1.3 
   max_glacier_snowd = 1100.0
   
   zsetot = sum(soil%zse) 
   ssnow%tggav = 0.

   DO k = 1, ms

      ssnow%tggav = ssnow%tggav  + soil%zse(k)*ssnow%tgg(:,k)/zsetot

   END DO


   IF( cable_runtime%offline .or. cable_runtime%mk3l ) ssnow%t_snwlr = 0.05_r_2

   ssnow%fwtop1 = 0.0
   ssnow%fwtop2 = 0.0
   ssnow%fwtop3 = 0.0
   ssnow%runoff = 0.0 ! initialise total runoff
   ssnow%rnof1  = 0.0 ! initialise surface runoff
   ssnow%rnof2  = 0.0 ! initialise deep drainage
   ssnow%smelt  = 0.0 ! initialise snowmelt
   ssnow%dtmlt  = 0.0 
   ssnow%osnowd = ssnow%snowd
   ! Scaling  runoff to kg/m^2/s (mm/s) to match rest of the model
   ssnow%sinfil = 0.0   

   IF( .NOT.cable_user%cable_runtime_coupled ) THEN
   
      IF( ktau_gl <= 1 ) THEN
         
         IF (cable_runtime%um) canopy%dgdtg = 0.0 ! RML added um condition
                                                  ! after discussion with BP
         ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
         ssnow%wbtot = 0.0
         ssnow%wb(:,:)  = MIN( soil%ssat_vec(:,:), MAX ( ssnow%wb(:,:), soil%swilt_vec(:,:) ) )   

         DO k = 1, ms
            
            WHERE( ssnow%tgg(:,k) <= C%TFRZ .AND. ssnow%wbice(:,k) <= 0.01 )   &
               ssnow%wbice(:,k) = 0.5 * ssnow%wb(:,k)

            WHERE( ssnow%tgg(:,k) < C%TFRZ)                                    &
               ssnow%wbice(:,k) = 0.8 * ssnow%wb(:,k)
            
         END DO

         WHERE ( soil%isoilm .eq. 9 ) 

            ! permanent ice: fix hard-wired number in next version
            ssnow%snowd = max_glacier_snowd
            ssnow%osnowd = max_glacier_snowd
            ssnow%tgg(:,1) = ssnow%tgg(:,1) - 1.0

         END WHERE

         WHERE ( spread(soil%isoilm,2,ms) .eq. 9 )

              ssnow%wb    = 0.95 * soil%ssat_vec
              ssnow%wbice = 0.90 * ssnow%wb

         END WHERE
         
         xx=soil%css * soil%rhosoil

         ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat_vec(:,1)) * soil%css * soil%rhosoil &
              & + (ssnow%wb(:,1) - ssnow%wbice(:,1) ) * cswat * C%denliq &
              & + ssnow%wbice(:,1) * csice * C%denice, xx ) * soil%zse(1)

      END IF

   ENDIF  ! if(.NOT.cable_runtime_coupled)

   !Start with wb and wbice.  Need wbliq, wmliq,wmice,wmtot
   !find the mass of ice and liq from the prognostic volumetric values
   ssnow%wbliq = ssnow%wb - ssnow%wbice                     !liquid volume
   ssnow%wmice = ssnow%wbice*real(C%denice*spread(soil%zse,1,mp),r_2) !ice mass
   ssnow%wmliq = ssnow%wbliq*real(C%denliq*spread(soil%zse,1,mp),r_2) !liquid mass
   ssnow%wmtot = ssnow%wmice + ssnow%wmliq                  !liq+ice mass

   xx=soil%css * soil%rhosoil
   IF (ktau <= 1)                                                              &
     ssnow%gammzz(:,1) = MAX( (1.0 - soil%ssat_vec(:,1)) * soil%css * soil%rhosoil      &
            & + ssnow%wbliq(:,1) * cswat * C%denliq           &
            & + ssnow%wbice(:,1) * csice * C%denice, xx ) * soil%zse(1) +   &
            & (1. - ssnow%isflag) * cgsnow * ssnow%snowd

   ssnow%wblf   = max(ssnow%wbliq/soil%ssat_vec,0.01_r_2)
   ssnow%wbfice = max(ssnow%wbice/soil%ssat_vec,0._r_2)

   !initial water in the soil column
   wbtot_ic  = sum(ssnow%wbliq(:,:)*C%denliq*spread(soil%zse,1,mp),2) + &
               sum(ssnow%wbice(:,:)*C%denice*spread(soil%zse,1,mp),2) + &
               ssnow%GWwb(:)*soil%GWdz*C%denliq
               
   GWwb_ic = ssnow%GWwb

   CALL snowcheck (dels, ssnow, soil, met )

   CALL snowdensity (dels, ssnow, soil)

   CALL snow_accum (dels, canopy, met, ssnow, soil )

   CALL snow_melting (dels, snowmlt, ssnow, soil )

   ! Add snow melt to global snow melt variable:
   ssnow%smelt = snowmlt

   ! Adjust levels in the snowpack due to snow accumulation/melting,
   ! snow aging etc...

   CALL snowl_adjust(dels, ssnow, canopy )

        
   CALL stempv(dels, canopy, ssnow, soil)

   !do the soil and snow melting, freezing prior to water movement
   ssnow%tss =  (1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   CALL snow_melting (dels, snowmlt, ssnow, soil )
   
   ! Add new snow melt to global snow melt variable: 
   ssnow%smelt = ssnow%smelt + snowmlt

   CALL remove_trans(dels, soil, ssnow, canopy, veg)        !transpiration loss per soil layer

   CALL  soilfreeze(dels, soil, ssnow)

   ssnow%fwtop = canopy%precis/dels + ssnow%smelt/dels   !water from canopy and snowmelt [mm/s]   
   !ssnow%rnof1 = ssnow%rnof1 + ssnow%smelt / dels          !adding snow melt directly to the runoff

   CALL iterative_wtd (ssnow, soil, veg, ktau, md_prin)  
   !CALL simple_wtd(ssnow, soil, veg, ktau, md_prin)

   CALL ovrlndflx (dels, ktau, ssnow, soil, veg, md_prin )         !surface runoff, incorporate ssnow%pudsto?

   ssnow%sinfil = ssnow%fwtop - canopy%segg  !canopy%fes/C%HL               !remove soil evap from throughfall
   !ssnow%pudsto = max(ssnow%pudsto - canopy%fesp/C%HL*dels,0._r_2)  !currently pudsto = 0.0 always

   CALL smoistgw (dels,ktau,ssnow,soil,veg,canopy, md_prin)               !vertical soil moisture movement. 
  
   ! correction required for energy balance in online simulations 
   IF( cable_runtime%um) THEN

      canopy%fhs_cor = ssnow%dtmlt(:,1)*ssnow%dfh_dtg
      canopy%fes_cor = ssnow%dtmlt(:,1)*(ssnow%cls*ssnow%dfe_ddq * ssnow%ddq_dtg)

      canopy%fhs = canopy%fhs+canopy%fhs_cor
      canopy%fes = canopy%fes+canopy%fes_cor

   ENDIF

   ssnow%smelt  = ssnow%smelt/dels    !change units to mm/s.  cable_driver then reverts back to mm
   ssnow%runoff = (ssnow%rnof1 + ssnow%rnof2)!*dels  !cable_driver converts from mm/s to mm
                                                     !rnof1 and rnof2 are already in mm/s

   ! Set weighted soil/snow surface temperature
   ssnow%tss=(1-ssnow%isflag)*ssnow%tgg(:,1) + ssnow%isflag*ssnow%tggsn(:,1)

   !total water mass at the end of the soilsnow_GW routine
   ssnow%wbtot = sum(ssnow%wbliq(:,:)*real(C%denliq,r_2)*real(spread(soil%zse,1,mp),r_2),2) + &
                 sum(ssnow%wbice(:,:)*real(C%denice,r_2)*real(spread(soil%zse,1,mp),r_2),2) + &
                 ssnow%GWwb(:)*soil%GWdz*real(C%denliq,r_2)
                 
   !for debug water balance.  del_wbtot = fluxes = infiltration [though-evap] - trans - qhorz drainage
   del_wbtot   = dels * (ssnow%sinfil - ssnow%rnof2 - canopy%fevc / C%HL)
   !set below to keep track of water imbalance within the GW module explicitly.  also must change cable_checks
   ssnow%wbtot = ssnow%wbtot-wbtot_ic

END SUBROUTINE soil_snow_gw


SUBROUTINE calc_equilibrium_water_content(ssnow,soil)

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    !local variables
    REAL(r_2), dimension(mp)    :: zaq      !node depth of the aquifer
    REAL(r_2), dimension(ms)    :: dzmm     !layer thickness for single tile
    REAL(r_2), dimension(mp)    :: GWdzmm   !aquifer thickness at each tile
    REAL(r_2), dimension(mp)    :: GWzimm   !aquifer layer interface depth
    REAL(r_2), dimension(0:ms)  :: zimm     !layer interface depth in mm  
    REAL(r_2), dimension(ms)    :: zmm      !node depths in mm
    REAL(r_2)                   :: tempi, temp0,voleq1,wbrat

    INTEGER :: k,i

    !make code cleaner define these here
    dzmm    = 1000.0_r_2 * real(soil%zse(:),r_2)
    zimm(0) = 0._r_2

    do k=1,ms
       zimm(k) = zimm(k-1) + dzmm(k)
       zmm(k)  = zimm(k-1) + 0.5_r_2*dzmm(k)
    end do 

    do i=1,mp
       GWdzmm(i) = real(soil%GWdz(i),r_2)*1000._r_2
       GWzimm(i) = zimm(ms)+GWdzmm(i)
       zaq(i)    = zimm(ms) + 0.5_r_2*GWdzmm(i)
    end do
    !equilibrium water content
    do k=1,ms
       do i=1,mp

          if ((ssnow%wtd(i) .le. zimm(k-1))) then         !fully saturated

             ssnow%wbeq(i,k) = soil%ssat_vec(i,k)
             
          elseif ((ssnow%wtd(i) .le. zimm(k)) .and. (ssnow%wtd(i) .gt. zimm(k-1))) then

             tempi = 1._r_2
             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k-1))/soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))               
             voleq1 = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
                       (1._r_2-1._r_2/soil%bch_vec(i,k))/(ssnow%wtd(i)-zimm(k-1))*(tempi-temp0)
             ssnow%wbeq(i,k) = (voleq1*(ssnow%wtd(i)-zimm(k-1)) + (soil%ssat_vec(i,k)-soil%watr(i,k))&
                            *(zimm(k)-ssnow%wtd(i)))/(zimm(k)-zimm(k-1)) + soil%watr(i,k)

          else

             tempi = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k))/soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))
             temp0 = (((soil%sucs_vec(i,k)+ssnow%wtd(i)-zimm(k-1))/soil%sucs_vec(i,k)))**(1._r_2-1._r_2/soil%bch_vec(i,k))   
             ssnow%wbeq(i,k) = -soil%sucs_vec(i,k)*(soil%ssat_vec(i,k)-soil%watr(i,k))/&
                               (1._r_2-1._r_2/soil%bch_vec(i,k))/(zimm(k)-zimm(k-1))*(tempi-temp0)+soil%watr(i,k)

          end if
          
          ssnow%wbeq(i,k) = min(max(ssnow%wbeq(i,k),soil%watr(i,k)),soil%ssat_vec(i,k))
          
          wbrat = min(max((ssnow%wbeq(i,k) - soil%watr(i,k))/(soil%ssat_vec(i,k)-soil%watr(i,k)),&
                          0.001_r_2),1._r_2)

          ssnow%zq(i,k) = max(-soil%sucs_vec(i,k)*(wbrat**(-soil%bch_vec(i,k))),sucmin)
          
          
       end do  !mp
    end do  !ms
 
    do i=1,mp
    !Aquifer Equilibrium water content
       if (ssnow%wtd(i) .le. zimm(ms)) then      !fully saturated

          ssnow%GWwbeq(i) = soil%GWssat_vec(i)-soil%GWwatr(i)

       elseif ((ssnow%wtd(i) .gt. GWzimm(i)))   then     !fully unsaturated

          tempi = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-GWzimm(i))/soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))
          temp0 = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))   
          ssnow%GWwbeq(i) = -soil%GWsucs_vec(i)*soil%GWssat_vec(i)/&
                          (1._r_2-1._r_2/soil%GWbch_vec(i))/(GWzimm(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)   

       else           

          tempi  = 1._r_2
          temp0  = (((soil%GWsucs_vec(i)+ssnow%wtd(i)-zimm(ms))/soil%GWsucs_vec(i)))**(1._r_2-1._r_2/soil%GWbch_vec(i))               
          voleq1 = -soil%GWsucs_vec(i)*(soil%GWssat_vec(i)-soil%GWwatr(i))/&
                (1._r_2-1._r_2/soil%GWbch_vec(i))/(ssnow%wtd(i)-zimm(ms))*(tempi-temp0) + soil%GWwatr(i)
          ssnow%GWwbeq(i) = (voleq1*(ssnow%wtd(i)-zimm(ms)) + (soil%GWssat_vec(i)-soil%GWwatr(i))*&
                         (GWzimm(i)-ssnow%wtd(i)))/(GWzimm(i)-zimm(ms)) + soil%GWwatr(i)

       end if

       ssnow%GWwbeq(i) = min(max(ssnow%GWwbeq(i),soil%GWwatr(i)),soil%GWssat_vec(i))

       ssnow%GWzq(i) = -soil%GWsucs_vec(i)*(max((ssnow%GWwbeq(i)-soil%GWwatr(i))/           &
                    (soil%GWssat_vec(i)-soil%GWwatr(i)),0.001_r_2))**(-soil%GWbch_vec(i))
       ssnow%GWzq(i) = max(sucmin, ssnow%GWzq(i))
       
    end do

END SUBROUTINE calc_equilibrium_water_content

SUBROUTINE calc_srf_wet_fraction(ssnow,soil)

  IMPLICIT NONE
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssnow  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(IN)    :: soil ! soil parameters

    !local variables
    REAL(r_2), DIMENSION(mp)           :: icef,satfrac_liqice,S
    REAL(r_2)                          :: fice,xx
    REAL(r_2)                          :: dzmm_one,liqmass,icemass,totmass
    INTEGER                            :: i,j,k
    REAL(r_2), parameter               :: pi=3.1415926535898
    REAL(r_2)                          :: wb_unsat,wb_lin,funcval,derv,slopeSTDmm,func_step
    REAL(r_2)                          :: wb_evap_threshold

    dzmm_one  = 1000._r_2 * real(soil%zse(1),r_2)

    do i = 1,mp
       icemass  = ssnow%wbice(i,1) * dzmm_one * dri
       liqmass  = (ssnow%wb(i,1)-ssnow%wbice(i,1)) * dzmm_one
       totmass  = max(liqmass+icemass,real(1e-2,r_2))
       icef(i)     = max(0._r_2,min(1._r_2, gw_params%IceBeta*icemass / totmass))
   end do

   !S(:) = 0._r_2
   !do k=1,2
   !   S(:) = S(:) + max(0.01,min(1.0, (ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%ssat_vec(:,k)-soil%watr(:,k)) ) )
   !end do
   !S(:) = S(:)/2._r_2

   S(:) = 0._r_2
   do k=1,ms
     S(:) = S(:) + max(0.01,min(1.0,(ssnow%wb(:,k)-ssnow%wbice(:,k)-soil%watr(:,k))/max(0.001,soil%ssat_vec(:,k)-soil%watr(:,k))) )*soil%zse(k)
   end do
   S(:) = S(:)/sum(soil%zse(1:ms),dim=1)

   !srf frozen fraction.  should be based on topography
   do i = 1,mp
      fice = (exp(gw_params%IceAlpha*(1._r_2-icef(i)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
      fice = min(1._r_2,max(0._r_2,fice))

      !Saturated fraction
       if (gw_params%MaxSatFraction .gt. 1e-7) then 
         slopeSTDmm = sqrt(min(max(gw_params%MaxSatFraction*soil%slope_std(i),1e-5),10000._r_2)) ! ensure some variability
         ssnow%satfrac(i)    = 1._r_2 - erf( slopeSTDmm / sqrt(2.0* S(i)) )
      else
         ssnow%satfrac(i)  = 0._r_2
      end if
      ssnow%satfrac(i)    = max(1e-6,min(0.95,ssnow%satfrac(i)))
      satfrac_liqice(i) = fice + (1._r_2-fice)*ssnow%satfrac(i)

      wb_unsat = ((ssnow%wb(i,1)-ssnow%wbice(i,1)) - ssnow%satfrac(i)*soil%ssat_vec(i,1))/(1.-ssnow%satfrac(i))
      wb_unsat = min(soil%ssat_vec(i,1),max(0.,wb_unsat))

      wb_evap_threshold = min( max( gw_params%SoilEvapAlpha*soil%sfc_vec(i,1), soil%swilt_vec(i,1) ), soil%ssat_vec(i,1) )

      !Sakguchi and Zeng 2009
      if (wb_unsat .ge. wb_evap_threshold) then
         xx = 1.
      else
         xx = 0.25 * (1._r_2 - cos(pi*wb_unsat/(wb_evap_threshold)))**2.0
      end if

      if (.not.cable_user%or_evap) then
         ssnow%wetfac(i) = max(0.0,min(1.0,satfrac_liqice(i) +&
                                     (1. - satfrac_liqice(i))*xx ) )
      else
         ssnow%wetfac(i) = ssnow%satfrac(i)
      end if

   end do


END SUBROUTINE calc_srf_wet_fraction

SUBROUTINE calc_soil_hydraulic_props(ssnow,soil,veg)
   USE cable_common_module
   TYPE(soil_parameter_type), INTENT(IN)     :: soil 
   TYPE(soil_snow_type)     , INTENT(INOUT)  :: ssnow
   TYPE(veg_parameter_type) , INTENT(IN)     :: veg

   INTEGER :: i,k,kk

   REAL(r_2), DIMENSION(mp) :: s1, &  !temporary variables for calculating hydraulic properties
                               s2, &
                               s_mid

   REAL(r_2), DIMENSION(0:ms) :: zimm  !depths at interface between layers

    !soil matric potential, hydraulic conductivity, and derivatives of each with respect to water (calculated using total (not liquid))

   zimm(0)  = 0.0_r_2                                          !depth of layer interfaces mm
   do k=1,ms
     zimm(k) = zimm(k-1) + soil%zse(k)*1000._r_2
   end do
   zimm(ms) = zimm(ms) + soil%GWdz(1)*1000._r_2


    do k=1,ms
       do i=1,mp
          ssnow%icefrac(i,k) = ssnow%wbice(i,k)/(max(ssnow%wb(i,k),0.01_r_2))
          ssnow%fracice(i,k) = (exp(gw_params%IceAlpha*(1._r_2-ssnow%icefrac(i,k)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
       end do
    end do

    ssnow%fracice(:,:) = max( min( ssnow%fracice, 1._r_2), 0._r_2)

    do k=1,ms-1
       kk=k+1
       do i=1,mp

          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + max(ssnow%wb(i,kk)-soil%watr(i,kk),0.)) / &
            (0.5_r_2*((soil%ssat_vec(i,k)-soil%watr(i,k)) + (soil%ssat_vec(i,kk)-soil%watr(i,kk))))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)

          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
          ssnow%hk(i,k)    =  (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))*s1(i)*s2(i)*&
                               exp(-gw_params%hkrz*(zimm(k)/1000.0_r_2-gw_params%zdepth)) 
          ssnow%dhkdw(i,k) = (1.0 - 0.5_r_2*(ssnow%fracice(i,k)+ssnow%fracice(i,kk)))* &
                             (2._r_2*soil%bch_vec(i,k)+3._r_2)*s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000.0_r_2-gw_params%zdepth))
       end do
    end do

    k = ms 
       do i=1,mp
          s1(i) = 0.5_r_2*(max(ssnow%wb(i,k)-soil%watr(i,k),0.) + max(ssnow%GWwb(i)-soil%GWwatr(i),0.)) / &
                  (0.5_r_2*(soil%ssat_vec(i,k)-soil%watr(i,k) + soil%GWssat_vec(i)-soil%GWwatr(i)))

          s1(i) = min(max(s1(i),0.01_r_2),1._r_2)

          s2(i) = soil%hyds_vec(i,k)*s1(i)**(2._r_2*soil%bch_vec(i,k)+2._r_2)
          ssnow%hk(i,k)    = s1(i)*s2(i)*(1.-ssnow%fracice(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000._r_2-gw_params%zdepth))
          ssnow%dhkdw(i,k) = (1.-ssnow%fracice(i,k))* (2._r_2*soil%bch_vec(i,k)+3._r_2)*&
                             s2(i)*0.5_r_2/(soil%ssat_vec(i,k)-soil%watr(i,k))*&
                             exp(-gw_params%hkrz*(zimm(k)/1000._r_2-gw_params%zdepth))
       end do
 
    do k=1,ms 
       do i=1,mp
          s_mid(i) = (ssnow%wb(i,k)-soil%watr(i,k))/&  !+dri*ssnow%wbice(:,k)
              (soil%ssat_vec(i,k)-soil%watr(i,k))

          s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)

          ssnow%smp(i,k) = -soil%sucs_vec(i,k)*s_mid(i)**(-soil%bch_vec(i,k))

          ssnow%smp(i,k) = max(min(ssnow%smp(i,k),-soil%sucs_vec(i,k)),sucmin)

          ssnow%dsmpdw(i,k) = -soil%bch_vec(i,k)*ssnow%smp(i,k)/&
                    (max(s_mid(i)*(soil%ssat_vec(i,k)-soil%watr(i,k)),0.001_r_2))       
       end do   
    end do

    do i=1,mp
       !Aquifer properties
       s_mid(i) = (ssnow%GWwb(i)-soil%GWwatr(i))/(soil%GWssat_vec(i)-soil%GWwatr(i))
       s_mid(i) = min(max(s_mid(i),0.001_r_2),1._r_2)
       s2(i)    = soil%GWhyds_vec(i)*s_mid(i)**(2._r_2*soil%bch_vec(i,ms)+2._r_2)

       ssnow%GWhk(i)     =s_mid(i)*s2(i)*(1._r_2-ssnow%fracice(i,ms))*&
                          exp(-gw_params%hkrz*(zimm(ms)/1000._r_2-gw_params%zdepth))

       ssnow%GWdhkdw(i)  = (1._r_2-ssnow%fracice(i,ms))* (2._r_2*soil%bch_vec(i,ms)+3._r_2)*&
                           s2(i)*0.5_r_2/(soil%GWssat_vec(i)-soil%GWwatr(i))*&
                           exp(-gw_params%hkrz*(zimm(ms)/1000._r_2-gw_params%zdepth))

       ssnow%GWsmp(i)    = -soil%sucs_vec(i,ms)*s_mid(i)**(-soil%bch_vec(i,ms))
       ssnow%GWsmp(i)    = max(min(ssnow%GWsmp(i),-soil%GWsucs_vec(i)),sucmin)
       ssnow%GWdsmpdw(i) = -soil%GWbch_vec(i)*ssnow%GWsmp(i)/(s_mid(i)*(soil%GWssat_vec(i)-soil%GWwatr(i)))
    end do

END SUBROUTINE calc_soil_hydraulic_props


  SUBROUTINE aquifer_recharge(ssnow,soil,veg,zaq,zmm,dzmm)
  USE cable_common_module

  IMPLICIT NONE
  
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssnow ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(IN)    :: soil  ! soil parameters
    TYPE (veg_parameter_type), INTENT(IN)     :: veg
    REAL(r_2), dimension(:), intent(in)       :: zaq
    REAL(r_2), dimension(:), intent(in)       :: zmm,dzmm

    integer :: i    

    !Doing the recharge outside of the soln of Richards Equation makes it easier to track total recharge amount.
    !Add to ssnow at some point 
    do i=1,mp
       if ((ssnow%wtd(i) .le. sum(dzmm,dim=1)) .or. (veg%iveg(i) .ge. 16) .or. (soil%isoilm(i) .eq. 9))  then
          ssnow%Qrecharge(i) = 0._r_2
       else
          ssnow%Qrecharge(i) = -ssnow%hk(i,ms)*((ssnow%GWsmp(i)-ssnow%smp(i,ms)) - (ssnow%GWzq(i)-ssnow%zq(i,ms)))/(zaq(i) - zmm(ms))
       end if
    end do


  END SUBROUTINE aquifer_recharge


END MODULE cable_soil_snow_gw_module
