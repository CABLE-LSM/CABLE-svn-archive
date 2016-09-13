MODULE large_scale_hydro

  USE sli_numbers,        ONLY:  thousand, params, rhow, rhoi,  &
       freezefac, topmodel, alpha, fsat_max, botbc
  USE sli_utils,          ONLY: dx, x

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp           

  USE cable_common_module, ONLY : gw_params,cable_user
 
  IMPLICIT NONE

  REAL(r_2), PARAMETER ::  wtd_max      = 1000000.0,    &! maximum wtd [mm]
                           wtd_min      = 0.0           ! minimum wtd [mm]


  !global variables need to be wtd,q_through,q_surface,sat_fraction,sat_ice_fraction
  PUBLIC overland_runoff, diagnose_watertable_depth, determine_subsurface_runoff,&
         aquifer_recharge

  contains

  SUBROUTINE overland_runoff(parin,ssnow,soil,veg)
     !must be called after hofS is called as that sets thetai or call hofs again here...
     !USE ==> rhoi,vars (var),params, thousand
     IMPLICIT NONE 
     TYPE(params), DIMENSION(:,:),            INTENT(IN)    :: parin         !soil  parameters
     TYPE(soil_snow_type),                    INTENT(INOUT) :: ssnow           !prog vars
     TYPE(soil_parameter_type), DIMENSION(:), INTENT(IN)    :: soil
     TYPE(veg_parameter_type),  DIMENSION(:), INTENT(IN)    :: veg
     !LOCAL
     REAL(r_2), DIMENSION(mp)    :: ice_frac
     REAL(r_2), DIMENSION(mp)    :: ice_mass
     REAL(r_2), DIMENSION(mp)    :: liq_mass
     REAL(r_2), DIMENSION(mp)    :: ice_factor
     REAL(r_2), DIMENSION(mp)    :: rel_S
     REAL(r_2), DIMENSION(mp)    :: effective_porosity

     REAL(r_2) :: slopeSTDmm
     REAL(r_2) :: tmpa,tmpb,qinmax
     integer :: i,k

     do i=1,mp
        effective_porosity(i) = max(0.001_r_2, parin(i,1)%the - ssnow%thetai(i,1))
        ice_mass(i) = ssnow%thetai(i,1) * dx(i,1) * rhoi
        liq_mass(i) = (ssnow%S(i,1)*(parin(i,1)%the-parin(i,1)%thr)-ssnow%thetai(i,1)) * dx(i,1) * rhow
     end do
     ice_frac(:)     = max(0._r_2,min(1._r_2,ice_mass(:)/ (ice_mass(:) + liq_mass(:))))

     rel_S(:) = 0._r_2
     do k=1,ms
       rel_S(:) = rel_S(:) + max(0.001,min(1.0,ssnow%S(:,k) ))*dx(:,k)
     end do
     rel_S(:) = rel_S(:)/sum(dx(:,1:ms),dim=2)

     !srf frozen fraction.  should be based on topography
     do i = 1,mp

      ice_factor(i) = (exp(gw_params%IceAlpha*(1._r_2-ice_frac(i)))-exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
      ice_factor(i)  = min(max(ice_factor(i),0._r_2),1._r_2)
      !Saturated fraction

       if (gw_params%MaxSatFraction .gt. 1e-7) then 
          slopeSTDmm = sqrt(max(real(gw_params%MaxSatFraction,r_2)*soil%slope_std(i),1e-5)) ! ensure some variability
          ssnow%satfrac(i)    = max(1e-6,min(0.95,1._r_2 - erf( slopeSTDmm / sqrt(2.0* rel_S(i)) ) ) )  
       else
          ssnow%satfrac(i) = 0. 
       end if
       ssnow%sat_ice_frac(i)   = max(1e-6,min(0.95,ice_factor(i) + (1._r_2-ice_factor(i))*ssnow%satfrac(i) ) )
     end do

     do i=1,mp
        tmpa = (ssnow%S(i,1)*(parin(i,1)%the-parin(i,1)%thr)-ssnow%thetai(i,1)) / effective_porosity(i)
        tmpb = max( (tmpa-ssnow%sat_ice_frac(i))/max(0.01_r_2,(1._r_2-ssnow%sat_ice_frac(i))), 0._r_2)
        tmpa = -2._r_2*abs(parin(i,1)%he)/(parin%lam(i,1)*dx(1))
        qinmax = (1._r_2 + tmpa*(tmpb-1._r_2))*par(i,1)%Ke*thousand !IN mm/s 

        ssnow%rnof1(i) = ssnow%sat_ice_frac(i) * ssnow%fwtop(i) + &
                         (1._r_2-ssnow%sat_ice_frac(i))*max((ssnow%fwtop(i)-qinmax) , 0._r_2)

        ssnow%fwtop(i) = ssnow%fwtop(i) - ssnow%rnof1(i)

     end do  !mp

     !add back to the lakes to keep saturated instead of drying
     !imporant for ACCESS
     where (veg%iveg .eq. 16)
         ssnow%fwtop(:) = ssnow%fwtop(:) + ssnow%rnof1(:)
         ssnow%fwtop(:) = 0._r_2
     endwhere
     


  END SUBROUTINE overland_runoff



  SUBROUTINE diagnose_watertable_depth(ssnow,soil,parin,veg)
    IMPLICIT NONE

     TYPE(soil_snow_type),          INTENT(INOUT)             :: ssnow            !degree of saturation
     TYPE(soil_parameter_type),     INTENT(IN)                :: soil            !degree of saturation
     TYPE(params), DIMENSION(:,:),  INTENT(IN)                :: parin            !soil  parameters
     TYPE (veg_parameter_type),     INTENT(IN)                :: veg
 
     !Local vars 
     REAL(r_2), DIMENSION(mp,ms)   :: dzmm_mp,tmp_def
     REAL(r_2), DIMENSION(mp,0:ms) :: zimm
     REAL(r_2), DIMENSION(ms)      :: zmm
     REAL(r_2), DIMENSION(mp)      :: GWzimm,temp
     REAL(r_2), DIMENSION(mp)      :: def,defc     

     REAL(r_2)                     :: deffunc,tempa,tempb,derv,calc,tmpc
     REAL(r_2), DIMENSION(mp)      :: Nsmpsat  !ensure this is positive
     REAL(r_2), DIMENSION(mp,ms)   :: theta_total
     INTEGER :: k,i,wttd,jlp

     !make code cleaner define these here 
     do i=1,mp
        Nsmpsat(i)   = abs(parin(i,ms)%he)*thousand                                !psi_saturated mm
     end do
     dzmm_mp   = dx * thousand
     zimm(:,0) = 0._r_2 !                                       !depth of layer interfaces mm

     do k=1,ms
       zimm(:,k) = zimm(:,k-1) + dx(:,k)*1000._r_2
     end do
     zimm(:,ms) = zimm(:,ms) +  soil%GWdz(:)*1000._r_2

     do k=1,ms
        do i=1,mp
           theta_total(i,k) = ssnow%S(i,k)*(parin(i,k)%the-parin(i,k)%thr) + parin(i,k)%thr
        end do
     end do
     !find the deficit if the water table is at the bottom of the soil column
     do i=1,mp
        defc(i) = (parin(i,ms)%the)*(zimm(i,ms)+Nsmpsat(i)/(1._r_2-parin(i,ms)%lam)*            &
                (1._r_2-((Nsmpsat(i)+zimm(i,ms))/Nsmpsat(i))**(1._r_2-parin(i,ms)%lam))) 
        defc(i) = max(0.1_r_2,defc(i)) 
     end do

     def(:) = 0._r_2
     do k=1,ms
        do i=1,mp

           if (parin(i,k)%the .gt. theta_total(i,k)) then
              def(i) = def(i) +                                                           &
                     (parin(i,k)%the-theta_total(i,k))*dzmm_mp(i,k)
           end if
        end do  !mp
     end do  !ms

     do i=1,mp
        def(i) = def(i) + max(0._r_2,soil%GWwatsat(i)-ssnow%GWwb(i))*soil%GWdz(i)*1000._r_2
     end do   

     ssnow%wtd(:) = zimm(:,ms)*def(:)/defc(:)

     do i=1,mp

        if ((veg%iveg(i) .ne. 16) .and. (soil%isoilm(i) .ne. 9)) then      

           if (defc(i) > def(i)) then                 !iterate tfor wtd

              jlp=0

              mainloop: DO

                 tempa   = 1.0_r_2
                 tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-parin(i,ms)%lam)
                 derv    = (parin(i,ms)%the)*(tempa-tempb) + &
                                       parin(i,ms)%the

                 if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

                 tempa   = 1.0_r_2
                 tempb   = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-parin(i,ms)%lam)
                 deffunc = (parin(i,ms)%the)*(ssnow%wtd(i) +&
                           Nsmpsat(i)/(1._r_2-parin(i,ms)%lam)* &
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

                 tmpc     = Nsmpsat(i)+ssnow%wtd(i)-zimm(i,ms)
                 tempa    = (abs(tmpc/Nsmpsat(i)))**(-parin(i,ms)%lam)
                 tempb    = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(-parin(i,ms)%lam)
                 derv     = (parin(i,ms)%the)*(tempa-tempb)
                 if (abs(derv) .lt. real(1e-8,r_2)) derv = sign(real(1e-8,r_2),derv)

                 tempa    = (abs((Nsmpsat(i)+ssnow%wtd(i)-zimm(i,ms))/Nsmpsat(i)))**(1._r_2-parin(i,ms)%lam)
                 tempb    = (1._r_2+ssnow%wtd(i)/Nsmpsat(i))**(1._r_2-parin(i,ms)%lam)
                 deffunc  = (parin(i,ms)%the)*(zimm(i,ms) +&
                             Nsmpsat(i)/(1._r_2-parin(i,ms)%lam)*(tempa-tempb))-def(i)
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

              ssnow%wtd(i) = zimm(i,ms)

           endif

       endif  !check veg and soils

     end do   !mp loop


     ssnow%wtd(:) = ssnow%wtd(:)/thousand
     !limit wtd to be within a psecified range
     where (ssnow%wtd(:) .gt. wtd_max) ssnow%wtd(:) = wtd_max
     where (ssnow%wtd(:) .lt. wtd_min) ssnow%wtd(:) = wtd_min



  END SUBROUTINE diagnose_watertable_depth



  SUBROUTINE determine_subsurface_runoff(parin,ssnow,soil,veg)

     IMPLICIT NONE 
     TYPE(params), DIMENSION(:,:),            INTENT(IN)     :: parin            !soil  parameters
     TYPE(soil_snow_type), DIMENSION(:),      INTENT(INOUT)  :: ssnow           !prog vars
     TYPE(soil_parameter_type), DIMENSION(:), INTENT(IN)     :: soil
     TYPE(veg_parameter_type),  DIMENSION(:), INTENT(IN)     :: veg
     !LOCAL
     REAL(r_2)                  :: ice_frac

     REAL(r_2) :: q_sub_total
     REAL(r_2) :: tmpa,tmpb,qinmax
     INTEGER :: k_drain,k,i



     do i=1,mp

        ssnow%qhz(i)  = max(tan(soil%slope(i)),0.001) *gw_params%MaxHorzDrainRate * &
                    exp(-ssnow%wtd(i)/gw_params%EfoldHorzDrainRate)

        !Keep "lakes" saturated forcing qhz = 0.  runoff only from lakes
        !overflowing
        if ((soil%isoilm(i) .eq. 9) .or. (veg%iveg(i) .eq. 16)) then
           ssnow%qhz(i) = 0._r_2
        end if
 
        !identify first no frozen layer.  drinage from that layer and below
        !drain from sat layers
        k_drain = ms+1
        do k=ms,2,-1
           if (ssnow%wtd(i) .le. sum(dx(i,1:k),dim=1)) then
               k_drain = k
           end if
        end do
        k_drain = max(k_drain,2)

        ssnow%qhlev(i,:) = 0._r_2
        sm_tot(i) = 0._r_2
        if (k_drain .le. ms) then
           do k=k_drain,ms
              sm_tot(i) = sm_tot(i) + max(ssnow%S(i,k)*(parin(i,k)%the-parin(i,k)%thr)+parin(i,k)%thr,0._r_2)
           end do
           sm_tot(i) = max(sm_tot(i),0.01_r_2)

           do k=k_drain,ms
              ssnow%qhlev(i,k) = q_sub_total(i)*max(ssnow%S(i,k)*(parin(i,k)%the-parin(i,k)%thr)+parin(i,k)%thr,0._r_2)/sm_tot(i)
           end do

         else
            ice_frac = (exp(gw_params%IceAlpha*(1._r_2-(ssnow%thetai(i,ms)-parin(i,ms)%thr)/(parin(i,ms)%the-parin(i,ms)%thr))) &
                        -exp(gw_params%IceAlpha))/(1._r_2-exp(gw_params%IceAlpha))
            ssnow%qhlev(i,ms+1) = max(1._r_2-ice_frac,0._r_2)*ssnow%qhz(i)
         end if

         !incase every layer is frozen very dry
         ssnow%qhz(i) = ssnow%qhlev(i,ms+1)
         do k=k_drain,ms
            ssnow%qhz(i) = ssnow%qhz(i) + ssnow%qhlev(i,k)
         end do

    end do  


  END SUBROUTINE determine_subsurface_runoff


  SUBROUTINE aquifer_recharge(delt,ssnow,parin,veg,soil)
     REAL(r_2),                               INTENT(IN)     :: delt
     TYPE(params), DIMENSION(:,:),            INTENT(IN)     :: parin            !soil  parameters
     TYPE(soil_snow_TYPE),                    INTENT(INOUT)  :: ssnow
     TYPE(soil_parameter_type),               INTENT(IN)     :: soil
     TYPE(veg_parameter_type),  DIMENSION(:), INTENT(IN)     :: veg

     REAL(r_2) :: mass_needed
     INTEGER :: i,k
  
     do i=1,mp
        if ((ssnow%wtd(i) .ge. sum(dx(i,:),dim=1)) .and. (veg%iveg(i) .ne. 17) .and. (soil%isoilm(i) .ne. 9))  then

           ssnow%q_recharge(i) = -soil%GW_Kaq*((soil%GW_he-var(i,ms)%h) - (ssnow%wtd(i) - x(i,ms)))/(ssnow%wtd(i) - x(i,ms))
        else

           ssnow%q_recharge(i) = 0._r_2
        end if
     end do

     !limit in case we don't have enough mass of liquid water in bottom soil layer
     mass_needed = ssnow%q_recharge(i)*delt
     if (ssnow%S(i,ms) .lt. (parin(i,ms)%the-parin(i,ms)%thr)*mass_needed/dx(i,ms)) then

        ssnow%q_recharge(i) = 0.9 * ssnow%S(i,ms)*(parin(i,ms)%the-parin(i,ms)%thr)*dx(i,ms)/delt

     end if

     ssnow%GWwb(i) = ssnow%GWwb(i) + (ssnow%q_recharge(i) -ssnow%qhlev(i,ms+1) )*delt/soil%GWdz(i)

     ssnow%S(i,ms) = ssnow%S(i,ms) - (parin(i,ms)%the-parin(i,ms)%thr)*(ssnow%q_recharge(i))*delt/dx(i,ms)

     !do we have too much water? add to runoff
     if (ssnow%GWwb(i) .gt. soil%GWwatsat(i)) then
        ssnow%qhlev(i,ms+1) = ssnow%qhlev(i,ms+1) + (soil%GWwatsat(i)-ssnow%GWwb(i))*soil%GWdz(i)/delt

        ssnow%GWwb(i) = soil%GWwatsat(i)

     end if



  END SUBROUTINE aquifer_recharge








END MODULE large_scale_hydro
