MODULE cable_psm

   USE cable_def_types_mod, only : r_2,ms,mp,air_type,met_type,soil_snow_type,&
                                  canopy_type,soil_parameter_type,veg_parameter_type,&
                                  roughness_type
   USE cable_common_module, only : cable_user
   use cable_io_vars_module, only: wlogn

implicit none


   REAL(r_2), parameter ::rt_Dff=2.5e-5, & !diffusivity in air
                      lm=1.73e-5, &       !converts units
                      !lm= 1.73e-9,&        ! MMY
                      c2 = 2.0,&                  !params
                      litter_thermal_diff=2.7e-5  !param based on vh thermal diffusivity

   real(r_2), parameter :: rtevap_max = 10000.0
   ! these precomputed values are taken by the sample code in Wikipedia,
   ! and the sample itself takes them from the GNU Scientific Library
   REAL(r_2), DIMENSION(0:8), parameter :: gamma_pre = &
         (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
         771.32342877765313, -176.61502916214059, 12.507343278686905, &
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)
   INTEGER, PARAMETER                      :: c_gamma = 7
    real(r_2),parameter :: pi_r_2=3.14159

PUBLIC  or_soil_evap_resistance,update_or_soil_resis,rtevap_max,rt_Dff

contains

  recursive function my_gamma(a) result(g)


    real(r_2), intent(in) :: a
    real(r_2) :: g

    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    !real(r_2), dimension(0:8), parameter :: gamma_pre = &
    !     (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
    !     771.32342877765313, -176.61502916214059, 12.507343278686905, &
    !     -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)

    real(r_2) :: t, w, x
    integer :: i

    x = a

    if ( x < 0.5 ) then
       g = (pi_r_2) / ( sin((pi_r_2)*x) * my_gamma(1.0-x) )
    else
       x = x - 1.0
       t = gamma_pre(0)
       do i=1, c_gamma+2
          t = t + gamma_pre(i-1)/(x+real(i,r_2))
       end do
       w = x + real(c_gamma,r_2) + 0.5
       g = sqrt(2.0*(pi_r_2)) * w**(x+0.5) * exp(-w) * t
    end if
  end function my_gamma

SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough)

   TYPE (air_type), INTENT(INOUT)       :: air
   TYPE (met_type), INTENT(INOUT)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy
   TYPE (soil_parameter_type), INTENT(INOUT)   :: soil
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   TYPE (roughness_type), INTENT(INOUT) :: rough



   REAL(r_2), DIMENSION(mp) :: eddy_shape,eddy_mod,soil_moisture_mod, &
                          soil_moisture_mod_sat, wb_liq, &
                          pore_size,pore_radius, rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m

   REAL(r_2), DIMENSION(mp) :: litter_dz

   INTEGER, DIMENSION(mp) :: int_eddy_shape

   !

   integer :: i,j,k

   canopy%sublayer_dz(:) = 0.005

   litter_dz(:) = 0.0
   if (cable_user%litter) then
      where (ssnow%isflag .eq. 0 .or. ssnow%snowd .le. 0.1)
         litter_dz(:) = veg%clitt*0.003
      endwhere
   endif

   pore_radius(:) = 0.148  / (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle ! MMY
!   pore_radius(:) = 0.148*0.707 / (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0) !should replace 0.148 with surface tension, unit coversion, and angle ! MMY

   pore_size(:) = pore_radius(:)*sqrt((pi_r_2))

   !scale ustar according to the exponential wind profile, assuming we are a mm from the surface

   eddy_mod(:) = 0.0
   eddy_shape(:) = 1.0
   int_eddy_shape(:) = 0

   do i=1,mp

     if (veg%iveg(i) .lt. 16) then

        eddy_shape(i) = 0.3*met%ua(i)/ max(1.0e-4,max(1.0e-3,canopy%us(i))*&
                        exp(-rough%coexp(i)*(1.0-canopy%sublayer_dz(i)/max(1e-2,rough%hruff(i)))))

        int_eddy_shape(i) = floor(eddy_shape(i))

        eddy_mod(i) = 2.2*sqrt(112.0*(pi_r_2)) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))

         if (int_eddy_shape(i) .gt. 0) then
            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
            do k=1,int_eddy_shape(i)
               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
            end do
         end if
         canopy%sublayer_dz(i) = min(0.05, max(eddy_mod(i) * air%visc(i) / max(1.0e-3,canopy%us(i))*&
                           exp(-rough%coexp(i)*(1.0-canopy%sublayer_dz(i)/max(1e-2,rough%hruff(i)))),1e-7)  )

      else

         canopy%sublayer_dz(i) = 0.0

      end if

   end do

   do i=1,mp
     if (veg%iveg(i) .lt. 16) then

       ! ________________________________ MMY _______________________________
       wb_liq(i) = real(max(0.0001,min((pi_r_2)/4.0, &
                    (ssnow%wb(i,1)-ssnow%wbice(i,1) - ssnow%satfrac(i)*soil%ssat_vec(i,1)) / &
                    max((1._r_2 - ssnow%satfrac(i)),1e-5) ) ) )
       !print *, "MMY ssnow%satfrac(i) is ", ssnow%satfrac(i)
       !wb_liq(i) = real(max(0.0001,min((pi_r_2)/4.0, &
       !             ((ssnow%wb(i,1)*0.022 + ssnow%wb(i,2)*0.058 + &
       !               ssnow%wb(i,3)*0.154 + ssnow%wb(i,4)*0.266)/0.5 &
       !             -(ssnow%wbice(i,1)*0.022 + ssnow%wbice(i,2)*0.058 + &
       !               ssnow%wbice(i,3)*0.154 + ssnow%wbice(i,4)*0.266)/0.5 &
       !             -ssnow%satfrac(i)*soil%ssat_vec(i,1)) / &
       !             max((1._r_2 - ssnow%satfrac(i)),1e-5) ) ) )
       !print *, "MMY top 50cm wb_liq(i) is ", wb_liq(i)
       ! ____________________________________________________________________
       rel_s(i) = real( max(wb_liq(i)-soil%watr(i,1),0._r_2)/(soil%ssat_vec(i,1)-soil%watr(i,1)) )
       !hk_zero(i) = max(0.001*soil%hyds_vec(i,1)*(min(max(rel_s(i),0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(i,1)+3._r_2) ),1e-12)
       !hk_zero_sat(i) = max(0.001*soil%hyds_vec(i,1),1e-12)
         
       hk_zero(i) = max(0.001*soil%hyds_vec(i,1)*(min(max(rel_s(i),0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(i,1)+3._r_2)),1e-8)
       hk_zero_sat(i) = max(0.001*soil%hyds_vec(i,1),1e-8)

       !print *, "MMY top 2cm soil%watr(i,1)  is ", soil%watr(i,1)
       !print *, "MMY top 2cm soil%ssat_vec(i,1) is ", soil%ssat_vec(i,1)
       !print *, "MMY top 2cm soil%hyds_vec(i,1) is ", soil%hyds_vec(i,1)
       ! _____________________________ MMY ADD ______________________________
       !wb_liq(i) = real(max(0.0001,min((pi_r_2)/4.0, &
       !             (ssnow%wb(i,1)-ssnow%wbice(i,1) -ssnow%satfrac(i)*soil%ssat_vec(i,1)) / &
       !             max((1._r_2 - ssnow%satfrac(i)),1e-5) ) ) )
       !rel_s(i) = real(max(wb_liq(i)-soil%watr(i,1),0._r_2)/(soil%ssat_vec(i,1)-soil%watr(i,1)) )
       !print *, "MMY top 2cm wb_liq(i) is ", wb_liq(i)
       !print *, "MMY top 2cm rel_s(i) is", rel_s(i)
       ! ____________________________________________________________________

       soil_moisture_mod(i)     = 1.0/(pi_r_2)/sqrt(wb_liq(i))* ( sqrt((pi_r_2)/(4.0*wb_liq(i)))-1.0)
       soil_moisture_mod_sat(i) = 1.0/(pi_r_2)/sqrt(soil%ssat_vec(i,1))* ( sqrt((pi_r_2)/(4.0*soil%ssat_vec(i,1)))-1.0)

       if (ssnow%isflag(i) .eq. 1 .or. (ssnow%snowd(i) .gt. 0.1) ) then

          hk_zero(i)     = 1.0e15
          hk_zero_sat(i) = 1.0e15
          soil_moisture_mod(i)     = 0.0
          soil_moisture_mod_sat(i) = 0.0

       end if

       canopy%sublayer_dz(i) = canopy%sublayer_dz(i) + litter_dz(i)

       if (canopy%sublayer_dz(i) .ge. 1.0e-7 .and. hk_zero(i) .lt. 1.0e14) then
! _____________________________________ MMY ___________________________________
          ssnow%rtevap_unsat(i) = min(rtevap_max, &
                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero(i)) +&
                                   (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / rt_Dff))
          ssnow%rtevap_sat(i)  = min(rtevap_max, &
                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero_sat(i)) + &
                                  (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / rt_Dff))

!          ssnow%rtevap_unsat(i) = min(rtevap_max, &
!                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero(i))))
!          ssnow%rtevap_sat(i)   = min(rtevap_max, &
!                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/(4.0*hk_zero_sat(i))))

! _____________________________________________________________________________
          ssnow%rt_qh_sublayer(i) = canopy%sublayer_dz(i) / litter_thermal_diff

       else
! _____________________________________ MMY  ________________________________________
          ssnow%rtevap_unsat(i) = min(rtevap_max, &
                               lm/ (4.0*hk_zero(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / rt_Dff)
          ssnow%rtevap_sat(i)  = min(rtevap_max, &
                             lm/ (4.0*hk_zero_sat(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / rt_Dff)
!          ssnow%rtevap_unsat(i) = min(rtevap_max,  lm/ (4.0*hk_zero(i)))
!          ssnow%rtevap_sat(i)   = min(rtevap_max, lm/ (4.0*hk_zero_sat(i)))
! ___________________________________________________________________________________
          ssnow%rt_qh_sublayer(i) = canopy%sublayer_dz(i) / litter_thermal_diff
       end if

     else
     !no additional evap resistane over lakes
        ssnow%rtevap_unsat(i) = 0.0
        ssnow%rt_qh_sublayer(i) = 0.0
        ssnow%satfrac(i) = 0.5
        if (veg%iveg(i) .eq. 16 .and. met%tk(i) .lt. 268.15 ) &
              ssnow%rtevap_sat(i) = 0.41*ssnow%rtsoil(i)

    end if

  end do

  print *, "# ",wb_liq," ", hk_zero ! MMY
  print *, "+ ",wb_liq," ", lm/ (4.0*hk_zero)  ! MMY
  print *, "MMY rel_s is ", rel_s  ! MMY
  print *, "MMY hk_zero is ", hk_zero ! MMY
  print *, "MMY hk_zero_sat is ", hk_zero_sat ! MMY
  print *, "MMY unsat rsv is ", lm/ (4.0*hk_zero) ! MMY
  print *, "MMY sat rsv is ", lm/ (4.0*hk_zero_sat) ! MMY
  print *, "MMY unsat rBL is ",(canopy%sublayer_dz+pore_size*soil_moisture_mod)/rt_Dff ! MMY
  print *, "MMY sat rBL is ",(canopy%sublayer_dz+pore_size*soil_moisture_mod_sat)/rt_Dff ! MMY
  print *, "MMY ssnow%rtevap_unsat is ", ssnow%rtevap_unsat ! MMY
  print *, "MMY ssnow%rtevap_sat is ", ssnow%rtevap_sat ! MMY
  print *, "MMY canopy%sublayer_dz is ", canopy%sublayer_dz ! MMY
END SUBROUTINE or_soil_evap_resistance


SUBROUTINE update_or_soil_resis(ssnow,canopy,veg,dq,dqu)

   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy
   TYPE (veg_parameter_type), INTENT(INOUT) :: veg
   REAL, DIMENSION(mp), INTENT(IN) :: dq,&
                                      dqu

   INTEGER :: i


   do i=1,mp

      if (veg%iveg(i) .lt. 16 .and. ssnow%snowd(i) .lt. 1e-7) THEN

         if (dq(i) .le. 0.0) THEN
            ssnow%rtevap_sat(i) = min(rtevap_max,canopy%sublayer_dz(i)/rt_Dff)
         end if

         if (dqu(i) .le. 0.0) THEN
            ssnow%rtevap_unsat(i) = min(rtevap_max,canopy%sublayer_dz(i)/rt_Dff)
         end if

      end if
   end do






END SUBROUTINE update_or_soil_resis

END MODULE cable_psm
