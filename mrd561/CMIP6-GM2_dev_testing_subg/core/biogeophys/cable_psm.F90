MODULE cable_psm

   USE cable_def_types_mod, only : r_2,ms,mp,air_type,met_type,soil_snow_type,&
                                  canopy_type,soil_parameter_type,veg_parameter_type,&
                                  roughness_type
   USE cable_common_module, only : cable_user

implicit none

   PRIVATE
   REAL(r_2), parameter :: Dff=2.5e-5, &  !diffusivity water vapor in air
                      lm=1.73e-5, &       !converts units
                      pi = 3.14159265358979324, &  !obvous
                      c2 = 2.0,&                  !params
                      subg_cover_Dsense=1.25e-4  !param based on vh thermal diffusivity

   real(r_2), parameter :: rtevap_max = 10000.0

   PUBLIC  or_soil_evap_resistance,update_or_soil_resis

contains

  recursive function my_gamma(a) result(g)

   
    real(r_2), intent(in) :: a 
    real(r_2) :: g 

    real(r_2), parameter :: pi = 3.14159265358979324
    integer, parameter :: cg = 7

    ! these precomputed values are taken by the sample code in Wikipedia,
    ! and the sample itself takes them from the GNU Scientific Library
    real(r_2), dimension(0:8), parameter :: p = &
         (/ 0.99999999999980993, 676.5203681218851, -1259.1392167224028, &
         771.32342877765313, -176.61502916214059, 12.507343278686905, &
         -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 /)

    real(r_2) :: t, w, x 
    integer :: i 

    x = a

    if ( x < 0.5 ) then 
       g = pi / ( sin(pi*x) * my_gamma(1.0-x) )
    else 
       x = x - 1.0
       t = p(0) 
       do i=1, cg+2 
          t = t + p(i-1)/(x+real(i,r_2))
       end do
       w = x + real(cg,r_2) + 0.5
       g = sqrt(2.0*pi) * w**(x+0.5) * exp(-w) * t
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



   REAL(r_2), DIMENSION(mp) :: sublayer_dz, eddy_shape,eddy_mod,soil_moisture_mod, &
                          soil_moisture_mod_sat, wb_liq, &
                          pore_size,pore_radius, rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m

   REAL(r_2), DIMENSION(mp) :: litter_dz

   INTEGER, DIMENSION(mp) :: int_eddy_shape

   !

   integer :: i,j,k 

   if (cable_user%litter) then
      litter_dz(:) = veg%clitt*0.003
   else
      litter_dz(:) = 0.0
   endif

   pore_radius(:) = 0.148  / (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
   pore_size(:) = pore_radius(:)*sqrt(pi)

      !scale ustar according to the exponential wind profile, assuming we are a mm from the surface
      eddy_shape = 0.3*met%ua/ max(1.0e-2,canopy%us*exp(-rough%coexp*(1.0-max(canopy%sublayer_dz,1.0e-6)/max(1e-2,rough%hruff))))
      int_eddy_shape = floor(eddy_shape)
      eddy_mod(:) = 0.0
      do i=1,mp  
        if (veg%iveg(i) .lt. 16) then 
         eddy_mod(i) = 2.2*sqrt(112.0*pi) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))

         if (int_eddy_shape(i) .gt. 0) then
            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
            do k=1,int_eddy_shape(i)
               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
            end do
         end if
         canopy%sublayer_dz(i) = max(eddy_mod(i) * air%visc(i) / max(1.0e-2,canopy%us(i)*&
                           exp(-rough%coexp(i)*(1.0-canopy%sublayer_dz(i)/max(1e-2,rough%hruff(i))))),1e-7) 

         else

         canopy%sublayer_dz(i) = 0.0

         end if

      end do 

   do i=1,mp
     if (veg%iveg(i) .lt. 16 ) then
 
       !when snow covered soil evap from pores not involved 
       if (ssnow%isflag(i) .eq. 1 .or. ssnow%snowd(i) .ge. 1.0e-4) then
          hk_zero(i) = 1.0e15       !make large so term -> 0, soil not important
          hk_zero_sat(i) = 1.0e15   
          soil_moisture_mod(i) = 0.    !snow covered so not soil moisture limited
          soil_moisture_mod_sat(i) = 0.
          canopy%sublayer_dz(i) = 0.

       else  !no snow so pore evap considered

          wb_liq(i) = real(max(0.0001,min(pi/4.0, &
                       (ssnow%wb(i,1)-ssnow%wbice(i,1) - ssnow%satfrac(i)*soil%ssat_vec(i,1)) / &
                       max((1._r_2 - ssnow%satfrac(i)),1e-5) ) ) )
       
          rel_s(i) = real( max(wb_liq(i)-soil%watr(i,1),0._r_2)/(soil%ssat_vec(i,1)-soil%watr(i,1)) )
          hk_zero(i) = max(0.001*soil%hyds_vec(i,1)*&
                       (min(max(rel_s(i),0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(i,1)+3._r_2) ),1e-8)
          hk_zero_sat(i) = max(0.001*soil%hyds_vec(i,1),1e-8)
       
          soil_moisture_mod(i)     = 1.0/pi/sqrt(wb_liq(i))* ( sqrt(pi/(4.0*wb_liq(i)))-1.0)
          soil_moisture_mod_sat(i) = 1.0/pi/sqrt(soil%ssat_vec(i,1))* ( sqrt(pi/(4.0*soil%ssat_vec(i,1)))-1.0)
   
          canopy%sublayer_dz(i) = canopy%sublayer_dz(i) + litter_dz(i)

       end if
    
    
       if (canopy%sublayer_dz(i) .ge. 1.0e-7) then
          ssnow%rtevap_unsat(i) = min(rtevap_max, &
                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero(i)) +&
                                   (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / Dff))
          ssnow%rtevap_sat(i)  = min(rtevap_max, &
                                   rough%z0soil(i)/canopy%sublayer_dz(i) * (lm/ (4.0*hk_zero_sat(i)) + &
                                  (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / Dff))
    
          ssnow%rt_qh_sublayer(i) = canopy%sublayer_dz(i) / subg_cover_Dsense
    
       else
          ssnow%rtevap_unsat(i) = min(rtevap_max, &
                               lm/ (4.0*hk_zero(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod(i)) / Dff)
          ssnow%rtevap_sat(i)  = min(rtevap_max, &
                             lm/ (4.0*hk_zero_sat(i)) + (canopy%sublayer_dz(i) + pore_size(i) * soil_moisture_mod_sat(i)) / Dff)
    
          ssnow%rt_qh_sublayer(i) = 0.0
       end if

     else
     !no additional evap resistane over lakes
        ssnow%rtevap_sat(i) = 0.0
        ssnow%rtevap_unsat(i) = 0.0
        ssnow%rt_qh_sublayer(i) = 0.0
        ssnow%satfrac(i) = 0.5
        !keep the same evap limits at defaul evap scheme over lakes
        !line below matches total resistance to equal defauult beta form
        if (veg%iveg(i) .eq. 16 .and. met%tk(i) .lt. 268.15 ) ssnow%rtevap_sat(i) = 0.41*ssnow%rtsoil(i)
  
    end if

  end do


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
            ssnow%rtevap_sat(i) = min(rtevap_max,canopy%sublayer_dz(i)/Dff)
         end if

         if (dqu(i) .le. 0.0) THEN
            ssnow%rtevap_unsat(i) = min(rtevap_max,canopy%sublayer_dz(i)/Dff)
         end if

      end if
   end do






END SUBROUTINE update_or_soil_resis

END MODULE cable_psm

