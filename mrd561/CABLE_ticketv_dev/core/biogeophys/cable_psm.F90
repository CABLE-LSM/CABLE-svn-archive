MODULE cable_psm
USE cable_def_types_mod
implicit none

PUBLIC  or_soil_evap_resistance

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


!SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough,snow_covered,tile_index,dz_litter)
!SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough,snow_covered,dz_litter)
!   USE cable_def_types_mod
!   USE cable_air_module
!   USE cable_common_module
!
!   TYPE (air_type), INTENT(IN)       :: air
!   TYPE (met_type), INTENT(IN)       :: met
!   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
!   TYPE (canopy_type), INTENT(INOUT)    :: canopy
!   TYPE (soil_parameter_type), INTENT(IN)   :: soil
!   TYPE (veg_parameter_type), INTENT(IN) :: veg
!   TYPE (roughness_type), INTENT(IN) :: rough
!   integer, dimension(:), intent(in)  :: snow_covered
!   !integer, intent(in), optional :: tile_index
!   !real(r_2), intent(in), optional :: dz_litter
!   real(r_2), dimension(:), intent(in) :: dz_litter
!
!
!   REAL(r_2), DIMENSION(mp) :: sublayer_dz, eddy_shape,eddy_mod,soil_moisture_mod, &
!                          soil_moisture_mod_sat, wb_liq, &
!                          pore_size,pore_radius,rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m
!
!   INTEGER, DIMENSION(mp) :: int_eddy_shape
!
!   REAL(r_2), parameter :: Dff=2.5e-5, &
!                      lm=1.73e-5, &
!                      pi = 3.14159265358979324, &
!                      c2 = 2.0
!
!   real(r_2), parameter :: rtevap_max = 100000.0
!   real(r_2),dimension(mp)  :: litter_thickness
!
!   integer :: i,j,k,tile_index
!   logical, save ::  first_call = .true.
!
!   !if (present(dz_litter)) then
!      litter_thickness = dz_litter
!   !else
!   !   litter_thickness = 0.0_r_2
!   !end if
!
!!   if (.not.present(tile_index)) then  !called from default cable canopy
!
!      if (first_call)  canopy%sublayer_dz(:) = 0.001
!   
!      pore_radius(:) =0.148/ (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
!      pore_size(:) = pore_radius(:)*sqrt(pi)
!   
!      !scale ustar according to the exponential wind profile, assuming we are a
!      !mm from the surface
!      eddy_shape = 0.3*met%ua/ max(1.0e-4,canopy%us*exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff))))
!      int_eddy_shape = floor(eddy_shape)
!      eddy_mod(:) = 0.0
!      do i=1,mp
!         eddy_mod(i) = 2.2*sqrt(112.0*pi) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))
!   
!         if (int_eddy_shape(i) .gt. 0) then
!            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
!            do k=1,int_eddy_shape(i)
!               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
!            end do
!         end if
!      end do
!      canopy%sublayer_dz = max(eddy_mod(:) * air%visc / max(1.0e-4,canopy%us*exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff)))),1e-7)  + litter_thickness
!   
!      if (first_call) then
!         wb_liq(:) = real(max(1.0e-7,min(pi/4.0, ssnow%wb(:,1)) ) )
!      else
!         wb_liq(:) = real(max(1.0e-7,min(pi/4.0, (ssnow%wb(:,1)-ssnow%wbice(:,1) -ssnow%satfrac(:)*soil%ssat_vec(:,1))/(1._r_2 - ssnow%satfrac(:)) ) ) )
!      end if
!   
!      rel_s = real( max(wb_liq(:)-soil%watr(:,1),0._r_2)/(soil%ssat_vec(:,1)-soil%watr(:,1)) )
!      hk_zero = max(0.001*soil%hyds_vec(:,1)*(min(max(rel_s,0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(:,1)+3._r_2)),1e-8)
!      hk_zero_sat = max(0.001*soil%hyds_vec(:,1),1e-8)
!   
!      soil_moisture_mod(:)     = 1.0/pi/sqrt(wb_liq)* ( sqrt(pi/(4.0*wb_liq))-1.0)
!      soil_moisture_mod_sat(:) = 1.0/pi/sqrt(soil%ssat_vec(:,1))* (sqrt(pi/(4.0*soil%ssat_vec(:,1)))-1.0)
!  
!      !if (present(snow_covered)) then 
!      where (snow_covered .ne. 0)
!            soil_moisture_mod = 0._r_2
!            soil_moisture_mod_sat = 0._r_2
!            canopy%sublayer_dz = canopy%sublayer_dz - litter_thickness
!      end where
!      !end if
!   
!      ssnow%rtevap_unsat(:) = min( (lm/(4.0*hk_zero) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod) / Dff),&
!                            rtevap_max )
!      ssnow%rtevap_sat(:)  =min( (lm/(4.0*hk_zero_sat) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod_sat)/ Dff),&
!                            rtevap_max )
!   
!      !no additional evap resistane over lakes
!      where(veg%iveg .eq. 16)
!         ssnow%rtevap_sat = 0.0
!         ssnow%rtevap_unsat = 0.0
!      endwhere
!   
!   
!!  else  !called point by point from SLI
!!
!!  do tile_index=1,mp
!!   
!!      canopy%sublayer_dz((tile_index)) = 0.001
!!   
!!      pore_radius((tile_index)) =0.148/ (1000.0*9.81*abs(soil%sucs_vec((tile_index),1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
!!      pore_size((tile_index)) = pore_radius((tile_index))*sqrt(pi)
!!   
!!      !scale ustar according to the exponential wind profile, assuming we are a
!!      !mm from the surface
!!      eddy_shape = 0.3*met%ua(tile_index)/ max(1.0e-4,canopy%us(tile_index)*exp(-rough%coexp(tile_index)*(1.0-canopy%sublayer_dz(tile_index)/max(1e-2,rough%hruff(tile_index)))))
!!      int_eddy_shape = floor(eddy_shape)
!!      eddy_mod(:) = 0.0
!!      do i=1,mp
!!         eddy_mod(i) = 2.2*sqrt(112.0*pi) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))
!!   
!!         if (int_eddy_shape(i) .gt. 0) then
!!            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
!!            do k=1,int_eddy_shape(i)
!!               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
!!            end do
!!         end if
!!      end do
!!      canopy%sublayer_dz(tile_index) = max(eddy_mod((tile_index)) * air%visc(tile_index) / max(1.0e-4,canopy%us(tile_index)*exp(-rough%coexp(tile_index)*(1.0-canopy%sublayer_dz(tile_index)/max(1e-2,rough%hruff(tile_index))))),1e-7)  + litter_thickness
!!   
!!      if (first_call) then
!!         wb_liq(tile_index) = real(max(1.0e-7,min(pi/4.0, ssnow%wb(tile_index,1)) ) )
!!      else
!!         wb_liq(tile_index) = real(max(1.0e-7,min(pi/4.0, (ssnow%wb(tile_index,1)-ssnow%wbice(tile_index,1) -ssnow%satfrac(tile_index)*soil%ssat_vec(tile_index,1))/(1._r_2 - ssnow%satfrac(tile_index)) ) ) )
!!      end if
!!   
!!      rel_s(tile_index) = real( max(wb_liq(tile_index)-soil%watr(tile_index,1),0._r_2)/(soil%ssat_vec(tile_index,1)-soil%watr(tile_index,1)) )
!!      hk_zero(tile_index) = max(0.001*soil%hyds_vec(tile_index,1)*(min(max(rel_s(tile_index),0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(tile_index,1)+3._r_2)),1e-8)
!!      hk_zero_sat(tile_index) = max(0.001*soil%hyds_vec(tile_index,1),1e-8)
!!   
!!      soil_moisture_mod(tile_index)     = 1.0/pi/sqrt(wb_liq(tile_index))* ( sqrt(pi/(4.0*wb_liq(tile_index)))-1.0)
!!      soil_moisture_mod_sat(tile_index) = 1.0/pi/sqrt(soil%ssat_vec(tile_index,1))* (sqrt(pi/(4.0*soil%ssat_vec(tile_index,1)))-1.0)
!!   
!!      if (present(snow_covered)) then
!!         if (snow_covered .ne. 0) then
!!
!!            soil_moisture_mod(tile_index) = 0._r_2
!!            soil_moisture_mod_sat(tile_index) = 0._r_2
!!
!!            canopy%sublayer_dz(tile_index) = canopy%sublayer_dz(tile_index) - litter_thickness
!!   
!!         ssnow%rtevap_unsat(tile_index) = min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* &
!!                                            ((canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod(tile_index)) / Dff),&
!!                            rtevap_max )
!!         ssnow%rtevap_sat(tile_index)  =min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* &
!!                                            ((canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod_sat(tile_index))/ Dff),&
!!                            rtevap_max )
!!   
!!         else
!!      ssnow%rtevap_unsat(tile_index) = min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* (lm/(4.0*hk_zero(tile_index)) + (canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod(tile_index)) / Dff),&
!!                            rtevap_max )
!!      ssnow%rtevap_sat(tile_index)  =min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* (lm/(4.0*hk_zero_sat(tile_index)) + (canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod_sat(tile_index))/ Dff),&
!!                            rtevap_max )
!!         end if
!!   
!!      else
!!   
!!         ssnow%rtevap_unsat(tile_index) = min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* (lm/(4.0*hk_zero(tile_index)) + (canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod(tile_index)) / Dff),&
!!                               rtevap_max )
!!         ssnow%rtevap_sat(tile_index)  =min(rough%z0soil(tile_index)/max(canopy%sublayer_dz(tile_index),1e-5)* (lm/(4.0*hk_zero_sat(tile_index)) + (canopy%sublayer_dz(tile_index) + pore_size(tile_index) * soil_moisture_mod_sat(tile_index))/ Dff),&
!!                               rtevap_max )
!!      end if
!!      !no additional evap resistane over lakes
!!      if (veg%iveg(tile_index) == 16) then
!!         ssnow%rtevap_sat(tile_index) = 0.0
!!         ssnow%rtevap_unsat(tile_index) = 0.0
!!      end if
!!
!!   end do
!!
!!   end if
!
!
!   first_call = .false.
!
!END SUBROUTINE or_soil_evap_resistance

SUBROUTINE or_soil_evap_resistance(soil,air,met,canopy,ssnow,veg,rough,snow_covered,dz_litter)
   USE cable_def_types_mod
   USE cable_air_module
   USE cable_common_module   

   TYPE (air_type), INTENT(IN)       :: air
   TYPE (met_type), INTENT(IN)       :: met
   TYPE (soil_snow_type), INTENT(INOUT) :: ssnow
   TYPE (canopy_type), INTENT(INOUT)    :: canopy
   TYPE (soil_parameter_type), INTENT(IN)   :: soil 
   TYPE (veg_parameter_type), INTENT(IN) :: veg
   TYPE (roughness_type), INTENT(IN) :: rough
   integer, dimension(:), intent(in)  :: snow_covered
   real(r_2), dimension(:), intent(in) :: dz_litter



   REAL(r_2), DIMENSION(mp) :: sublayer_dz, eddy_shape,eddy_mod,soil_moisture_mod, &
                          soil_moisture_mod_sat, wb_liq, &
                          pore_size,pore_radius, rel_s,hk_zero,hk_zero_sat,time_scale  !note pore_size in m

   INTEGER, DIMENSION(mp) :: int_eddy_shape

   REAL(r_2), parameter :: Dff=2.5e-5, &
                      lm=1.73e-5, &
                      pi = 3.14159265358979324, &
                      c2 = 2.0

   real(r_2), parameter :: rtevap_max = 10000.0

   integer :: i,j,k 
   logical, save ::  first_call = .true.

   if (first_call) canopy%sublayer_dz(:) = 0.001

   pore_radius(:) = 0.148  / (1000.0*9.81*abs(soil%sucs_vec(:,1))/1000.0)  !should replace 0.148 with surface tension, unit coversion, and angle
   pore_size(:) = pore_radius(:)*sqrt(pi)

      !scale ustar according to the exponential wind profile, assuming we are a mm from the surface
      eddy_shape = 0.3*met%ua/ max(1.0e-4,canopy%us*exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff))))
      int_eddy_shape = floor(eddy_shape)
      eddy_mod(:) = 0.0
      do i=1,mp   
         eddy_mod(i) = 2.2*sqrt(112.0*pi) / (2.0**(eddy_shape(i)+1.0) * sqrt(eddy_shape(i)+1.0))

         if (int_eddy_shape(i) .gt. 0) then
            eddy_mod(i) = eddy_mod(i) / my_gamma(eddy_shape(i)+1.0) * (2.0*eddy_shape(i)+1.0)
            do k=1,int_eddy_shape(i)
               eddy_mod(i) = eddy_mod(i) * (2.0*(eddy_shape(i) - k) + 1.0)
            end do
         end if
      end do
      canopy%sublayer_dz = max(eddy_mod(:) * air%visc / max(1.0e-4,canopy%us*exp(-rough%coexp*(1.0-canopy%sublayer_dz/max(1e-2,rough%hruff)))),1e-7)              !exp(-canopy%vlaiw)), 1e-7)



   if (first_call) then
      wb_liq(:) = real(max(0.0001,min(pi/4.0, ssnow%wb(:,1)) ) )
   else
      wb_liq(:) = real(max(0.0001,min(pi/4.0, (ssnow%wb(:,1)-ssnow%wbice(:,1) - ssnow%satfrac(:)*soil%ssat_vec(:,1))/(1._r_2 - ssnow%satfrac(:)) ) ) )
   end if

   rel_s = real( max(wb_liq(:)-soil%watr(:,1),0._r_2)/(soil%ssat_vec(:,1)-soil%watr(:,1)) )
   hk_zero = max(0.001*soil%hyds_vec(:,1)*(min(max(rel_s,0.001_r_2),1._r_2)**(2._r_2*soil%bch_vec(:,1)+3._r_2) ),1e-8)
   hk_zero_sat = max(0.001*soil%hyds_vec(:,1),1e-8)

   soil_moisture_mod(:)     = 1.0/pi/sqrt(wb_liq)* ( sqrt(pi/(4.0*wb_liq))-1.0)
   soil_moisture_mod_sat(:) = 1.0/pi/sqrt(soil%ssat_vec(:,1))* ( sqrt(pi/(4.0*soil%ssat_vec(:,1)))-1.0)

   where(snow_covered .ne. 0)
      soil_moisture_mod = 0.
      soil_moisture_mod_sat = 0.
   elsewhere
      canopy%sublayer_dz = canopy%sublayer_dz + dz_litter
   endwhere


   where(canopy%sublayer_dz .ge. 1.0e-7) 
      ssnow%rtevap_unsat(:) = min( rough%z0soil/canopy%sublayer_dz * (lm/ (4.0*hk_zero) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod) / Dff),&  
                         rtevap_max )
      ssnow%rtevap_sat(:)  = min( rough%z0soil/canopy%sublayer_dz * (lm/ (4.0*hk_zero_sat) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod_sat) / Dff),& 
                         rtevap_max )

   elsewhere
      ssnow%rtevap_unsat(:) = min( lm/ (4.0*hk_zero) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod) / Dff,rtevap_max)
      ssnow%rtevap_sat(:)  = min( lm/ (4.0*hk_zero_sat) + (canopy%sublayer_dz + pore_size(:) * soil_moisture_mod_sat) / Dff,rtevap_max)
   endwhere
   !no additional evap resistane over lakes
   where(veg%iveg .eq. 16) 
      ssnow%rtevap_sat = 0.0
      ssnow%rtevap_unsat = 0.0
   endwhere

   first_call = .false.


END SUBROUTINE or_soil_evap_resistance

END MODULE cable_psm

