subroutine c_cable(n_params,n_timesteps,param_values,qh_fluxes,qle_fluxes) bind(c)

   use iso_c_binding
   use cable_offline_driver_module, only: cable_offline_driver


   integer(c_int) :: n_params,n_timesteps
   real(c_double), dimension(n_params) :: param_values
   real(c_double), dimension(n_timesteps) :: qh_fluxes
   real(c_double), dimension(n_timesteps) :: qle_fluxes

   !call cable_offline_driver(param_values,qh_fluxes,qle_fluxes)
   call cable_offline_driver(real(param_values,kind=8),real(qh_fluxes,kind=8),real(qle_fluxes,kind=8))
end subroutine 
