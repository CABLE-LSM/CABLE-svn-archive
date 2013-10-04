      
         !___cable vars de-referenced by UM
         !----------------------------------------------------------------------------   
         ! snow depth (liquid water), factor for latent heat
         real, intent(in), dimension(mp) :: ssoil_snowd, ssoil_cls
         ! surface wind speed (m/s)
         real, intent(in), dimension(mp) :: met_ua 
         ! latent heat for water (j/kg), dry air density (kg m-3)
         real, intent(in), dimension(mp) :: air_rlam, air_rho 
         ! frac SW diffuse transmitted thru canopy, rad. temp. (soil and veg)
         real, intent(in), dimension(mp) :: rad_trad,rad_transd 
         ! total latent heat (W/m2), total sensible heat (W/m2)
         real, intent(in), dimension(mp) :: canopy_fe, canopy_fh  
         ! fraction of canopy wet
         real, intent(in), dimension(mp) :: canopy_fwet, canopy_wetfac_cs
         ! friction velocity, drag coefficient for momentum
         real, intent(in), dimension(mp) :: canopy_us, canopy_cdtq
         ! net rad. absorbed by surface (W/m2), total potential evaporation 
         real, intent(in), dimension(mp) :: canopy_rnet, canopy_epot        
         ! stability correction
         real, intent(in), dimension(mp,niter) :: canopy_zetar
         ! roughness length, Reference height for met forcing
         real, intent(in), dimension(mp) :: rough_z0m, rough_zref_tq 
         !----------------------------------------------------------------------------   
 
