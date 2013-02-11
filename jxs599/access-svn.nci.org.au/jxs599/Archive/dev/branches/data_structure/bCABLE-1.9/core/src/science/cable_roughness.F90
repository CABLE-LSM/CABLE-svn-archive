module roughness_module
  use define_dimensions
  implicit none
  private
  public ruff_resist

   contains


    ! m.r. raupach, 24-oct-92
    ! see: Raupach, 1992, BLM 60 375-395
    !      MRR notes "Simplified wind model for canopy", 23-oct-92
    !      MRR draft paper "Simplified expressions...", dec-92
    ! modified to include resistance calculations by Ray leuning 19 Jun 1998  
   subroutine ruff_resist(i, o, a, &          !roughness type in/out/auto
                           ssnow_snowd, &    !snow depth - from ssnow type 
                           ssnow_ssdnn, &    !snow depth - from ssnow type
                           veg_hc, & 
                           veg_vlai, & 
                           veg_iveg, & 
                           canopy_vlaiw, &
                           canopy_rghlai, &
                           cp ) 
    
      use cable_common_module, only : cable_runtime, cable_user
      use cable_diag_module,   only : cable_stat
      use cable_data_module ,  only : rough_in_type, rough_outtype, rough_auto_type, &
                                      physical_constants
      implicit none
      
      type (rough_in_type), intent(in) :: i 
      type (rough_out_type), intent(out) :: o
      type (rough_auto_type) :: a

      real, intent(in), dimension(mp) :: ssnow_snow, & ! snow depth (liquid water)
                                         ssnow_ssdnn   ! average snow density
 
      real, intent(in), dimension(mp) :: veg_hc, &  ! roughness height of canopy (veg - snow)
                                         veg_vlai   ! leaf area index

      integer, intent(in), dimension(mp) :: veg_iveg ! vegetation type
 
      real, intent(out), dimension(mp) :: canopy_vlaiw, &  ! lai adj for snow depth for calc of resistances
                                         canopy_rghlai  
      
      type (physical_constants) :: cp  
     
       
      !___ local variables 
      real, dimension(:), allocatable :: xx, &  ! =ccd*LAI; working variable 
                                         dh, &  ! d/h where d is zero-plane displacement
                                         hmax,& ! maximum height of canopy from
                                                ! tiles belonging to the same grid
                                         term1, term2, term3, term4, & ! for aerodyn resist. calc.
                                         z0soil ! roughness length of bare soil surface
     
      real, parameter :: 
             !jhan:Eva has changed 0.01 to 0.001
            min_canopy_ht = 0.01, &               ! Set canopy height above snow level:
     
      
      allocate( xx(mp), dh(mp), hmax(mp), term1(mp), term2(mp), &
                term3(mp), term4(mp), z0soil(mp) )
 
      if( cable_user%RUN_DIAG_LEVEL == 'BASIC' ) &    
         call cable_stat('ruff_resist')


      ! Set canopy height above snow level:
      a%rough_hruff = MAX( min_canopy_ht, veg_hc-1.2*ssnow_snowd/MAX(ssnow_ssdnn,100.) ) 
      ! maximum height of canopy from tiles belonging to the same grid
      hmax = i%rough_hruff_grmx
   
      
      ! LAI decreases due to snow and vegetation fraction:
      !jhan: check min_canopy_height here is desired
      canopy_vlaiw = veg_vlai * a%rough_hruff/MAX(min_canopy_ht, veg_hc )
      canopy_rghlai = canopy_vlaiw
  

      where(ssoil%snowd.lt.0.001.and.veg%iveg.ne.1) &
         canopy%rghlai = min( 3.,canopy%vlaiw)
    
    ! Roughness length of bare soil (m):
    z0soil = 1.e-6
    a%rough_z0soilsn = max(z0soil-0.5e-7*min(ssnow_snowd,20.),0.1e-7)


    WHERE (canopy_vlaiw.LT.0.01 .OR. a%rough_hruff.LT. a%rough_z0soilsn) ! BARE SOIL SURFACE
       o%rough_z0m = a%rough_z0soilsn
       a%rough_hruff = 0.0
       a%rough_rt0us = 0.0  
       a%rough_disp = 0.0
    ! Reference height zref is height above the displacement height
       !a%rough_zref_uv = max(3.5,i%rough_za_uv - a%rough_disp + hmax -  &
                                             min(1.,ssnow_snowd/max(ssnow_ssdnn,100.)))
       !o%rough_zref_tq = max(3.5,i%rough_za_tq - a%rough_disp + hmax - &
        !                                     min(1.,ssnow_snowd/max(ssnow_ssdnn,100.)))
!       rough_zref_uv = max(3.5,i%rough_za_uv + rough_hruff)
       a%rough_zref_uv = max(3.5,i%rough_za_uv )
       o%rough_zref_tq = max(3.5,i%rough_za_tq )

       a%rough_zruffs = 0.0
       a%rough_rt1usa = 0.0 
       a%rough_rt1usb = 0.0
       ! Friction velocity/windspeed at canopy height
       ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
       ! (usuhm set in physical_constants module):
       a%rough_usuh = MIN(SQRT(cp%csd+cp%crd*(canopy_vlaiw*0.5)), cp%usuhm)
       ! xx is ccd (see physical_constants) by LAI
       xx = SQRT(cp%ccd*MAX((canopy_vlaiw*0.5),0.0005))
       ! Displacement height/canopy height:
       ! eq.8 Raupach 1994, BLM, vol 71, p211-216
       dh = 1.0 - (1.0 - EXP(-xx))/xx
       ! Extinction coefficient for wind profile in canopy:
       ! eq. 3.14, SCAM manual (CSIRO tech report 132)
       a%rough_coexp = a%rough_usuh / (cp%vonk*cp%ccw_c*(1.0 - dh))
    ELSEWHERE ! VEGETATED SURFACE
       ! Friction velocity/windspeed at canopy height
       ! eq. 7 Raupach 1994, BLM, vol 71, p211-216
       ! (usuhm set in physical_constants module):
!jhan:Eva uses         
        !rough_usuh = MIN(SQRT(csd+crd*(canopy_vlaiw*0.5)), usuhm)
        a%rough_usuh = MIN(SQRT(cp%csd+cp%crd*(canopy_rghlai*0.5)), cp%usuhm)
       ! xx is ccd (see physical_constants) by LAI:
!jhan:Eva uses         
       !xx = SQRT(ccd*MAX((canopy_vlaiw*0.5),0.0005))
       xx = SQRT(cp%ccd*MAX((canopy_rghlai*0.5),0.0005))
       ! eq.8 Raupach 1994, BLM, vol 71, p211-216:
       dh = 1.0 - (1.0 - EXP(-xx))/xx
       ! Calculate zero-plane displacement:
       a%rough_disp = dh*a%rough_hruff
    ! Reference height zref is height above the displacement height
       !a%rough_zref_uv = max(3.5,i%rough_za_uv -  & 
       !                max(a%rough_disp,min(1.,ssnow_snowd/max(ssnow_ssdnn,100.))) + hmax)
       !o%rough_zref_tq = max(3.5,i%rough_za_tq -  &
       !                max(a%rough_disp,min(1.,ssnow_snowd/max(ssnow_ssdnn,100.))) + hmax)
       a%rough_zref_uv = max(3.5,i%rough_za_uv )
       o%rough_zref_tq = max(3.5,i%rough_za_tq )
!       rough_zref_uv = max(3.5,rough_za_uv + rough_hruff)
!       rough_zref_tq = max(3.5,rough_za_tq + rough_hruff)
       ! Calcualte roughness length:
       o%rough_z0m = ( (1.0 - dh) * EXP( LOG(cp%ccw_c) - 1. + 1./cp%ccw_c &
            & - cp%vonk/a%rough_usuh ) ) * a%rough_hruff
       ! find coexp: see notes "simplified wind model ..." eq 34a
       ! Extinction coefficient for wind profile in canopy:
       ! eq. 3.14, SCAM manual (CSIRO tech report 132)
       a%rough_coexp = a%rough_usuh / (cp%vonk*cp%ccw_c*(1.0 - dh))
!jhan:Eva uses         
       !jh term1  = EXP(2*csw*min(veg_vlaiw,3.)*(1-rough_disp/rough_hruff))
       !jh term2  = a33**2*ctl*2*csw*min(veg_vlaiw,3.)

       term1  = EXP(2*cp%csw*canopy_rghlai*(1-a%rough_disp/a%rough_hruff))
       term2  = cp%a33**2*cp%ctl*2*cp%csw*canopy_rghlai
       term3  = MAX((2./3.)*a%rough_hruff/a%rough_disp, 1.0)
       term4 =  exp(3.* a%rough_coexp*( a%rough_disp/ a%rough_hruff -1.))
       ! eq. 3.54, SCAM manual (CSIRO tech report 132)
       a%rough_rt0us  = term3*(cp%zdlin*LOG(cp%zdlin* a%rough_disp/ a%rough_z0soilsn) &
!jhan:Eva uses         
            !jh + (1-zdlin))*(EXP(2*csw*min(veg_vlaiw,3.)) - term1)/term2  ! &
            + (1-cp%zdlin))*(EXP(2*cp%csw*canopy_rghlai) - term1)/term2  ! &
!              / term4
       !        rt1 = turbulent resistance from canopy (z1 = disp) to
       !        reference level zref (from disp as origin). Normalisation:
       !        rt1us = us*rt1 = rt1usa + rt1usb + rt1usc
       !        with term a = resistance from disp to hruf
       !        term b = resistance from hruf to zruffs (or zref if zref<zruffs)
       !        term c = resistance from zruffs to zref (= 0 if zref<zruffs)
       !        where zruffs = SCALAR roughness sublayer depth (ground=origin)
       !        xzrufs = xdisp + xhruf*a33**2*ctl/vonk
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.49:
        a%rough_zruffs = a%rough_disp + a%rough_hruff*cp%a33**2*cp%ctl/cp%vonk/term3
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.51:
        a%rough_rt1usa = term3*(term1 - 1.0)/term2
        a%rough_rt1usb = term3*(MIN(o%rough_zref_tq+ a%rough_disp,  a%rough_zruffs) -  a%rough_hruff)/ &
            (cp%a33**2*cp%ctl*  a%rough_hruff)
        a%rough_rt1usb = MAX( a%rough_rt1usb,0.0)       ! in case zrufs < rough_hruff
    END WHERE


      deallocate( xx, dh, hmax, term1, term2, term3, term4, z0soil)

      return
  END SUBROUTINE ruff_resist


END MODULE roughness_module














