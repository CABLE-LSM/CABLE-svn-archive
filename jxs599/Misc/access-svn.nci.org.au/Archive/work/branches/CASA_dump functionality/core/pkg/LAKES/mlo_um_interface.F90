module mlo_um_interface

! This is the interface between CABLE and the mlo code (for lakes) in
! the ACCESS model.

! To do:
!  - U and V wind components fed to model
!  - Include coriols term
!  - Store free surface height in D1 array

implicit none

public  intialize_mlo_vars

integer, parameter :: veg_type_lake = 16  ! 17 tile case

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This subroutine initalises the lake routines.  This must be
   ! called before the mloeval.

   subroutine intialize_mlo_vars(fktau         & 
         & ,land_pts, ntiles, l_tile_pts       &
#include "arg_fields_mlo.h"         
         & )
   
   use mlo  
   use define_dimensions, only : mp
   use cable_variables, only : veg,ssoil, rad
!   use cable_common_module, only : ktau_gl
   
   implicit none

   integer, intent(in) :: fktau 
   integer :: land_pts, ntiles
   logical, dimension(land_pts, ntiles) :: l_tile_pts
#include "typ_fields_mlo2.h"

   integer k
   real, dimension(mp) :: dep,ochgt
   real, dimension(mp,wlev,4) :: mlodwn
   real, dimension(mp,10) :: icedwn
   

   if (mp.eq.0) return

	 
   if (fktau == 0) then
     ! Need lake bathymetry.  Currently assume all lakes are 32.2m
     ! deep and wlev=6 for compatibility with CABLE.  Great lakes
     ! are approx 200m deep which requires wlev=12 in mlo.f90.
     
     dep=0.
     where (veg%iveg.eq. veg_type_lake)
       ! keep 32.2m lakes for now ---------
       !where (rad%longitude < 0.)
          dep=32.2
       !elsewhere
       !   dep=60.
       !end where
     end where
     !if(maxval(dep) > 0.) write(6,*) "iveg = 7 (inland water), lake depth=", dep 
     
     write(6, *) "mp     = ", mp
     write(6, *) "dep xn = ",maxval(dep),minval(dep)
     
     call mloinit(mp,dep,0)

     !check water temperature at first level to decide initialisation
     if(MAXVAL(MLO_WATER_TEMP(:,:,1)) > 200.0) then 
!    if (ktau_gl.gt.1) then ! import data from dump
       write(6,*) "mlo init from dump file"

       ! Assume lake depth is unchanged in restart file since
       ! UM is not relocatable with a warm start
       !call mloregrid(ocndwn,mlodwn)
       do k=1,wlev
         mlodwn(:,k,1)=pack(MLO_WATER_TEMP(:,:,k), l_tile_pts)     ! water temperature (K)
         mlodwn(:,k,2)=pack(MLO_WATER_SALINITY(:,:,k), l_tile_pts) ! water salinity (PSU)
         mlodwn(:,k,3)=pack(MLO_WATER_UCUR(:,:,k), l_tile_pts)     ! water U (m/s)
         mlodwn(:,k,4)=pack(MLO_WATER_VCUR(:,:,k), l_tile_pts)     ! water V (m/s)
       end do
       icedwn(:,1)=pack(MLO_SNOW_SURF_TEMP, l_tile_pts)            ! ice temperature surf (K)
       icedwn(:,2)=pack(MLO_SNOW_LEV1_TEMP, l_tile_pts)            ! ice temperature t0 (K)
       icedwn(:,3)=pack(MLO_ICE_LEV2_TEMP, l_tile_pts)             ! ice temperature t1 (K)
       icedwn(:,4)=pack(MLO_ICE_LEV3_TEMP, l_tile_pts)             ! ice temperature t2 (K)
       icedwn(:,5)=pack(MLO_ICE_FRAC, l_tile_pts)                  ! ice fraction (0-1)
       icedwn(:,6)=pack(MLO_ICE_DEPTH, l_tile_pts)                 ! ice depth (m)
       icedwn(:,7)=pack(MLO_SNOW_DEPTH, l_tile_pts)                ! snow depth (m)
       icedwn(:,8)=pack(MLO_ICE_ENERGY, l_tile_pts)                ! Internal energy store (Ws/m^2)
       icedwn(:,9)=0.                                              ! Ice U (m/s)
       icedwn(:,10)=0.                                             ! Ice V (m/s)
       ochgt=0.                                                    ! Free surface height (m)
     else ! new initial conditions
       write(6,*) "mlo init default due to zero values in dump file"

       do k=1,wlev
         mlodwn(:,k,1)=max(ssoil%tgg(:,1),272.)
         mlodwn(:,k,2)=0. ! salinity=0 PSU
         mlodwn(:,k,3)=0.
         mlodwn(:,k,4)=0.
       end do
       icedwn(:,1:4)=273.16
       icedwn(:,5:10)=0.
       ochgt=0.
     end if
     call mloload(mlodwn,ochgt,icedwn,0)
   end if
	 
    return
    end subroutine intialize_mlo_vars


    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine updates the lake prognostic variables and
    ! returns fluxes back to CABLE

    subroutine mlorun(TIMESTEP                 &
         & ,land_pts, ntiles, l_tile_pts       &
	 & ,calprog                            &
#include "arg_fields_mlo.h"
         & )

    use cable_variables
    use define_dimensions, only : mp
    use mlo
    
    implicit none

    real, intent(in) :: TIMESTEP
   integer :: land_pts, ntiles
   logical, dimension(land_pts, ntiles) :: l_tile_pts
   logical, intent(in) :: calprog
#include "typ_fields_mlo2.h"


    real, dimension(mp) :: f,vwnd,ua_10m,dep,vnratio,epan 
    real, dimension(mp) :: sicefrac,sicedepth,snowdepth,inflow
    real, dimension(mp,wlev,4) :: mlodwn
    real, dimension(mp,10) :: icedwn
    real, dimension(mp) :: ochgt
    integer :: k  !MJT lake    
    real :: miss   
 
    inflow = 0.
    miss = 0.0

    if (mp.eq.0) return
    
    !if(UBOUND(mlodwn, 1) /= mp) then 
    !  write(6, *) UBOUND(mlodwn, 1), ", mp is different from that of initialisation: ", mp
    !end if
    
    if(any(veg%iveg.eq. veg_type_lake)) then     
       call mloalb4(1,mp,met%coszen,rad%reffbm(:,1),rad%reffdf(:,1), &
                    rad%reffbm(:,2),rad%reffdf(:,2),0)
       f=0.    ! actually a dummy for f (currently not used)
       vwnd=0. ! actually a dummy for v wind component (need to fix this)
       epan=0.
       sicefrac=0.
       sicedepth=0.
       snowdepth=0.
       vnratio=0.
       where (met%fsd(:,3).ne.0.)
         vnratio=met%fsd(:,1)/met%fsd(:,3)
       endwhere
       call mloeval(rad%trad,rough%z0m,canopy%cduv,canopy%fh,               &  !OUTPUT
                    canopy%fe,ssoil%wetfac,ssoil%potev,epan,                &  !OUTPUT
                    sicefrac,sicedepth,snowdepth,TIMESTEP,                  &  !INPUT
                    rough%za_uv,rough%za_tq,met%fsd(:,3),met%fld,           &  !INPUT
                    met%precip/TIMESTEP,met%ua,vwnd,met%tk,met%qv,          &  !INPUT
                    met%pmb*100.,f,vnratio,                                 &  !INPUT
                    rad%fbeam(:,1),rad%fbeam(:,2),inflow,0,calprog)            !INPUT
       call mloscrnout(canopy%tscrn,canopy%qscrn,canopy%uscrn,ua_10m,0)
       call mlosave(mlodwn,dep,ochgt,icedwn,0)

       do k=1,wlev
MLO_WATER_TEMP(:,:,k)     = unpack(mlodwn(:,k,1), l_tile_pts, miss) ! Water temperature (K)
MLO_WATER_SALINITY(:,:,k) = unpack(mlodwn(:,k,2), l_tile_pts, miss) ! Water salinity (PSU)
MLO_WATER_UCUR(:,:,k)     = unpack(mlodwn(:,k,3), l_tile_pts, miss) ! Water U current (m/s)
MLO_WATER_VCUR(:,:,k)     = unpack(mlodwn(:,k,4), l_tile_pts, miss) ! Water V current (m/s)
       end do
MLO_SNOW_SURF_TEMP = unpack(icedwn(:,1), l_tile_pts, miss)  ! Ice/snow surface surface (K)
MLO_SNOW_LEV1_TEMP = unpack(icedwn(:,2), l_tile_pts, miss)  ! Snow Level 1 (if exists) temperature
MLO_ICE_LEV2_TEMP  = unpack(icedwn(:,3), l_tile_pts, miss)  ! Ice Level 2 (if exists) temperature
MLO_ICE_LEV3_TEMP = unpack(icedwn(:,4), l_tile_pts, miss)   ! Ice level 3 (if exists) temperature
MLO_ICE_FRAC      = unpack(icedwn(:,5), l_tile_pts, miss)   ! Ice fraction (0-1) 
MLO_ICE_DEPTH     = unpack(icedwn(:,6), l_tile_pts, miss)   ! Ice depth (m) 
MLO_SNOW_DEPTH    = unpack(icedwn(:,7), l_tile_pts, miss)   ! Snow (on ice) depth (m)        
MLO_ICE_ENERGY    = unpack(icedwn(:,8), l_tile_pts, miss)   ! Stored energy in ice (brine pockets) (Ws/m2)
!MLO_ICE_???       = unpack(icedwn(:,9), l_tile_pts, miss)  ! Ice U current (m/s)
!MLO_ICE_???       = unpack(icedwn(:,10), l_tile_pts, miss) ! Ice V current (m/s)
!MLO_WATER_???     = unpack(ochgt, l_tile_pts, miss)        ! Free surface height (m)

       where (veg%iveg.eq. veg_type_lake)
        rad%albedo(:,1)=rad%fbeam(:,1)*rad%reffbm(:,1)+rad%fbeam(:,2)*rad%reffdf(:,1) ! approx
        rad%albedo(:,2)=rad%fbeam(:,1)*rad%reffbm(:,2)+rad%fbeam(:,2)*rad%reffdf(:,2) ! approx
       end where
    end if
    
    return
    end subroutine mlorun
    
    
end module mlo_um_interface
