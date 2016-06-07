module cable_TwoDim_GW


   !add ssnow%GWconvergence(i) and soil%dx in soil_define_types
   !also need ssnow%elv

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp,mlat,mlon
  USE cable_data_module, ONLY : issnow_type, point2constants
  USE cable_IO_vars_module, ONLY: landpt, patch, max_vegpatches, parID_type,  &
                          metGrid, land_x, land_y, logn, output,               &
                          xdimsize, ydimsize, check

  USE cable_2dgw_types

  IMPLICIT NONE

  PRIVATE
 
  TYPE ( issnow_type ) :: C 

  PUBLIC lateral_fluxes

  contains


  subroutine lateral_fluxes(dels,ssnow,soil,latitude,map_indices,northern_halo_parm,southern_halo_parm, &
                                                        northern_halo_var,southern_halo_var)
   
    implicit none

    real, intent(in)                         :: dels
    type(soil_snow_type),      intent(inout) :: ssnow
    type(soil_parameter_type), intent(in)    :: soil
    real, dimension(:),        intent(in)    :: latitude
    integer, dimension(:,:),   intent(in)    :: map_indices

    type(gw_halo_param_type), optional,intent(in)  :: northern_halo_parm
    type(gw_halo_param_type), optional,intent(in)  :: southern_halo_parm
    type(gw_halo_var_type)  , optional,intent(in)  :: northern_halo_var
    type(gw_halo_var_type)  , optional,intent(in)  :: southern_halo_var


    real(r_2), dimension(mlon,0:mlat+1) ::  &
        head              ! head = elev - wtd

    integer ::  istep,i,j,k,kk,klev,ii,jj,ib,ie,jb,je
! Local arrays:
    !all need to be allocatable for mpi, mlon,mlat are the worker values
    real(r_2),  dimension(0:mlon+1,0:mlat+1)   :: Qlateral    ! 
    real(r_2),  dimension(0:mlon+1,0:mlat+1)   :: hksat_lateral    ! 
    real(r_2),  dimension(0:mlon+1,0:mlat+1) :: trans_lateral,fdepth
    real(r_2)  :: soil_z_depth,delta_mass,Qtemp
    real(r_2), dimension(mlon,0:mlat+1) :: hkfactor_lateral,area
    real(r_2) :: radius_earth,delta_radians,dXY,factor

    integer, dimension(0:mlon+1,0:mlat+1)  :: map_indices_wrap

    integer :: j_b,j_e

    !if we are not using mpi map_indices only goes from 1:mlon,1:mlat
    !need it to go from 0:n+1, so create it this way
    !otherwise just copy it
    if ((size(map_indices,dim=1) .eq. mlon) .and. (size(map_indices,dim=2).eq.mlat)) then
       map_indices_wrap(1:mlon,1:mlat) = map_indices
       map_indices_wrap(0,1:mlat) = map_indices(mlon,1:mlat)
       map_indices_wrap(mlon+1,1:mlat) = map_indices(1,1:mlat)
       map_indices_wrap(:,0) = -1
       map_indices_wrap(:,mlat+1) = -1
    else
       map_indices_wrap(:,:) = map_indices
    end if
    

    radius_earth = 6.371e6
    hkfactor_lateral(:,:) = 30.0
    delta_radians = 0.25*3.14159/180.0
    dXY = delta_radians*radius_earth

   soil_z_depth = real(sum(soil%zse,dim=1),r_2)

   Qlateral(:,:) = 0._r_2

   !need to pass soil%hksat,hk_factor_lateral,slope,wtd .... all one dim~!!!! yeah!!!!!!!
   !use a temp variable to recieve from other workers, and normal to send

   do j=1,mlat
      do i=1,mlon
         !for mpi version this will have to be the worker map_index, which is global_mask_index - points for lower workers
         !map indices is passed in, not the global map index unless not using mpi
         k = map_indices_wrap(i,j)

         if (k .gt. 0) then
            !hard coded for gswp3 forcing
            area(i,j) = dXY*dXY*cos(latitude(k)*3.14159/180.0)

            fdepth(i,j) = 100._r_2 / (1._r_2 + 150._r_2 * soil%slope(k))

            if (ssnow%wtd(k) .gt. soil_z_depth*1000._r_2) then
               trans_lateral(i,j) = 0.001*soil%hksat(k,ms)*hkfactor_lateral(i,j)*fdepth(i,j)*exp((ssnow%wtd(k)/1000._r_2-soil_z_depth)/fdepth(i,j))
            else
               trans_lateral(i,j) = 0.001*soil%hksat(k,ms)*(hkfactor_lateral(i,j)*fdepth(i,j) + (soil_z_depth - ssnow%wtd(k)/1000._r_2))
            end if
            !no movement when frozen
            if (ssnow%tgg(k,ms) .lt. 270.) trans_lateral(i,j) = 0.

            head(i,j) = soil%elev(k) - ssnow%wtd(k)/1000._r_2

         else

            head(i,j) = 0.

         end if
      end do

  end do


  !fill edges with halo data
  if (present(northern_halo_parm) .and. present(northern_halo_var)) then
     j = 0
     do i=1,mlon
         !for mpi version this will have to be the worker map_index, which is global_mask_index - points for lower workers
         !map indices is passed in, not the global map index unless not using mpi
         k = map_indices_wrap(i,j)

         if (k .gt. 0) then
            !hard coded for gswp3 forcing
            area(i,j) = dXY*dXY*cos(northern_halo_parm%latitude(k)*3.14159/180.0)

            fdepth(i,j) = 100._r_2 / (1._r_2 + 150._r_2 * soil%slope(k))

            if (northern_halo_var%wtd(k) .gt. soil_z_depth*1000._r_2) then
               trans_lateral(i,j) = 0.001*northern_halo_parm%hksat(k)*hkfactor_lateral(i,j)*fdepth(i,j)*exp((northern_halo_var%wtd(k)/1000._r_2-soil_z_depth)/fdepth(i,j))
            else
               trans_lateral(i,j) = 0.001*northern_halo_parm%hksat(k)*(hkfactor_lateral(i,j)*fdepth(i,j) + (soil_z_depth - northern_halo_var%wtd(k)/1000._r_2))
            end if
            !no movement when frozen
            if (northern_halo_var%tgg_ms(k) .lt. 270.) trans_lateral(i,j) = 0.

            head(i,j) = northern_halo_parm%elev(k) - northern_halo_var%wtd(k)/1000._r_2

         else

            head(i,j) = 0.

         end if
     end do

  end if

  if (present(southern_halo_var) .and. present(southern_halo_parm)) then

     j=mlat+1
     do i=1,mlon
         !for mpi version this will have to be the worker map_index, which is global_mask_index - points for lower workers
         !map indices is passed in, not the global map index unless not using mpi
         k = map_indices_wrap(i,j)

         if (k .gt. 0) then
            !hard coded for gswp3 forcing
            area(i,j) = dXY*dXY*cos(southern_halo_parm%latitude(k)*3.14159/180.0)

            fdepth(i,j) = 100._r_2 / (1._r_2 + 150._r_2 * soil%slope(k))

            if (southern_halo_var%wtd(k) .gt. soil_z_depth*1000._r_2) then
               trans_lateral(i,j) = 0.001*southern_halo_parm%hksat(k)*hkfactor_lateral(i,j)*fdepth(i,j)*exp((southern_halo_var%wtd(k)/1000._r_2-soil_z_depth)/fdepth(i,j))
            else
               trans_lateral(i,j) = 0.001*southern_halo_parm%hksat(k)*(hkfactor_lateral(i,j)*fdepth(i,j) + (soil_z_depth - southern_halo_var%wtd(k)/1000._r_2))
            end if
            !no movement when frozen
            if (northern_halo_var%tgg_ms(k) .lt. 270.) trans_lateral(i,j) = 0.

            head(i,j) = southern_halo_parm%elev(k) - southern_halo_var%wtd(k)/1000._r_2

         else

            head(i,j) = 0.

         end if
     end do

   end if

   !fill eastern edge and wester edge with ghots values
   head(0,:)    = head(mlon,:)
   head(mlon+1,:) = head(1,:)
   trans_lateral(0,:)      = trans_lateral(mlon,:)
   trans_lateral(mlon+1,:) = trans_lateral(1,:)

   do j = 1,mlat

      jb = j - 1
      je = j + 1

      do i=1,mlon

         ib = i - 1
         ie = i + 1

         Qtemp = 0.

         if (map_indices_wrap(i,j) .ge. 1) then  !same as checking for land mask

            do jj=jb,je
               do ii=ib,ie
 
                  
                  if ((abs(jj-j) + abs(ii-i)) .eq. 2) then
                     factor = dXY/sqrt(2.0)
                  else
                     factor = dXY
                  end if
                  if (jj .ne. j) then
                     factor = factor*cos(latitude(k)*2.0*3.14159/360.0)
                  end if

                  if (map_indices_wrap(ii,jj) .ge. 1) then
                     Qtemp = Qtemp + 0.5*factor*(head(ii,jj)-head(i,j))*(trans_lateral(ii,jj) + trans_lateral(i,j))
                  end if
               end do
            end do

            Qlateral(i,j) = Qtemp  / area(i,j)

         end if
      end do ! i longitude
   end do  ! j latitude

   !Qlateral is in m/s
   head(1:mlon,1:mlat) = head(1:mlon,1:mlat) + Qlateral(1:mlon,1:mlat)*dels

   !translate into change in aquifer or soil water
   do j=1,mlat

      do i=1,mlon

         k = map_indices_wrap(i,j)

         if (k .gt. 0) then

            ssnow%wtd(k) = 1000._r_2*(soil%elev(k)-head(i,j))

            delta_mass = Qlateral(i,j)*dels*1000._r_2

            if (delta_mass .ge. 0.) then
               if (delta_mass .lt. soil%GWdz(k)*1000._r_2*(ssnow%GWwb(k)-soil%GWwatsat(k))) then
                  ssnow%GWwb(k) = ssnow%GWwb(k) + delta_mass/soil%GWdz(k)*1000._r_2
                  delta_mass = 0.
               else
                  delta_mass = delta_mass - (soil%GWwatsat(k)-ssnow%GWwb(k))*soil%GWdz(k)*1000._r_2
                  ssnow%GWwb(k) = soil%GWwatsat(k)
               end if

               if (delta_mass .gt. 0._r_2) then

                  do kk=ms,1,-1
                     if (delta_mass .lt. soil%zse(kk)*1000._r_2*(ssnow%wbliq(k,kk)-soil%watsat(k,kk))) then
                        ssnow%wbliq(k,kk) = ssnow%wbliq(k,kk) + delta_mass/soil%zse(kk)*1000._r_2
                        delta_mass = 0.
                     else
                        delta_mass = delta_mass - (soil%watsat(k,kk)-ssnow%wbliq(k,kk))*soil%zse(kk)*1000._r_2
                        ssnow%wbliq(k,kk) = soil%watsat(k,kk)
                     end if
                   end do

               end if

            else
  
               if (ssnow%wtd(k) .gt. 1000._r_2*sum(soil%zse,dim=1)) then

                  ssnow%GWwb(k) = ssnow%GWwb(k) + delta_mass/soil%GWdz(k)*1000._r_2
                  if (ssnow%GWwb(k) .lt. 0._r_2) then
                     delta_mass = (ssnow%GWwb(k)-1e-8)*soil%GWdz(k)*1000._r_2
                     ssnow%GWwb(k) = 1.e-8
                  else
                     delta_mass = 0.
                  end if

                  if (delta_mass .lt. 0.) then
                     do kk=ms,1,-1
                        if (delta_mass .lt. -soil%zse(kk)*1000._r_2*(ssnow%wbliq(k,kk))) then
                           ssnow%wbliq(k,kk) = ssnow%wbliq(k,kk) + delta_mass/soil%zse(kk)*1000._r_2
                           delta_mass = 0.
                        else
                           delta_mass = delta_mass + (ssnow%wbliq(k,kk)-1e-7)*soil%zse(kk)*1000._r_2
                           ssnow%wbliq(k,kk) = 1.e-7
                        end if
                     end do
                   end if
 
               else
                  klev = ms
                  do kk=2,ms-1
                     if (ssnow%wtd(kk) .ge. sum(soil%zse(1:kk-1),dim=1) .and. &
                         ssnow%wtd(kk) .le. sum(soil%zse(1:kk-1),dim=1) ) then
                         klev = kk
                     end if
                  end do

                  do kk=ms,klev,-1
                     if (delta_mass .lt. -soil%zse(kk)*1000._r_2*(ssnow%wbliq(k,kk))) then
                        ssnow%wbliq(k,kk) = ssnow%wbliq(k,kk) + delta_mass/soil%zse(kk)*1000._r_2
                        delta_mass = 0.
                     else
                        delta_mass = delta_mass + (ssnow%wbliq(k,kk)-1e-7)*soil%zse(kk)*1000._r_2
                        ssnow%wbliq(k,kk) = 1.e-7
                     end if
                  end do

               end if

            end if

         end if

      end do

   end do


  end subroutine lateral_fluxes
   
      
end module cable_TwoDim_GW
