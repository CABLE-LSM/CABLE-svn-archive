module cable_subgrid_transfer
   
   USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                                   r_2, ms, mp, mland
   USE cable_data_module, ONLY : issnow_type, point2constants

   USE cable_common_module, ONLY : gw_params

   IMPLICIT NONE

   PRIVATE

   TYPE ( issnow_type ) :: C 
   
contains

SUBROUTINE subgrid_sm_transfer(dels,ktau,ssnow,soil)
  USE cable_IO_vars_module , ONLY : landpt
  
  implicit none
  real, intent(in)                              :: dels
  integer, intent(in)                           :: ktau
  type(soil_snow_type)     , intent(inout)      :: ssnow
  type(soil_parameter_type), intent(in)         :: soil
  
  !local variables
  integer :: i,ii,k,kk,ib,ie
  real(r_2), dimension(mp)     :: qsrf_store_n,qsrf_flow_n, q_tot, wb_avg
  real(r_2), dimension(mp,ms+1) :: q_lev
  
  real(r_2) :: def, theta, ani,f, slope  !ani = anisotopic factor (2e3 for chen kumar 2001)

  logical :: verbose_debug

  verbose_debug = .false.
  ani = 2.0e3
  f = 1. / sum(soil%zse)
  
  do i=1,mland
  
    ib = landpt(i)%cstart
    ie = landpt(i)%cend
    do ii=ib+1,ie
    
      !slope = (soil%elevation(ii) - soil%elevation(ii-1)) / (soil%distance(ii)-soil%distance(ii-1))
      def   = sum((soil%watsat(ii,:)-(ssnow%wbliq(ii,:)+dri*ssnow%wbice(ii,:))))!*soil%zse(:))
      def   = def + (soil%GWwatsat(ii) - ssnow%GWwb(ii))!*soil%GWdz(ii)
      q_tot(ii) = ani*soil%hksat(ii,ms)*sin(soil%slope(ii))/soil%distance(ii)*def*soil%area(ii)
      ! sin(slope)*(exp(-f*sum(soil%zse)*cos(slope)) - exp(-f*def*cos(slope)))

    end do
    
    ii = ib
    !slope = (soil%elevation(ii+1) - soil%elevation(ii)) / (soil%distance(ii+1)-soil%distance(ii))
    def   = sum(max(soil%watsat(ii,:)-(ssnow%wbliq(ii,:)+dri*ssnow%wbice(ii,:)),0.))!*soil%zse(:))
    def   = def + (soil%GWwatsat(ii) - ssnow%GWwb(ii))!*soil%GWdz(ii)
    q_tot(ii) = ani*soil%hksat(ii,ms)*sin(soil%slope(ii))/soil%distance(ii)*def*soil%area(ii)
    !sin(slope)* (exp(-f*sum(soil%zse)*cos(slope)) - exp(-f*def*cos(slope)))
 
      if (ii .eq. 388 .and. verbose_debug) then
         write(*,*) soil%hksat(ii,ms)
         write(*,*) slope
         write(*,*) sin(soil%slope(ii))
         write(*,*) sum(soil%zse)
         write(*,*) cos(slope)
         write(*,*) exp(-f*sum(soil%zse)*cos(slope))
         write(*,*) 'def is ',def
         write(*,*) exp(-f*def*cos(slope))
         write(*,*) soil%area(ii)
       end if

  end do
  
  wb_avg(:) = 0.
  do i=1,mland
    ib = landpt(i)%cstart
    ie = landpt(i)%cend
    do ii=ib,ie
      wb_avg(ii) =max( 1e-6, &
                     sum( max(ssnow%wbliq(ii,:) - 0.2*soil%sfc(ii),0.)*soil%zse(:)*1000.0*ssnow%fracice(ii,:) ) + &
                        (ssnow%GWwb(ii) - 0.2*soil%sfc(ii))*soil%GWdz(ii)*1000.0*ssnow%fracice(ii,ms)    )
    end do
  end do
         
  q_lev(:,:) = 0.
  do i=1,mland
    ib = landpt(i)%cstart
    ie = landpt(i)%cend
    do ii=ib,ie
      q_lev(ii,1:ms) = q_tot(ii) * max(ssnow%wb(ii,:) - 0.2*soil%sfc(ii),0.)*soil%zse(:)*1000.0*ssnow%fracice(ii,:) / wb_avg(ii)
      q_lev(ii,ms+1) = q_tot(ii) * max(ssnow%GWwb(ii) - 0.2*soil%sfc(ii),0.)*soil%GWdz(ii)*1000.0*ssnow%fracice(ii,ms) / wb_avg(ii)
    end do
  end do

    
  do k=1,ms
    do i=1,mland
      ib = landpt(i)%cstart
      ie = landpt(i)%cend
      do ii=ib,ie
        if (ii .lt. ie) then
          ssnow%wbliq(ii,k) = ssnow%wbliq(ii,k) + (q_lev(ii+1,k) - q_lev(ii,k))*dels / (soil%zse(k)*1000.0) / soil%area(ii)   !mm^3
        else
          ssnow%wbliq(ii,k) = ssnow%wbliq(ii,k) - q_lev(ii,k)*dels / (soil%zse(k)*1000.0) / soil%area(ii)
        end if
      end do  !ib,ie
    end do  !mland
  end do  !ms
  
  k = ms + 1
    do i=1,mland
      ib = landpt(i)%cstart
      ie = landpt(i)%cend
      do ii=ib,ie
        if (ii .lt. ie) then
          ssnow%GWwb(ii) = ssnow%GWwb(ii) + (q_lev(ii+1,k) - q_lev(ii,k))*dels / (soil%GWdz(ii)*1000.0) / soil%area(ii)   !mm^3
        else
          ssnow%GWwb(ii) = ssnow%GWwb(ii) - q_lev(ii,k)*dels / (soil%GWdz(ii)*1000.0) / soil%area(ii)
        end if
      end do  !ib,ie
    end do  !mland
  ! end GW layer

  ssnow%wb = ssnow%wbice + ssnow%wbliq
  
  !total subsurface flow from the lowest tile elevation goes to the river
  do i=1,mland
    ib = landpt(i)%cstart
    ie = landpt(i)%cend
    !do ii=ib,ie 
    ii = ib
       ssnow%rnof2(ii) = sum(q_lev(ii,:))/soil%area(ii)
    do ii=ib+1,ie
       ssnow%rnof2(ii) = 0.0
    end do
    !end do
  end do
  
  theta = 0.4
  qsrf_store_n(:) = 0.  !both _n in kg instead of kg/m2 like the non _n variables
  qsrf_flow_n(:) = 0.
  
  if (ktau .le. 1) then
     ssnow%qsrf_store(:) = 50*soil%area(:)
     ssnow%qsrf_flow(:) = 0.
  end if
  
  do i=1,mland
    ib = landpt(i)%cstart
    ie = landpt(i)%cend
    
    do ii=ie,ib,-1
    
      qsrf_store_n(ii) = ((1.-theta)*ssnow%qsrf_store(ii) + ssnow%qsrf_flow(ii) + ssnow%qsrf_gen(ii)*soil%area(ii))
      if (ii .gt. ib) then
        qsrf_flow_n(ii-1) = qsrf_flow_n(ii-1) + theta * ssnow%qsrf_store(ii)
      end if
      
    end do  !ib,ie
  end do  !mland
  
  do i=1,mland
    ib = landpt(i)%cstart
    
    ssnow%rnof1(ib) = ssnow%qsrf_store(ib)/dels*theta / soil%area(ib)  !grid cell surface runoff to river
    ssnow%rnof1(ib+1:ie) = 0.
  end do
  
  !do i=1,mp
  !  ssnow%qsrf_store(i) = qsrf_store_n(i) / soil%area(i)
  !  ssnow%qsrf_flow(i)  = qsrf_flow_n(i) / soil%area(i)
  !end do
  
  
  
END SUBROUTINE  subgrid_sm_transfer


end module cable_subgrid_transfer
