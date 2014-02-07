module cable_TwoDim_GW


   !add ssnow%GWconvergence(i) and soil%dx in soil_define_types
   !also need ssnow%elv

  USE cable_def_types_mod, ONLY : soil_snow_type, soil_parameter_type,        &
                             veg_parameter_type, canopy_type, met_type,        &
                             balances_type, r_2, ms, mp,mlat,mlon
  USE cable_data_module, ONLY : issnow_type, point2constants
  USE cable_IO_vars_module, ONLY: landpt, patch, max_vegpatches, parID_type,  &
                          metGrid, land_x, land_y, logn, output,               &
                          xdimsize, ydimsize, check, mask

  IMPLICIT NONE

  PRIVATE
 
  TYPE ( issnow_type ) :: C 

  PUBLIC gwstep, update_GW

  contains
   
   
  subroutine gwstep(dels,ssnow,soil)! dx,              &
            !ltype, elev, bot,        &
            !hycond, poros, compres,  &
            !ho, h, convgw,           &
            ! ebot, eocn,              &
            ! dt, istep)

! taken from WRF-HYDRO initially
  ! Steps ground-water hydrology (head) through one timestep.
  ! Modified from Prickett and Lonnquist (1971), basic one-layer aquifer 
  ! simulation program, with mods by Zhongbo Yu(1997).
  ! Solves S.dh/dt = d/dx(T.dh/dx) + d/dy(T.dh/dy) + "external sources"
  ! for a single layer, where h is head, S is storage coeff and T is 
  ! transmissivity. 3-D arrays in main program (hycond,poros,h,bot)
  ! are 2-D here, since only a single (uppermost) layer is solved.
  ! Uses an iterative time-implicit ADI method.

  ! use module_hms_constants


    real, intent(in)                      :: dels
    type(soil_snow_type), intent(inout)   :: ssnow
    type(soil_parameter_type), intent(in) :: soil
    !LOCAL variables to map ssnow to previous code
    !note assume continantal??
    real(r_2), dimension(mlon,mlat) ::  &
        elev,           &  ! elev/bathymetry of sfc rel to sl (m) (supp)
        bot,            &  ! elev. aquifer bottom rel to sl (m)   (supp)
        hycond,         &  ! hydraulic conductivity (m/s per m/m) (supp)
        poros,          &  ! porosity (m3/m3)                     (supp)
        compres,        &  ! compressibility (1/Pa)               (supp)
        ho                 ! head at start of timestep (m)        (supp)

    real(r_2), dimension(mlon,mlat) ::  &
        h,              &  ! head, after ghmcompute (m)           (ret)
        convgw             ! convergence due to gw flow (m/s)     (ret)

    real(r_2)  :: ebot, eocn
    integer ::  istep,i
!       eocn  = mean spurious sink for h_ocn = sealev fix (m/s)(ret)
!               This equals the total ground-water flow across 
!               land->ocean boundaries.
!       ebot  = mean spurious source for "bot" fix (m/s) (returned)
!       time  = elapsed time from start of run (sec)
!       dels = timestep length (sec)
!       istep = timestep counter

! Local arrays:
    real(r_2), dimension(mlon,mlat)   :: sf2    ! storage coefficient (m3 of h2o / bulk m3)
    real(r_2), dimension(mlon,mlat,2) ::   t    ! transmissivity (m2/s)..1 for N-S,..2 for E-W
    real(r_2), dimension(0:mlon+mlat) :: b,g    ! work arrays


    real(r_2), parameter    :: botinc = 0.01  ! re-wetting increment to fix h < bot
!   parameter (botinc = 0.  )  ! re-wetting increment to fix h < bot
                                 ! (m); else no flow into dry cells
    real(r_2), parameter    :: delskip = 0.005 ! av.|dhead| value for iter.skip out(m)
    integer, parameter :: itermax = 10    ! maximum number of iterations
    integer, parameter :: itermin = 3     ! minimum number of iterations
    real(r_2), parameter    :: sealev = -1.     ! sea-level elevation (m)
      
    logical  :: ContinueLoop

    integer ::                &
        iter,                   &
        j,                      &
        jp,                     &
        ip,                     &
        ii,                     &
        n,                      &
        jj,                     &
        ierr,                   &
        ier
        
    real(r_2) ::                &
        dy,                     &
        dx,                     &
        e,                      &
        su,                     &
        sc,                     &
        shp,                    &
        bb,                     &
        dd,                     &
        aa,                     &
        cc,                     &
        w,                      &
        ha,                     &
        delcur,                 &
        dtot,                   &
        dtoa,                   &
        darea,                  &
        tareal,                 &
        zz
        
        
    !map the 1D vector of ssnow to the 2D matrix that is longitude x latitude
    !hold ocean points to h=0
      
      !procedure
      !DO i = 1, mland ! over all land grid points
            ! Write to temporary variable (area weighted average across all
            ! patches):
            !otmp3xyt(land_x(i), land_y(i), 1) 
            
            !land_x(i),land_y(i) does the mapping.  make sure these are available.  use io_vars
    elev(:,:)   = 0._r_2
    poros(:,:)  = 0._r_2
    ho(:,:)     = 0._r_2
    hycond(:,:) = 0._r_2
    do i=1,mp
      elev(land_x(i),land_y(j)) = ssnow%elevation(i)
      poros(land_x(i),land_y(j)) = soil%GWwatsat(i)
      ho(land_x(i),land_y(j))    = ssnow%elevation(i) - ssnow%wtd(i)/1000._r_2
      hycond(land_x(i),land_y(j)) = ssnow%hk(i,ms)
    end do
    compres = 0.0_r_2
    h       = ho  !initial set it
    convgw  = 0._r_2
    bot     = elev - 100.0

    dx = ssnow%delx                           !assumed constant
    dy = ssnow%dely                           !need to add dx to global params
    darea = dx*dy
!       Top of iterative loop for ADI solution
    ContinueLoop = .true.
    iter = 0


    do while(ContinueLoop)

      iter = iter+1


      e    = 0.       ! absolute changes in head (for iteration control)
!       Set storage coefficient (sf2)
       tareal = 0.


!#OMP PARALLEL DO PRIVATE(j,i,su,sc,shp)
      do j=1,mlat
        do i=1,mlon
         tareal = tareal + darea
          su = poros(i,j)                    ! new (volug)
          sc = 1._r_2
 
          if      (ho(i,j).le.elev(i,j) .and. h(i,j).le.elev(i,j)) then
            sf2(i,j) = su
          else if (ho(i,j).ge.elev(i,j) .and. h(i,j).ge.elev(i,j)) then
            sf2(i,j) = sc
          else if (ho(i,j).le.elev(i,j) .and. h(i,j).ge.elev(i,j)) then
            shp = sf2(i,j) * (h(i,j) - ho(i,j))
            sf2(i,j) = shp * sc / (shp - (su-sc)*(elev(i,j)-ho(i,j)))
          else if (ho(i,j).ge.elev(i,j) .and. h(i,j).le.elev(i,j)) then
            shp = sf2(i,j) * (ho(i,j) - h(i,j))
            sf2(i,j) = shp * su / (shp + (su-sc)*(ho(i,j)-elev(i,j)))
          endif

        enddo
      enddo
!#OMP END PARALLEL DO
!==========================
!       Column calculations
!==========================

!       Set transmissivities. Use min(h,elev)-bot instead of h-bot,
!       since if h > elev, thickness of groundwater flow is just
!       elev-bot.

!#OMP PARALLEL DO PRIVATE(j,jp,i,ip)
      do j=1,mlat
        jp = min (j+1,mlat)
        do i=1,mlon
          ip = min (i+1,mlon)

          t(i,j,2) = sqrt( abs(                                           &
                        hycond(i, j)*(min(h(i ,j),elev(i ,j))-bot(i ,j))  &
                       *hycond(ip,j)*(min(h(ip,j),elev(ip,j))-bot(ip,j))  &
                         )    )                                           &
                   * (0.5*(dy+dy)) & ! in WRF the dx and dy are usually equal
                   / (0.5*(dx+dx))

          t(i,j,1) = sqrt( abs(                                           &
                        hycond(i,j )*(min(h(i,j ),elev(i,j ))-bot(i,j ))  &
                       *hycond(i,jp)*(min(h(i,jp),elev(i,jp))-bot(i,jp))  &
                         )    )                                           &
                   * (0.5*(dx+dx))  &
                   / (0.5*(dy+dy))
        enddo
      enddo
!#OMP END PARALLEL DO

      b = 0.
      g = 0.

      do ii=1,mlon
      
        i=ii
        if (mod(istep+iter,2).eq.1) i=mlon-i+1

!          calculate b and g arrays
        do j=1,mlat
       
          bb = (sf2(i,j)/dels) * darea
          dd = ( ho(i,j)*sf2(i,j)/dels) * darea
          aa = 0.0
          cc = 0.0
          
          if (j .gt. 1) then
            aa = -t(i,j-1,1)
            bb = bb + t(i,j-1,1)
          end if

          if (j .lt. mlat) then
            cc = -t(i,j,1)
            bb = bb + t(i,j,1)
          end if

          if (i .gt. 1) then
            bb = bb + t(i-1,j,2)
            dd = dd + h(i-1,j)*t(i-1,j,2)
          end if

          if (i .lt. mlon) then
            bb = bb + t(i,j,2)
            dd = dd + h(i+1,j)*t(i,j,2)
          end if

          w = bb - aa*b(j-1)
          b(j) = cc / w
          g(j) = (dd - aa*(g(j-1)))/w

        end do  !inner loop

!          re-estimate heads
        e = e + abs(h(i,mlat)-g(mlat))
        h(i,mlat) = g(mlat)
        n = mlat-1
        
        do while (n .ne. 0)
          ha = g(n) - b(n)*h(i,n+1)
          e  = e + abs(ha-h(i,n))
          h(i,n) = ha
          n = n - 1
        end do

      end do  !outer loop

!=======================
!       Row calculations
!=======================

      do j=1,mlat
        jp = min (j+1,mlat)
        do i=1,mlon
          ip = min (i+1,mlon)
          t(i,j,2) = sqrt( abs(                                         &
                   hycond(i, j)*(min(h(i ,j),elev(i ,j))-bot(i ,j))    &
                  *hycond(ip,j)*(min(h(ip,j),elev(ip,j))-bot(ip,j))    &
                    )    )                                             &
!                    * (0.5*(dy(i,j)+dy(ip,j)))                  &
!                    / (0.5*(dx(i,j)+dx(ip,j)))
                   * (0.5*(dy+dy))                               &
                   / (0.5*(dx+dx))

          t(i,j,1) = sqrt( abs(                                          &
                    hycond(i,j )*(min(h(i,j ),elev(i,j ))-bot(i,j ))    &
                   *hycond(i,jp)*(min(h(i,jp),elev(i,jp))-bot(i,jp))    &
                         )    )                                         &
                   * (0.5*(dx+dx))                               &
                   / (0.5*(dy+dy))
        enddo
      enddo

      b = 0.
      g = 0.

      do jj=1,mlat

        j=jj
        if (mod(istep+iter,2).eq.1) j = mlat-j+1
!          calculate b and g arrays
        do i=1,mlon
       
          bb = (sf2(i,j)/dels) * darea
          dd = ( ho(i,j)*sf2(i,j)/dels) * darea
          aa = 0.0
          cc = 0.0

          if (j .gt. 1) then
            bb = bb + t(i,j-1,1)
            dd = dd + h(i,j-1)*t(i,j-1,1)
          end if

          if (j .lt. mlat) then
            dd = dd + h(i,j+1)*t(i,j,1)
            bb = bb + t(i,j,1)
          end if

          if (i .gt. 1) then
            bb = bb + t(i-1,j,2)
            aa = -t(i-1,j,2)
          end if
          
          if (i .lt. mlon) then
            bb = bb + t(i,j,2)
            cc = -t(i,j,2)
          end if

          w = bb - aa*b(i-1)
          b(i) = cc/w
          g(i) = (dd-aa*g(i-1))/w
          
        end do   !inner loop

!          re-estimate heads
        e = e + abs(h(mlon,j)-g(mlon))
        h(mlon,j) = g(mlon)
        n = mlon-1
        
        do while (n .ne. 0)
          ha = g(n)-b(n)*h(n+1,j)
          e = e + abs(h(n,j)-ha)
          h(n,j) = ha
          n = n-1
        end do

      end do     !outer loop


      do j=1,mlat
        do i=1,mlon
          if (h(i,j).le.bot(i,j)+botinc) then
            e = e +  bot(i,j) + botinc - h(i,j)
            ebot = ebot + (bot(i,j)+botinc-h(i,j))*sf2(i,j)*darea
            h(i,j) = bot(i,j) + botinc
          endif
        enddo
      enddo
!        maintain head = sea level for ocean (only for adjacent ocean,
!        rest has hycond=0)
      do j=1,mlat
        do i=1,mlon
            eocn = eocn + (h(i,j)-sealev)*sf2(i,j)*darea
            h(i,j) = sealev
        enddo
      enddo

!        Loop back for next ADI iteration

      delcur = e/(mlon*mlat)


      if ( (delcur.gt.delskip*dels .and. iter.lt.itermax)      &
           .or. iter.lt.itermin ) then
         ContinueLoop = .true.
      else
         ContinueLoop = .false.  
      endif
          
    end do   !main iteration loop. leaves when iter >= itermax or delcur is small enough

!        Compute convergence rate due to ground water flow (returned)

    do j=1,mlat
      do i=1,mlon
        convgw(i,j) = sf2(i,j) * (h(i,j)-ho(i,j)) / dels
      enddo
    enddo

      !map the 2D convergence of water to the 1D cable array
    do i=1,mp
      ssnow%GWconvergence(i) = convgw(land_x(i),land_y(i)) / 1000._r_2 * darea * dels         ![mm]
      ssnow%wtd(i)           = (elev(land_x(i),land_y(i)) - h(land_x(i),land_y(i)))* 1000.0  ![mm]
    end do

!        Diagnostic water conservation check for this timestep
    dtot = 0.     ! total change in water storage (m3)
    dtoa = 0.

    do j=1,mlat
      do i=1,mlon
         dtot = dtot + sf2(i,j) *(h(i,j)-ho(i,j)) * darea
         dtoa = dtoa + sf2(i,j) * abs(h(i,j)-ho(i,j)) * darea
      enddo
    enddo

    dtot = (dtot/tareal)/dels   ! convert to m/s, rel to land area
    dtoa = (dtoa/tareal)/dels
    eocn = (eocn/tareal)/dels
    ebot = (ebot/tareal)/dels

      
    return
  end subroutine gwstep



!!!!!----------------------------------------------!!!!!
  subroutine update_GW(dels,ssnow,soil)

    real, intent(in)                         :: dels
    TYPE(soil_snow_type),      INTENT(INOUT) :: ssnow 
    TYPE(soil_parameter_type), INTENT(IN)    :: soil
    

    !local variables
    INTEGER :: i,j,k
    REAL(r_2), DIMENSION(mp) :: xx,vol_av,mss_av


    CALL point2constants( C ) 

   
    vol_av(:) = 0._r_2
    mss_av(:) = 0._r_2
    xx(:) = 0._r_2

    do i=1,mp
      xx(i) = ssnow%GWconvergence(i)  !convergence of water in mm
      !try to add to GWwb

      vol_av(i) = max(soil%GWwatsat(i) - ssnow%GWwb(i),0._r_2)
      mss_av(i) = vol_av(i)*soil%GWdz(i)*C%denliq

      if (xx(i) .le. mss_av(i)) then
        ssnow%GWwb(i) = ssnow%GWwb(i) + xx(i)/(soil%GWdz(i)*C%denliq)
        xx(i) = 0._r_2
      elseif (xx(i) .gt. mss_av(i)) then
        xx(i) = xx(i) - mss_av(i)
        ssnow%GWwb(i) = soil%GWwatsat(i)
      end if


      do k=ms,1,-1
        if (xx(i) .gt. 1e-7) then
          vol_av(i) = soil%watsat(i,k) - ssnow%wb(i,k)
          mss_av(i) = vol_av(i) * soil%zse(k) * C%denliq

          if (xx(i) .le. mss_av(i)) then
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xx(i) / (soil%zse(k)*C%denliq)
            xx(i) = 0._r_2
          elseif (xx(i) .gt. mss_av(i)) then
            xx(i) = xx(i) - mss_av(i)
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + mss_av(i) / (soil%zse(k)*C%denliq)
          end if
        end if
      end do   !soil layer loop

    end do  !mp patch loop



    where(xx(:) .gt. 1e-7)   !if column is fully saturated gets added to surface runoff
      ssnow%rnof1 = ssnow%rnof1 + xx(:) / dels
      ssnow%runoff = ssnow%rnof2 + ssnow%rnof1  !adjust the total runoff as well
    end where
      


  end subroutine update_GW

      
      
      
end module cable_TwoDim_GW
