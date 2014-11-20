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

  PUBLIC gwstep, update_GW,mpi_step_gw_model

  contains
   
   
  subroutine gwstep(dels,ssnow,soil)! dx,              &
            !ltype, evl, bot,        &
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
    implicit none

    real, intent(in)                         :: dels
    type(soil_snow_type), intent(inout)      :: ssnow
    type(soil_parameter_type), intent(in)    :: soil
    !LOCAL variables to map ssnow to previous code
    !note assume continantal??
    real(r_2),  dimension(mlon,mlat) ::  &
        evl,           &  ! evl/bathymetry of sfc rel to sl (m) (supp)
        bot,            &  ! evl. aquifer bottom rel to sl (m)   (supp)
        hycond,         &  ! hydraulic conductivity (m/s per m/m) (supp)
        poros,          &  ! porosity (m3/m3)                     (supp)
        compres,        &  ! compressibility (1/Pa)               (supp)
        ho                 ! head at start of timestep (m)        (supp)

    real(r_2), dimension(mlon,mlat) ::  &
        h,              &  ! head, after ghmcompute (m)           (ret)
        convgw             ! convergence due to gw flow (m/s)     (ret)

    real(r_2)  :: ebot, eocn
    integer ::  istep,i
!       eocn  = mean spurious sink for h_ocn = sealev fmlon (m/s)(ret)
!               This equals the total ground-water flow across 
!               land->ocean boundaries.
!       ebot  = mean spurious source for "bot" fmlon (m/s) (returned)
!       time  = elapsed time from start of run (sec)
!       dels = timestep length (sec)
!       istep = timestep counter

! Local arrays:
    real(r_2),  dimension(mlon,mlat)   :: sf2    ! storage coefficient (m3 of h2o / bulk m3)
    real(r_2),  dimension(mlon,mlat,2) ::   t    ! transmissivity (m2/s)..1 for N-S,..2 for E-W
    real(r_2),  dimension(0:mlon+mlat) :: b,g    ! work arrays


    real(r_2), parameter    :: botinc = 0.01  ! re-wetting increment to fmlon h < bot
!   parameter (botinc = 0.  )  ! re-wetting increment to fmlon h < bot
                                 ! (m); else no flow into dry cells
    real(r_2), parameter    :: delskip = 0.005 ! av.|dhead| value for iter.skip out(m)
    integer, parameter :: itermax = 25    ! maximum number of iterations
    integer, parameter :: itermin = 3     ! minimum number of iterations
    real(r_2), parameter    :: sealev = -1.     ! sea-level evlation (m)
      
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

   LOGICAL :: debug,KeepLooping
      
   debug = .false.  

   !allocate(evl(mlon,mlat))
   !allocate(bot(mlon,mlat))
   !allocate(hycond(mlon,mlat))
   !allocate(poros(mlon,mlat))
   !allocate(compres(mlon,mlat))
   !allocate(ho(mlon,mlat))
   !allocate(h(mlon,mlat))
   !allocate(convgw(mlon,mlat))
   !allocate(sf2(mlon,mlat))
   !allocate(t(mlon,mlat,2))
   !allocate(b(0:mlon+mlat))
   !allocate(g(0:mlon+mlat))
        
    !map the 1D vector of ssnow to the 2D matrmlon that is longitude x latitude
    !hold ocean points to h=0
      
      !procedure
      !DO i = 1, mland ! over all land grid points
            ! Write to temporary variable (area weighted average across all
            ! patches):
            !otmp3xyt(land_x(i), land_y(i), 1) 
            
            !land_x(i),land_y(i) does the mapping.  make sure these are available.  use io_vars

    if (debug) write(*,*) 'about to initialize'
    if (debug) write(*,*) mlon,mlat
    evl(:,:)     = 0._r_2
    poros(:,:)  = 1._r_2
    ho(:,:)     = 0._r_2
    hycond(:,:) = 0._r_2
    if (debug) write(*,*) 'init loop through mp'

    do i=1,mp
      evl(land_x(i),land_y(i))    = soil%elevation(i)
      poros(land_x(i),land_y(i))  = soil%GWwatsat(i)
      ho(land_x(i),land_y(i))     = soil%elevation(i) - ssnow%wtd(i)/1000._r_2
      hycond(land_x(i),land_y(i)) = 0.000005*soil%GWhksat(i)/1000._r_2   !m/s
      if (ssnow%GWwb(i) .le. 1e-10) hycond(land_x(i),land_y(i)) = 0._r_2
    end do
    if (debug) write(*,*) 'done with mp init'
    compres = 0.0_r_2
    convgw  = 0._r_2

    bot(1,:) = 0.!evl(1,:) - 50.
    bot(:,1) = 0.!evl(:,1) - 50.
    bot(mlat,:) = 0.!evl(mlat,:) - 50.
    bot(:,mlon) = 0.!evl(:,mlon) - 50.
    do j=2,mlat-1
       do i=2,mlon-1
         bot(i,j) = 0.!sum(evl(i-1:i+1,j-1:j+1))/9._r_2 - 50._r_2
       end do
     end do

    where(bot .lt. 0.0) bot  = 0.0
    where(evl .le. 0.1) evl = 10.0

    where (poros .lt. 0.0) poros = 0.01
    where (poros .ge. 1.0) poros = 0.9
    where (ho .lt. evl) ho = 0.99*evl
    where(hycond .lt. 0.0) hycond = 0._r_2

    h = ho

    if (debug) write(*,*) 'ELEVATION MAX IS ',maxval(soil%elevation)

    if (debug) write(*,*) 'ELEVATION MAX IS ',maxval(evl)
    if (debug) write(*,*) 'ELEVATION MIN IS ',minval(evl)
    if (debug) write(*,*) 'ELEVATION AVG IS ',sum(evl)/real(xdimsize*ydimsize)
    if (debug) write(*,*) 'hycond max ', maxval(hycond)
    if (debug) write(*,*) 'hycond min ', minval(hycond)


    dx = soil%delx                           !assumed constant
    dy = soil%dely                           !need to add dx to global params
    darea = dx*dy
!       Top of iterative loop for ADI solution
    ContinueLoop = .true.
    iter = 0


      iter = 0
      KeepLooping = .True.
!~~~~~~~~~~~~~
!   80 continue
    MainLoop: do while (KeepLooping)  
!~~~~~~~~~~~~~
      iter = iter+1


      e  = 0.0       ! absolute changes in head (for iteration control)
!       Set storage coefficient (sf2)
      tareal = 0.


!#OMP PARALLEL DO PRIVATE(j,i,su,sc,shp)
      do j=1,mlat
        do i=1,mlon
         tareal = tareal + darea
          su = poros(i,j)                    ! new (volug)
          sc = 1._r_2
 
          if      (ho(i,j).le.evl(i,j) .and. h(i,j).le.evl(i,j)) then
            sf2(i,j) = su
          else if (ho(i,j).ge.evl(i,j) .and. h(i,j).ge.evl(i,j)) then
            sf2(i,j) = sc
          else if (ho(i,j).le.evl(i,j) .and. h(i,j).ge.evl(i,j)) then
            shp = sf2(i,j) * (h(i,j) - ho(i,j))
            sf2(i,j) = shp * sc / (shp - (su-sc)*(evl(i,j)-ho(i,j)))
          else if (ho(i,j).ge.evl(i,j) .and. h(i,j).le.evl(i,j)) then
            shp = sf2(i,j) * (ho(i,j) - h(i,j))
            sf2(i,j) = shp * su / (shp + (su-sc)*(ho(i,j)-evl(i,j)))
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
                        hycond(i, j)*(min(h(i ,j),evl(i ,j))-bot(i ,j))  &
                       *hycond(ip,j)*(min(h(ip,j),evl(ip,j))-bot(ip,j))  &
                         )    )                                           &
                   * (0.5*(dy+dy)) & 
                   / (0.5*(dx+dx))

          t(i,j,1) = sqrt( abs(                                           &
                        hycond(i,j )*(min(h(i,j ),evl(i,j ))-bot(i,j ))  &
                       *hycond(i,jp)*(min(h(i,jp),evl(i,jp))-bot(i,jp))  &
                         )    )                                           &
                   * (0.5*(dx+dx))  &
                   / (0.5*(dy+dy))

        enddo
      enddo
!#OMP END PARALLEL DO

      b(:) = 0.0
      g(:) = 0.0

!-------------------
      LonLoopOne: do ii=1,mlon
!-------------------
        i=ii
        if (mod(istep+iter,2).eq.1) i=mlon-i+1

!          calculate b and g arrays

!>>>>>>>>>>>>>>>>>>>>
        LatLoopOne: do j=1,mlat
!>>>>>>>>>>>>>>>>>>>>
          bb = (sf2(i,j)/dels) * darea
          dd = ( ho(i,j)*sf2(i,j)/dels ) * darea
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

          w    = bb - aa*b(j-1)
          b(j) = cc / w
          g(j) = (dd -aa*g(j-1))/w
!>>>>>>>>>>>>>>>
        end do LatLoopOne
!>>>>>>>>>>>>>>>

        !calc the heads again
        e = e + abs(h(i,mlat)-g(mlat))
        h(i,mlat) = g(mlat)

        do n=mlat-1,1,-1
          ha = g(n) - b(n)*h(i,n+1)
          e = e + abs(ha - h(i,n))
          h(i,n) = ha
        end do


!-------------
   end do LonLoopOne
!-------------
!=======================
!       Row calculations
!=======================

!       set transmissivities (same as above)
      do j=1,mlat
        jp = min (j+1,mlat)
        do i=1,mlon
          ip = min (i+1,mlon)
          t(i,j,2) = sqrt( abs(                                             &
                        hycond(i, j)*(min(h(i ,j),evl(i ,j))-bot(i ,j))    &
                       *hycond(ip,j)*(min(h(ip,j),evl(ip,j))-bot(ip,j))    &
                         )    )                                             &
!                    * (0.5*(dy(i,j)+dy(ip,j)))                               &
!                    / (0.5*(dx(i,j)+dx(ip,j)))
                   * (0.5*(dy+dy))                               &
                   / (0.5*(dx+dx))

          t(i,j,1) = sqrt( abs(                                             &
                        hycond(i,j )*(min(h(i,j ),evl(i,j ))-bot(i,j ))    &
                       *hycond(i,jp)*(min(h(i,jp),evl(i,jp))-bot(i,jp))    &
                         )    )                                             &
                   * (0.5*(dx+dx))                               &
                   / (0.5*(dy+dy))
        enddo
      enddo

      b(:) = 0.0
      g(:) = 0.0

!-------------------
      LatLoopTwo: do jj=1,mlat
!-------------------
        j=jj
        if (mod(istep+iter,2).eq.1) j = mlat-j+1
!         calculate b and g arrays
!>>>>>>>>>>>>>>>>>>>>
        LonLoopTwo: do i=1,mlon
!>>>>>>>>>>>>>>>>>>>>
          bb = (sf2(i,j)/dels) * darea
          dd = ( ho(i,j)*sf2(i,j)/dels ) * darea
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
          b(i) = cc / w
          g(i) = (dd - aa * g(i-1))/w

!>>>>>>>>>>>>>>>
        end do LonLoopTwo
!>>>>>>>>>>>>>>>
!          re-estimate heads
        e = e + abs(h(mlon,j)-g(mlon))
        h(mlon,j) = g(mlon)
        n = mlon-1

        e = e + abs(h(mlon,j)-g(mlon))
        h(mlon,j) = g(mlon)
        do n=mlon-1,1,-1
          ha = g(n) - b(n)*h(n+1,j)
          e = e + abs(h(n,j)-ha)
          h(n,j) = ha
        end do

!-------------
      end do LatLoopTwo
!-------------
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

      if (delcur .le. delskip*dels .or. iter .gt. itermax) then
          KeepLooping = .False.
      end if
      
  end do MainLoop


      if (debug) write(*,*) 'done with iterations'

!  Rate of convergence from ground water flow (returned)

    do j=1,mlat
      do i=1,mlon
        convgw(i,j) = sf2(i,j) * (h(i,j)-ho(i,j)) / dels
      enddo
    enddo

      if (debug) write(*,*) 'map to vectors'
      !map the 2D convergence of water to the 1D cable array
    do i=1,mp
      ssnow%GWconvergence(i) = convgw(land_x(i),land_y(i)) / 1000._r_2 * darea * dels         ![mm]
      ssnow%wtd(i)           = max((evl(land_x(i),land_y(i)) - h(land_x(i),land_y(i)))* 1000.0,0.01)  ![mm]
    end do



      
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

      if (xx(i) .ge. 0._r_2) then
         vol_av(i) = max(soil%GWwatsat(i) - ssnow%GWwb(i),0._r_2)
      else
         vol_av(i) = 0.9*ssnow%GWwb(i)
      end if 
      mss_av(i) = vol_av(i)*soil%GWdz(i)*C%denliq

      if (abs(xx(i)) .le. mss_av(i)) then
        ssnow%GWwb(i) = ssnow%GWwb(i) + xx(i)/(soil%GWdz(i)*C%denliq)
        xx(i) = 0._r_2
      elseif (abs(xx(i)) .gt. mss_av(i)) then
        if (xx(i) .ge. 0._r_2) then
           xx(i) = xx(i) - mss_av(i)
           ssnow%GWwb(i) = soil%GWwatsat(i)
        else
           xx(i) = xx(i) + mss_av(i)
           ssnow%GWwb(i) = (ssnow%GWwb(i)-vol_av(i))
        end if
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

        elseif (xx(i) .le. 1e-7) then

          vol_av(i) = 0.9*ssnow%wb(i,k)
          mss_av(i) = vol_av(i) * soil%zse(k) * C%denliq

          if (abs(xx(i)) .le. mss_av(i)) then
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k) + xx(i) / (soil%zse(k)*C%denliq)
            xx(i) = 0._r_2
          elseif (abs(xx(i)) .gt. mss_av(i)) then
            xx(i) = xx(i) + mss_av(i)
            ssnow%wbliq(i,k) = ssnow%wbliq(i,k) - mss_av(i) / (soil%zse(k)*C%denliq)
          end if


        end if
      end do   !soil layer loop

    end do  !mp patch loop



    where(xx(:) .gt. 0._r_2)   !if column is fully saturated gets added to surface runoff
      ssnow%rnof1 = ssnow%rnof1 + xx(:) / dels
      ssnow%runoff = ssnow%rnof2 + ssnow%rnof1  !adjust the total runoff as well
    end where

      


  end subroutine update_GW

!!======================================================================!!
!!                   step_gw_model                                      !!
!!======================================================================!!

  subroutine mpi_step_gw_model(dels,ssnow,soil)!,LonBg,LonEd,LatBg,LatEd)! dx,              &
            !ltype, evl, bot,        &
            !hycond, poros, compres,  &
            !ho, h, convgw,           &
            ! ebot, eocn,              &
            ! dt, istep)

    implicit none

    real, intent(in)                         :: dels
    type(soil_snow_type), intent(inout)      :: ssnow
    type(soil_parameter_type), intent(inout) :: soil
    !integer, intent(in)                      :: LonBg,&
    !                                            LonEd,&
    !                                            LatBg,&
    !                                            LatEd
    !LOCAL variables to map ssnow to previous code
    !note assume continantal??
    real(r_2),  dimension(mlon,mlat) ::  &
        evl,           &  ! evl/bathymetry of sfc rel to sl (m) (supp)
        bot,            &  ! evl. aquifer bottom rel to sl (m)   (supp)
        hycond,         &  ! hydraulic conductivity (m/s per m/m) (supp)
        poros,          &  ! porosity (m3/m3)                     (supp)
        compres,        &  ! compressibility (1/Pa)               (supp)
        ho                 ! head at start of timestep (m)        (supp)

    real(r_2), dimension(mlon,mlat) ::  &
        h,              &  ! head, after ghmcompute (m)           (ret)
        hn,             &  ! head at time t+1 (m)
        hs,             &  ! head at time t+0.5 (m)
        convgw             ! convergence due to gw flow (m/s)     (ret)

    real(r_2)  :: ebot, eocn
    integer ::  istep
!       eocn  = mean spurious sink for h_ocn = sealev fmlon (m/s)(ret)
!               This equals the total ground-water flow across 
!               land->ocean boundaries.
!       ebot  = mean spurious source for "bot" fmlon (m/s) (returned)
!       time  = elapsed time from start of run (sec)
!       dels = timestep length (sec)
!       istep = timestep counter

! Local arrays:
    real(r_2),  dimension(mlon,mlat)   :: sf2    ! storage coefficient (m3 of h2o / bulk m3)
    real(r_2),  dimension(mlon,mlat,2) ::   t    ! transmissivity (m2/s)..1 for N-S,..2 for E-W
    real(r_2),  dimension(0:mlon+mlat) :: b,g    ! work arrays
    real(r_2),  dimension(mlon,mlat)   :: xs_runoff  !runoff generated from water table less than zero

    real(r_2), parameter    :: botinc = 0.01  ! re-wetting increment to fmlon h < bot
!   parameter (botinc = 0.  )  ! re-wetting increment to fmlon h < bot
                                 ! (m); else no flow into dry cells
    real(r_2), parameter    :: delskip = 0.005 ! av.|dhead| value for iter.skip out(m)
    integer, parameter      :: itermax = 10    ! maximum number of iterations
    integer, parameter      :: itermin = 3     ! minimum number of iterations
    real(r_2), parameter    :: sealev = -1.     ! sea-level evlation (m)
      
    logical  :: ContinueLoop

    integer ::                &
        iter,                   &
        j,i,                    &
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

   LOGICAL :: debug,KeepLooping

   integer                                   :: LonBg,&
                                                LonEd,&
                                                LatBg,&
                                                LatEd

   LonBg=1;LatBg=1;
   LonEd=mlon;LatEd=mlat 
   debug = .false.  
    !map the 1D vector of ssnow to the 2D matrmlon that is longitude x latitude
    !hold ocean points to h=0
      
      !procedure
      !DO i = 1, mland ! over all land grid points
            ! Write to temporary variable (area weighted average across all
            ! patches):
            !otmp3xyt(land_x(i), land_y(i), 1) 
            
            !land_x(i),land_y(i) does the mapping.  make sure these are available.  use io_vars

    if (debug) write(*,*) 'about to initialize'
    if (debug) write(*,*) mlon,mlat
    evl(:,:)     = 0._r_2
    poros(:,:)  = 1._r_2
    ho(:,:)     = 0._r_2
    hycond(:,:) = 0._r_2
    if (debug) write(*,*) 'init loop through mp'

    !evl and poros should really be set outside this 
    !routine as they are invariant
    do i=1,mp
      evl(land_x(i),land_y(i))    = soil%elevation(i)
      poros(land_x(i),land_y(i))  = soil%GWwatsat(i)
      ho(land_x(i),land_y(i))     = soil%elevation(i) - ssnow%wtd(i)/1000._r_2
      hycond(land_x(i),land_y(i)) = soil%GWhksat(i)/1000._r_2   !m/s
    end do
    if (debug) write(*,*) 'done with mp init'
    compres   = 0.0_r_2
    convgw    = 0._r_2
    bot       = evl - 250.0
    xs_runoff = 0._r_2

    where(bot .lt. 0.0) bot  = 0.0
    where(evl .le. 0.1) evl = 10.0

    where (poros .lt. 0.0) poros = 0.01
    where (poros .ge. 1.0) poros = 0.9
    where (ho .lt. bot) ho = 1.01*evl
    where(hycond .lt. 0.0) hycond = 0._r_2

    h = ho

    if (debug) write(*,*) 'ELEVATION MAX IS ',maxval(soil%elevation)
    if (debug) write(*,*) 'ELEVATION MAX IS ',maxval(evl)
    if (debug) write(*,*) 'ELEVATION MIN IS ',minval(evl)
    if (debug) write(*,*) 'ELEVATION AVG IS ',sum(evl)/real(xdimsize*ydimsize)
    if (debug) write(*,*) 'hycond max ', maxval(hycond)
    if (debug) write(*,*) 'hycond min ', minval(hycond)


    dx = soil%delx*1000.                           !assumed constant
    dy = soil%dely*1000.                           !need to add dx to global params
    darea = dx*dy

    call compute_storage(h,ho,poros,evl,sf2)


    call TwoDGW_XDirCalc(LatBg,LatEd,dels,hycond,evl,dx,dy,sf2,h,hs)

    call compute_storage(hs,h,poros,evl,sf2)

    call TwoDGW_YDirCalc(LonBg,LonEd,dels,hycond,evl,dx,dy,sf2,h,hs,hn)


    do j=1,mlat
      do i=1,mlon
        convgw(i,j) = sf2(i,j) * (h(i,j)-ho(i,j)) / dels
      enddo
    enddo

    do i=1,mp
      ssnow%GWconvergence(i) = convgw(land_x(i),land_y(i)) / 1000._r_2 * darea * dels         ![mm]
      ssnow%wtd(i)           = max((evl(land_x(i),land_y(i)) - h(land_x(i),land_y(i)))* 1000.0,0.01)  ![mm]
    end do




  end subroutine mpi_step_gw_model


!!======================================================================!!
!!                     compute storage coefficient                      !!
!!======================================================================!!

subroutine compute_storage(h,ho,poros,evl,sf2) 
   real(r_2), intent(in)     :: h(:,:),ho(:,:),evl(:,:),poros(:,:)
   real(r_2), intent(inout)  :: sf2(:,:)

    where(ho .lt. evl .and. h .lt. evl)  &
          sf2 = poros

    where(ho .ge. evl .and. h .ge. evl)  &
          sf2 = 1._r_2

    where(ho .lt. evl .and. h .ge. evl)  & 
          sf2 = sf2*(h-ho)*1._r_2/ (sf2*(h-ho) - (poros-1._r_2)*(evl-ho))

    where(ho .ge. evl .and. h .lt. evl)  &
          sf2 = sf2*(ho-h)*poros / (sf2*(ho-h) + (poros - 1._r_2)*(ho-evl))

end subroutine compute_storage    

!!======================================================================!!
!!                   column calculation subroutine                      !!
!!======================================================================!!

  subroutine TwoDGW_XDirCalc(LatBg,LatEd,dels,K,evl,dx,dy,sf2,h,hs)
    implicit none

    integer,      intent(in)                  :: LatBg,LatEd
    real,         intent(in)                  :: dels  
    real(r_2),    intent(in),  dimension(:,:) :: K
    real(r_2),    intent(in),  dimension(:,:) :: h
    real(r_2),    intent(in),  dimension(:,:) :: evl   
    real(r_2),    intent(in)                  :: dx,dy
    real(r_2),    intent(in),  dimension(:,:) :: sf2     
    real(r_2),    intent(out), dimension(:,:) :: hs



    !Local variables
    integer                      :: i,im,ip
    integer                      :: j,jm,jp
    real(r_2), dimension(mlon)   :: at,bt,ct,rt

    at(:) = 0._r_2
    bt(:) = 0._r_2
    ct(:) = 0._r_2
    rt(:) = 1._r_2
    do j=LatBg,LatEd
      jp = min(j+1,mlat)
      jm = max(j-1,1)

      do i=1,mlon
        ip = min(i+1,mlon)
        im = max(i-1,1)

        at(i) = -trans(K(im,j),K(i,j),h(im,j),h(i,j)) !-0.5*(K(i,j)*h(i,j) + K(im,j)*h(im,j))/dx/dx
        ct(i) = -trans(K(i,j),K(ip,j),h(i,j),h(ip,j)) !-0.5*(K(ip,j)*h(ip,j) + K(i,j)*h(i,j))/dx/dx
        !bt(i) = 2.0*sf2(i,j)/dels - 0.5*( -(K(ip,j)*h(ip,j)+K(i,j)*h(i,j))/dx - &
        !        (K(i,j)*h(i,j)+K(im,j)*h(im,j))/dx)

        bt(i) = 2.0*sf2(i,j)/dels - (-trans(K(i,j),K(im,j),h(i,j),h(im,j))/dx - &
                trans(K(i,j),K(ip,j),h(i,j),h(ip,j))/dx)      

        !rt(i) = 0.5*( (K(i,jp)*h(i,jp)+K(i,j)*h(i,j))*(h(i,jp)-h(i,j))/dy - &
        !             (K(i,j)*h(i,j)+K(i,jm)*h(i,jm))*(h(i,j)-h(i,jm))/dy)/dy + &
        !            2.0*sf2(i,j)*h(i,j)/dels
        rt(i) = 0.5*( trans(K(i,jp),K(i,j),h(i,jp),h(i,j))*(h(i,jp)-h(i,j))/dy - &
                      trans(K(i,jm),K(i,j),h(i,jm),h(i,j))*(h(i,j)-h(i,jm))/dy)/dy + &
                    2.0*sf2(i,j)*h(i,j)/dels

      end do  !loop over the longitudes

      call solve_1d_tridiag(at,bt,ct,rt,hs(:,j),mlon)

    end do   !loop over a chunk of the latitudes

    !bounds checking to check h < 0, h > elevation
  



  end subroutine TwoDGW_XdirCalc  
!==================================================================!

!!======================================================================!!
!!                      row calculation subroutine                      !!
!!======================================================================!!

  subroutine TwoDGW_YDirCalc(LonBg,LonEd,dels,K,evl,dx,dy,sf2,h,hs,hn)
    implicit none

    integer,      intent(in)                  :: LonBg,LonEd
    real,         intent(in)                  :: dels  
    real(r_2),    intent(in),  dimension(:,:) :: K
    real(r_2),    intent(in),  dimension(:,:) :: h
    real(r_2),    intent(in),  dimension(:,:) :: evl
    real(r_2),    intent(in)                  :: dx,dy        
    real(r_2),    intent(in),  dimension(:,:) :: sf2
    real(r_2),    intent(in),  dimension(:,:) :: hs
    real(r_2),    intent(out), dimension(:,:) :: hn

    !Local variables
    integer                      :: i,im,ip
    integer                      :: j,jm,jp
    real(r_2), dimension(mlat)   :: at,bt,ct,rt

    at(:) = 0._r_2
    bt(:) = 0._r_2
    ct(:) = 0._r_2
    rt(:) = 1._r_2


    do i=LonBg,LonEd
      ip = min(i+1,mlon)
      im = max(i-1,1)

      do j=1,mlat
        jp = min(j+1,mlat)
        jm = max(j-1,1)

        at(j) = -trans(K(i,j),K(i,jm),hs(i,j),hs(i,jm))/dy/dy! -0.5*(K(i,j)*hs(i,j)+K(i,jm)*hs(i,jm))/dy/dy
        ct(j) = -trans(K(i,j),K(i,jp),hs(i,j),hs(i,jp))/dy/dy! -0.5*(K(i,jp)*hs(i,jp)+K(i,j)*hs(i,j))/dy/dy
        !bt(j) = 2.0*sf2(i,j)/dels +  &
        !        0.5* (K(i,jp)*hs(i,jp)+K(i,j)*hs(i,j) + &
        !              K(i,j)*hs(i,j) + K(i,jm)*hs(i,jm))/dy/dy  

        bt(j) = 2.0*sf2(i,j)/dels +  &
                (trans(K(i,j),K(i,jp),hs(i,j),hs(i,jp)) + &
                 trans(K(i,j),K(i,jm),hs(i,j),hs(i,jm)))/dy/dy

        !rt(j) = 2.0*sf2(i,j)*hs(i,j)/dels  + 0.5/dx/dx * &
        !        ( (K(ip,j)*h(ip,j)+K(i,j)*h(i,j))*(hs(ip,j)-hs(i,j)) - &
        !          (K(i,j)*h(i,j) + K(im,j)*h(im,j))*(hs(i,j)-hs(im,j)))


        rt(j) = 2.0*sf2(i,j)*hs(i,j)/dels  + &
                ( trans(K(ip,j),K(i,j),h(ip,j),h(i,j))*(hs(ip,j)-hs(i,j)) - &
                  trans(K(im,j),K(i,j),h(im,j),h(i,j))*(hs(i,j)-hs(im,j)))/dx/dx


      end do  !loop over the longitudes

      call solve_1d_tridiag(at,bt,ct,rt,hn(i,:),mlat)

    end do   !loop over a chunk of the latitudes


    !bounds checking to check h < 0, h > elevation
  


  end subroutine TwoDGW_YdirCalc  
!=====================================================================!


  function trans(ka,kb,ha,hb) result(tr)
    implicit none

    real(r_2), intent(in) :: ka,kb,ha,hb
    real(r_2)             :: tr

    tr = sqrt(abs(ka*ha*kb*hb))
    !tr = 0.5_r_2 * (k(1)*h(1) + k(2)*h(2))

  end function trans



!=====================================================================!
! SUBROUTINE solve_1d_tridiag
! solves tridiagonal set of linear equations.  returns the solution
!  
  subroutine solve_1d_tridiag (at, bt, ct, rt, ut,n)

! !ARGUMENTS:
    implicit none
    real(r_2), intent(in)    :: at(:)    ! 1 left of diagonal 
    real(r_2), intent(in)    :: bt(:)    ! diagonal 
    real(r_2), intent(in)    :: ct(:)    ! 1 right of diagonal 
    real(r_2), intent(in)    :: rt(:)    ! right hand side
    real(r_2), intent(inout) :: ut(:)    ! solution to the system of eqs
    integer, intent(in)      :: n
!   local variables
    integer  :: k
    REAL(r_2), DIMENSION(n) ::  gam  
    REAL(r_2)                  ::  bet 



    ! Solve the matrix
    if (abs(bt(1)) .gt. 1e-7) then
      bet = bt(1)
    else
      bet = bt(1) + sign(1e-7,bt(1))
    end if

    ut(1) = rt(1) / bet
    do k = 2,n
       gam(k) = ct(k-1) / bet
       bet    = max(bt(k) - at(k) * gam(k),0.00001_r_2)
       ut(k)  = (rt(k) - at(k)*ut(k-1)) / bet
    end do

    do k = n-1,1,-1
       ut(k) = ut(k) - gam(k+1) * ut(k+1)
    end do

  end subroutine solve_1d_tridiag
!=====================================================================!





      
      
end module cable_TwoDim_GW
