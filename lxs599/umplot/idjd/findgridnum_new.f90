      program findgridnum
!
!     read header from netcdf generated by xconv to link grid numbers with
!     lat/long
!     for any input lat/long, find nearest grid number and report land fraction
!
      !use netcdf
      implicit none

      character*80, parameter:: &
            !filein ='/cs/datastore/u/csdar/ste69f/idjd/landfraction_N48.nc', &
            !filein ='/cs/datastore/u/csdar/ste69f/plot/um_scripts/ncl_scripts/landfrac_ACCESS_N96.nc', &
            !filein ='/cs/datastore/u/csdar/ste69f/plot/um_scripts/ncl_scripts/landfrac_v2.nc', &
            !filein ='/home/ste69f/idjd/landfrac_N320.nc', &
            filein ='/home/ste69f/idjd/qrparm.landfrac_vn85_n96.nc', &
            !filein  ='/cs/datastore/u/csdar/ste69f/idjd/landfraction_N96.nc', &
            sitelist='/cs/datastore/u/csdar/ste69f/idjd/sites_um_moses', &
            !fileout='/cs/datastore/u/csdar/ste69f/idjd/sites_um_moses_grid.N48'
            !fileout ='/cs/datastore/u/csdar/ste69f/idjd/sites_um_moses_grid.N96_new'
            !fileout='/cs/datastore/u/csdar/ste69f/idjd/sites_um_moses_grid.N96_v1'
            !fileout='/cs/datastore/u/csdar/ste69f/idjd/sites_um_moses_grid.N96_minland0p9'
            !fileout='/home/ste69f/idjd/sites_um_moses_grid.N320_minland0p8'
            fileout='/home/ste69f/idjd/sites_um_moses_grid.N96_vn85'
      real, parameter :: minland=0.8

      integer ok,ncid,londim,latdim,nlon,nlat,lonid,latid
      integer landid
      integer i,j,nsite

      character*10 :: tempname
      real :: inlat,inlon
      integer :: ix(1),jy(1),inalt,inland,templon,templat,count
      logical :: landcheck
      
      real, dimension(:), allocatable :: lon,lat,diflat,diflon
      real, dimension(:,:), allocatable :: land
!
      include 'netcdf.inc'
!
      ok = nf_open(trim(filein),0,ncid)
      ok = nf_inq_dimid(ncid,'longitude',londim)
      ok = nf_inq_dimid(ncid,'latitude',latdim)
      ok = nf_inq_dimlen(ncid,londim,nlon)
      ok = nf_inq_dimlen(ncid,latdim,nlat)
      allocate(lon(nlon),lat(nlat))
      ok = nf_inq_varid(ncid,'longitude',lonid)
      ok = nf_inq_varid(ncid,'latitude',latid)
      ok = nf_get_var_real(ncid,lonid,lon)
      ok = nf_get_var_real(ncid,latid,lat)
      allocate(land(nlon,nlat))
      ok = nf_inq_varid(ncid,'lsm',landid)
      ok = nf_get_var_real(ncid,landid,land)

      allocate(diflat(nlat),diflon(nlon))
!
      where(land > 2) land=0.
!
      open(unit=1,file=trim(sitelist),form='formatted')
      open(unit=2,file=trim(fileout),form='formatted')
      read(1,*) tempname
      read(1,*) nsite

      do i=1,nsite
         read(1,*) tempname,inlat,inlon,inalt,inland
         if (inlon .lt. 0) inlon = inlon + 360
         diflat = lat-inlat
         diflon = lon-inlon

         ix = minloc(abs(diflon))
         jy = minloc(abs(diflat))

         count = 1.
         templon = ix(1)
         templat = jy(1)
         landcheck = (land(templon,templat) .lt. minland)
         do while (landcheck .and. (inland .eq. 1))
            print*,'Site Shift: ',tempname
            templon = ix(1) + count           
            if (land(templon,jy(1)) .gt. minland) then
               landcheck = .FALSE.
               templat = jy(1)
               goto 100
            endif

            templon = ix(1) - count
            if (land(templon,jy(1)) .gt. minland) then
               landcheck = .FALSE.
               templat = jy(1)
               goto 100
             endif
 
            templat = jy(1) + count
            if (land(ix(1),templat) .gt. minland) then
               landcheck = .FALSE.
               templon = ix(1)
               goto 100
            endif
        
            templat = jy(1) - count
            if (land(ix(1),templat) .gt. minland) then
               landcheck = .FALSE.
               templon = ix(1)
               goto 100
            endif

            count = count + 1     
            100 continue

         enddo

         38 format(1x,i3,1x,i3,1x,a3,1x,f10.7,1x,i2,1x,f8.3,1x,f8.3,1x,f8.3,1x,f8.3) 
         !write(2,38) templon, templat, tempname,land(templon,templat),inland, &
         !print *,'Order now row,col not col,row'
         write(2,38) templat, templon, tempname,land(templon,templat),inland, &
                    lon(templon),inlon,lat(templat),inlat
      enddo

       close(unit=1)
       close(unit=2)

      end      