  subroutine landuse_data(mlon,mlat,landmask,arealand,luc_atransit,luc_fharvw,luc_xluh2cable)
  use netcdf
  use cable_abort_module,   ONLY: nc_abort
  use cable_common_module,  ONLY: filename
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax,mharvw
  IMPLICIT NONE

  integer     mlon,mlat
  real(r_2), dimension(mland,mvmax,mvmax)         :: luc_atransit
  real(r_2), dimension(mland,mharvw)              :: luc_fharvw
  real(r_2), dimension(mland,mvmax,mstate)        :: luc_xluh2cable
  integer,    dimension(mlon,mlat)                :: landmask
  real(r_2),  dimension(mland)                    :: arealand
  ! "mland" variables
  real(r_2),  dimension(:,:),      allocatable    :: areax    
  !
  integer ivt,ee,hh,np,p,q,np1
  integer ncid,ok,xID,yID,varID,i,j,m,mpx

    ! get " mlon mlat landmask" from "gridinfo"
    ok = NF90_OPEN(filename%type, 0, ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error opening grid info file.')

    ok = NF90_INQ_DIMID(ncid, 'longitude', xID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'x', xID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring x dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, xID, LEN=mlon)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting x dimension.')
    ok = NF90_INQ_DIMID(ncid, 'latitude', yID)
    IF (ok /= NF90_NOERR) ok = NF90_INQ_DIMID(ncid, 'y', yID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error inquiring y dimension.')
    ok = NF90_INQUIRE_DIMENSION(ncid, yID, LEN=mlat)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error getting y dimension.')

    allocate(areax(mlon,mlat))
     
    ok = NF90_INQ_VARID(ncid, 'area', varID)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                  'Error finding variable area')
    ok = NF90_GET_VAR(ncid, varID, areax)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok,                                    &
                  'Error reading variable longitude.')

    ok = NF90_CLOSE(ncid)
    IF (ok /= NF90_NOERR) CALL nc_abort(ok, 'Error closing grid info file.')

     m=0; arealand= 0.0
     do i=1,mlon
     do j=1,mlat
        if(areax(i,j) >0.01) then
           landmask(i,j) = 1
           m=m+1
           arealand(m) = areax(i,j)
        else
           landmask(i,j) =0
        endif
     enddo
     enddo
     if(m/=mland) then
        print *, 'mland not consistent: check gridinof area'
        stop
     endif   

     ! get the mapping matrix (landuse type to PFT)
     call landuse_getxluh2(mlat,mlon,landmask,filename%fxluh2cable,luc_xluh2cable)    !"xluh2cable"
     call landuse_getdata(mlat,mlon,landmask,filename%fxpft,luc_atransit,luc_fharvw)
  end subroutine landuse_data


SUBROUTINE landuse_getxluh2(mlat,mlon,landmask,fxluh2cable,luc_xluh2cable)
! get data: luc%fprimary; luc%fsecondary
  USE netcdf
  use cable_abort_module,   ONLY: nc_abort
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax,mharvw
  IMPLICIT NONE
  character*500 fxluh2cable
  integer   mlat,mlon
  integer,  dimension(mlon,mlat)                 :: landmask
  real(r_2), dimension(mland,mvmax,mvmax)       :: luc_atransit
  real(r_2), dimension(mland,mharvw)            :: luc_fharvw
  real(r_2),dimension(mland,mvmax,mstate)        :: luc_xluh2cable
  ! local variables
  real(r_2),   dimension(:,:,:,:), allocatable   :: xluh2cable
  integer ok,ncid2,varxid
  integer i,j,m,v,s

    allocate(xluh2cable(mlon,mlat,mvmax,mstate))
    ok = nf90_open(fxluh2cable,nf90_nowrite,ncid2)
    ok = nf90_inq_varid(ncid2,"xluh2cable",varxid)
    ok = nf90_get_var(ncid2,varxid,xluh2cable)
    ok = nf90_close(ncid2)
    ! assig the values of luc%variables
    luc_xluh2cable(:,:,:) = 0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          luc_xluh2cable(m,:,:) = xluh2cable(i,j,:,:)
          do s=1,mstate
             do v=1,mvmax
                luc_xluh2cable(m,v,s) = luc_xluh2cable(m,v,s)/sum(luc_xluh2cable(m,1:mvmax,s))
             enddo
          enddo
       endif
    enddo
    enddo

    deallocate(xluh2cable)

 END SUBROUTINE landuse_getxluh2

SUBROUTINE landuse_getdata(mlat,mlon,landmask,fxpft,luc_atransit,luc_fharvw)
! get LUC data
  USE netcdf
  use cable_abort_module,   ONLY: nc_abort
  USE cable_def_types_mod,  ONLY: mland,r_2
  use landuse_constant,     ONLY: mstate,mvmax,mharvw
  IMPLICIT NONE
  character*500 fxpft
  integer mlat,mlon
  integer,   dimension(mlon,mlat)               :: landmask 
  real(r_2), dimension(mland,mvmax,mvmax)       :: luc_atransit
  real(r_2), dimension(mland,mharvw)            :: luc_fharvw
  ! local variables
  real(r_2),  dimension(:,:,:),   allocatable   :: fracharvw
  real(r_2),  dimension(:,:,:,:), allocatable   :: transitx
  integer  ok,ncid1,varxid
  integer  i,j,m,k,ivt

    allocate(fracharvw(mlon,mlat,mharvw))
    allocate(transitx(mlon,mlat,mvmax,mvmax))

    ok = nf90_open(fxpft,nf90_nowrite,ncid1)
    ok = nf90_inq_varid(ncid1,"harvest",varxid)
    ok = nf90_get_var(ncid1,varxid,fracharvw)
    ok = nf90_inq_varid(ncid1,"transition",varxid)
    ok = nf90_get_var(ncid1,varxid,transitx)
    ok = nf90_close(ncid1)

    ! assig the values of luc%variables
    luc_fharvw(:,:) =0.0; luc_atransit(:,:,:)=0.0
    m = 0
    do i=1,mlon
    do j=1,mlat
       if(landmask(i,j) ==1) then
          m= m +1
          luc_atransit(m,:,:)   = transitx(i,j,:,:)
          luc_fharvw(m,:)       = fracharvw(i,j,:)
       endif
    enddo
    enddo
    
    deallocate(fracharvw)
    deallocate(transitx)
END SUBROUTINE landuse_getdata
