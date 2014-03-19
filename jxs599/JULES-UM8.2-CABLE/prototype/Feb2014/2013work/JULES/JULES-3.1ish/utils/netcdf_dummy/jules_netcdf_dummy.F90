#if !defined(UM_JULES) && defined(NCDF_DUMMY)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/utils/netcdf_dummy/jules_netcdf_dummy.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

MODULE netcdf

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   This MODULE provides a dummy INTERFACE to NetCDF constants and PROCEDUREs
!   that can be used to link against. If any attempt is made to use the
!   PROCEDUREs, an error is given
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER ::                                                       &
    nf90_int    = 4,                                                          &
    nf90_float  = 5,                                                          &
    nf90_double = 6

  INTEGER, PARAMETER ::                                                       &
    nf90_nowrite   = 0,                                                       &
    nf90_write     = 1,                                                       &
    nf90_clobber   = 0,                                                       &
    nf90_noclobber = 4,                                                       &
    nf90_fill      = 0,                                                       &
    nf90_nofill    = 256

  INTEGER, PARAMETER :: nf90_unlimited = 0

  INTEGER, PARAMETER :: nf90_global = 0

  INTEGER, PARAMETER :: nf90_noerr = 0

! Overloaded variable FUNCTIONs
  INTERFACE nf90_def_var
    MODULE PROCEDURE nf90_def_var_Scalar, nf90_def_var_oneDim, nf90_def_var_ManyDims
  END INTERFACE ! nf90_def_var

! Overloaded attribute FUNCTIONs
  INTERFACE nf90_put_att
    MODULE PROCEDURE nf90_put_att_text, nf90_put_att_int, nf90_put_att_real
  END INTERFACE

  INTERFACE nf90_get_att
    MODULE PROCEDURE nf90_get_att_text, nf90_get_att_int, nf90_get_att_real
  END INTERFACE

! Overloaded variable FUNCTIONs
  INTERFACE nf90_put_var
    MODULE PROCEDURE nf90_put_var_text, nf90_put_var_int, nf90_put_var_real
    MODULE PROCEDURE nf90_put_var_1D_text, nf90_put_var_1D_int,               &
                     nf90_put_var_1D_real
    MODULE PROCEDURE nf90_put_var_2D_text, nf90_put_var_2D_int,               &
                     nf90_put_var_2D_real
    MODULE PROCEDURE nf90_put_var_3D_text, nf90_put_var_3D_int,               &
                     nf90_put_var_3D_real
    MODULE PROCEDURE nf90_put_var_4D_text, nf90_put_var_4D_int,               &
                     nf90_put_var_4D_real
    MODULE PROCEDURE nf90_put_var_5D_text, nf90_put_var_5D_int,               &
                     nf90_put_var_5D_real
    MODULE PROCEDURE nf90_put_var_6D_text, nf90_put_var_6D_int,               &
                     nf90_put_var_6D_real
    MODULE PROCEDURE nf90_put_var_7D_text, nf90_put_var_7D_int,               &
                     nf90_put_var_7D_real
  END INTERFACE ! nf90_put_var

  INTERFACE nf90_get_var
    MODULE PROCEDURE nf90_get_var_text, nf90_get_var_int, nf90_get_var_real
    MODULE PROCEDURE nf90_get_var_1D_text, nf90_get_var_1D_int,               &
                     nf90_get_var_1D_real
    MODULE PROCEDURE nf90_get_var_2D_text, nf90_get_var_2D_int,               &
                     nf90_get_var_2D_real
    MODULE PROCEDURE nf90_get_var_3D_text, nf90_get_var_3D_int,               &
                     nf90_get_var_3D_real
    MODULE PROCEDURE nf90_get_var_4D_text, nf90_get_var_4D_int,               &
                     nf90_get_var_4D_real
    MODULE PROCEDURE nf90_get_var_5D_text, nf90_get_var_5D_int,               &
                     nf90_get_var_5D_real
    MODULE PROCEDURE nf90_get_var_6D_text, nf90_get_var_6D_int,               &
                     nf90_get_var_6D_real
    MODULE PROCEDURE nf90_get_var_7D_text, nf90_get_var_7D_int,               &
                     nf90_get_var_7D_real
  END INTERFACE ! nf90_get_var

CONTAINS

  FUNCTION nf90_strerror(ncerr)
    INTEGER, INTENT( in) :: ncerr
    CHARACTER(len = 80)  :: nf90_strerror

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  END FUNCTION nf90_strerror

  FUNCTION nf90_enddef(ncid, h_minfree, v_align, v_minfree, r_align)
    INTEGER,           INTENT( in) :: ncid
    INTEGER, OPTIONAL, INTENT( in) :: h_minfree, v_align, v_minfree, r_align
    INTEGER                        :: nf90_enddef

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  END FUNCTION nf90_enddef

  FUNCTION nf90_close(ncid)
    INTEGER, INTENT( in) :: ncid
    INTEGER              :: nf90_close

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  END FUNCTION nf90_close

  FUNCTION nf90_inquire(ncid, nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum)
    INTEGER,           INTENT( in) :: ncid
    INTEGER, OPTIONAL, INTENT(out) :: nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum
    INTEGER                        :: nf90_inquire

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  END FUNCTION nf90_inquire

  function nf90_open(path, mode, ncid, chunksize, cache_size, cache_nelems,   &
                                                 cache_preemption, comm, info)
    implicit none
    character (len = *), intent(in) :: path
    integer, intent(in) :: mode
    integer, intent(out) :: ncid
    integer, optional, intent(inout) :: chunksize
    integer, optional, intent(in) :: cache_size, cache_nelems
    real, optional, intent(in) :: cache_preemption
    integer, optional, intent(in) :: comm, info
    integer :: nf90_open

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_open

  function nf90_create(path, cmode, ncid, initialsize, chunksize, cache_size, &
                                   cache_nelems, cache_preemption, comm, info)
    implicit none
    character (len = *), intent(in) :: path
    integer, intent(in) :: cmode
    integer, intent(out) :: ncid
    integer, optional, intent(in) :: initialsize
    integer, optional, intent(inout) :: chunksize
    integer, optional, intent(in) :: cache_size, cache_nelems
    integer, optional, intent(in) :: cache_preemption
    integer, optional, intent(in) :: comm, info
    integer :: nf90_create

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_create

  function nf90_def_dim(ncid, name, len, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent( in) :: len
    integer,             intent(out) :: dimid
    integer                          :: nf90_def_dim

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_def_dim

  function nf90_inq_dimid(ncid, name, dimid)
    integer,             intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer,             intent(out) :: dimid
    integer                          :: nf90_inq_dimid

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_inq_dimid

  function nf90_inquire_dimension(ncid, dimid, name, len)
    integer,                       intent( in) :: ncid, dimid
    character (len = *), optional, intent(out) :: name
    integer,             optional, intent(out) :: len
    integer                                    :: nf90_inquire_dimension

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_inquire_dimension

  function nf90_put_att_text(ncid, varid, name, values)
    integer,                          intent( in) :: ncid, varid
    character(len = *),               intent( in) :: name
    character(len = *),               intent( in) :: values
    integer                                       :: nf90_put_att_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_att_text

  function nf90_get_att_text(ncid, varid, name, values)
    integer,                          intent( in) :: ncid, varid
    character(len = *),               intent( in) :: name
    character(len = *),               intent(out) :: values
    integer                                       :: nf90_get_att_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_att_text

  function nf90_put_att_int(ncid, varid, name, values)
    integer,            intent( in) :: ncid, varid
    character(len = *), intent( in) :: name
    integer,            intent( in) :: values
    integer :: nf90_put_att_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_att_int

  function nf90_get_att_int(ncid, varid, name, values)
    integer,            intent( in) :: ncid, varid
    character(len = *), intent( in) :: name
    integer,            intent(out) :: values
    integer :: nf90_get_att_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_att_int

  function nf90_put_att_real(ncid, varid, name, values)
    integer,            intent( in) :: ncid, varid
    character(len = *), intent( in) :: name
    real,               intent( in) :: values
    integer :: nf90_put_att_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_att_real

  function nf90_get_att_real(ncid, varid, name, values)
    integer,            intent( in) :: ncid, varid
    character(len = *), intent( in) :: name
    real,               intent(out) :: values
    integer :: nf90_get_att_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_att_real

  function nf90_def_var_Scalar(ncid, name, xtype, varid)
    integer, intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer, intent(in) :: xtype
    integer, intent(out) :: varid
    integer :: nf90_def_var_Scalar

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_def_var_Scalar

  function nf90_def_var_oneDim(ncid, name, xtype, dimids, varid, contiguous,  &
                  chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                  cache_size, cache_nelems, cache_preemption)
    integer, intent( in) :: ncid
    character (len = *), intent( in) :: name
    integer, intent(in) :: xtype
    integer, intent(in) :: dimids
    integer, intent(out) :: varid
    logical, optional, intent(in) :: contiguous
    integer, optional, intent(in) :: chunksizes
    integer, optional, intent(in) :: deflate_level
    logical, optional, intent(in) :: shuffle, fletcher32
    integer, optional, intent(in) :: endianness
    integer, optional, intent(in) :: cache_size, cache_nelems, cache_preemption
    integer :: nf90_def_var_oneDim

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_def_var_oneDim

  function nf90_def_var_ManyDims(ncid, name, xtype, dimids, varid, contiguous,&
       chunksizes, deflate_level, shuffle, fletcher32, endianness, cache_size,&
       cache_nelems, cache_preemption)
    integer, intent(in) :: ncid
    character (len = *), intent(in) :: name
    integer, intent( in) :: xtype
    integer, dimension(:), intent(in) :: dimids
    integer, intent(out) :: varid
    logical, optional, intent(in) :: contiguous
    integer, optional, dimension(:), intent(in) :: chunksizes
    integer, optional, intent(in) :: deflate_level
    logical, optional, intent(in) :: shuffle, fletcher32
    integer, optional, intent(in) :: endianness
    integer, optional, intent(in) :: cache_size, cache_nelems, cache_preemption
    integer :: nf90_def_var_ManyDims

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_def_var_ManyDims

  function nf90_inq_varid(ncid, name, varid)
    integer, intent(in) :: ncid
    character (len = *), intent( in) :: name
    integer, intent(out) :: varid
    integer :: nf90_inq_varid

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_inq_varid

  function nf90_inquire_variable(ncid, varid, name, xtype, ndims, dimids,     &
                                 nAtts, contiguous, chunksizes, deflate_level,&
                                 shuffle, fletcher32, endianness, cache_size, &
                                 cache_nelems, cache_preemption)
    integer, intent(in) :: ncid, varid
    character (len = *), optional, intent(out) :: name
    integer, optional, intent(out) :: xtype, ndims
    integer, dimension(:), optional, intent(out) :: dimids
    integer, optional, intent(out) :: nAtts
    logical, optional, intent(out) :: contiguous
    integer, optional, dimension(:), intent(out) :: chunksizes
    integer, optional, intent(out) :: deflate_level
    logical, optional, intent(out) :: shuffle, fletcher32
    integer, optional, intent(out) :: endianness
    integer, optional, intent(out) :: cache_size, cache_nelems, cache_preemption
    integer :: nf90_inquire_variable

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_inquire_variable

  function nf90_put_var_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *),             intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_text

  function nf90_get_var_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *),             intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_text

  function nf90_put_var_1D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_1D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_1D_text

  function nf90_put_var_2D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_2D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_2D_text

  function nf90_put_var_3D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_3D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_3D_text

  function nf90_put_var_4D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_4D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_4D_text

  function nf90_put_var_5D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_5D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_5D_text

  function nf90_put_var_6D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_6D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_6D_text

  function nf90_put_var_7D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_7D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_7D_text

  function nf90_get_var_1D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_1D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_1D_text

  function nf90_get_var_2D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_2D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_2D_text

  function nf90_get_var_3D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_3D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_3D_text

  function nf90_get_var_4D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_4D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_4D_text

  function nf90_get_var_5D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_5D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_5D_text

  function nf90_get_var_6D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_6D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_6D_text

  function nf90_get_var_7D_text(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    character (len = *), dimension(:, :, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_7D_text

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_7D_text


  function nf90_put_var_int(ncid, varid, values, start)
    integer, intent( in) :: ncid, varid
    integer, intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer :: nf90_put_var_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_int

  function nf90_put_var_real(ncid, varid, values, start)
    integer, intent( in) :: ncid, varid
    real, intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer :: nf90_put_var_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_real

  function nf90_get_var_int(ncid, varid, values, start)
    integer, intent( in) :: ncid, varid
    integer, intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer :: nf90_get_var_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_int

  function nf90_get_var_real(ncid, varid, values, start)
    integer, intent( in) :: ncid, varid
    real, intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start
    integer :: nf90_get_var_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_real

  function nf90_put_var_1D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_1D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_1D_int

  function nf90_put_var_2D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_2D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_2D_int

  function nf90_put_var_3D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_3D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_3D_int

  function nf90_put_var_4D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_4D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_4D_int

  function nf90_put_var_5D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_5D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_5D_int

  function nf90_put_var_6D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_6D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_6D_int

  function nf90_put_var_7D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_put_var_7D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_7D_int

  function nf90_put_var_1D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_1D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_1D_real

  function nf90_put_var_2D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_2D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_2D_real

  function nf90_put_var_3D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_3D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_3D_real

  function nf90_put_var_4D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_4D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_4D_real

  function nf90_put_var_5D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_5D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_5D_real

  function nf90_put_var_6D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_6D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_6D_real

  function nf90_put_var_7D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :, :, :), intent( in) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer  :: nf90_put_var_7D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_put_var_7D_real

  function nf90_get_var_1D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_1D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_1D_int

  function nf90_get_var_2D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_2D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_2D_int

  function nf90_get_var_3D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_3D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_3D_int

  function nf90_get_var_4D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_4D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_4D_int

  function nf90_get_var_5D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_5D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_5D_int

  function nf90_get_var_6D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_6D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_6D_int

  function nf90_get_var_7D_int(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    integer, dimension(:, :, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_7D_int

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_7D_int

  function nf90_get_var_1D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_1D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_1D_real

  function nf90_get_var_2D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_2D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_2D_real

  function nf90_get_var_3D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_3D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_3D_real

  function nf90_get_var_4D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_4D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_4D_real

  function nf90_get_var_5D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_5D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_5D_real

  function nf90_get_var_6D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_6D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_6D_real

  function nf90_get_var_7D_real(ncid, varid, values, start, count, stride, map)
    integer, intent( in) :: ncid, varid
    real, dimension(:, :, :, :, :, :, :), intent(out) :: values
    integer, dimension(:), optional, intent( in) :: start, count, stride, map
    integer :: nf90_get_var_7D_real

    WRITE(*,*) 'FATAL ERROR: Attempt to use dummy NetCDF procedure'
    WRITE(*,*) 'To use NetCDF, recompile linking the NetCDF library'
    STOP
  end function nf90_get_var_7D_real

END MODULE netcdf
#endif
