! cable_abort.f90
!
! Source file containing abort routines for CABLE LSM
!
! Development by Gab Abramowitz, Harvey Davies
!
! bugs to bernard.pak@csiro.au
!
! This file contains subroutines:
!   abort
!   nc_abort


MODULE abort_module
  IMPLICIT NONE
CONTAINS
  !========================================================================
  SUBROUTINE abort(message)
    CHARACTER(LEN=*), INTENT(IN) :: message
    WRITE(*, *) message
    STOP 1
  END SUBROUTINE abort
  !========================================================================
  SUBROUTINE nc_abort(ok,message)
    USE netcdf
    ! Error subroutine in case of fatal error:
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER, INTENT(IN) :: ok
    WRITE(*,*) message ! error from subroutine
    WRITE(*,*) NF90_STRERROR(ok) ! netcdf error details
    STOP
  END SUBROUTINE nc_abort
  !========================================================================
  SUBROUTINE range_abort(message,ktau,met,value,var_range,i,xx,yy)
    USE define_types, ONLY: met_type
    USE define_dimensions, ONLY: r_1,i_d
    USE io_variables, ONLY: latitude,longitude,landpt,lat_all,lon_all
    ! Abort subroutine which also prints internal netcdf error message.
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER(i_d), INTENT(IN) :: ktau ! time step
    TYPE(met_type),INTENT(IN) :: met  ! met data
    REAL(r_1),INTENT(IN) :: value ! value deemed to be out of range
    INTEGER(i_d),INTENT(IN) :: i ! landpt number of erroneous grid square
    INTEGER(i_d),INTENT(IN),OPTIONAL:: xx ! coordinates of erroneous grid square
    INTEGER(i_d),INTENT(IN),OPTIONAL ::yy ! coordinates of erroneous grid square
    REAL(r_1),DIMENSION(2),INTENT(IN) :: var_range ! appropriate var range 
    WRITE(*,*) message ! error from subroutine
    IF(PRESENT(yy)) THEN ! i.e. using rectangular land/sea grid
      WRITE(*,*) 'Site lat, lon:',lat_all(xx,yy),lon_all(xx,yy)
      WRITE(*,*) 'Output timestep',ktau, &
         ', or ', met%hod(landpt(i)%cstart),' hod, ',&
         INT(met%doy(landpt(i)%cstart)),'doy, ',INT(met%year(landpt(i)%cstart))
    ELSE ! i.e. using compressed land only grid
      WRITE(*,*) 'Site lat, lon:',latitude(i),longitude(i)
      WRITE(*,*) 'Output timestep',ktau, &
         ', or ', met%hod(landpt(i)%cstart),' hod, ',&
         INT(met%doy(landpt(i)%cstart)),'doy, ',INT(met%year(landpt(i)%cstart))
    END IF
    WRITE(*,*) 'Specified acceptable range (checks.f90):', var_range(1), &
         'to',var_range(2)
    WRITE(*,*) 'Value:',value 
    STOP
  END SUBROUTINE range_abort
  !========================================================================
END MODULE abort_module
