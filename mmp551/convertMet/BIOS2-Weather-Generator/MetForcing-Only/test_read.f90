! Jatin: program to test read one of the binary files
PROGRAM test_read
IMPLICIT NONE
INTEGER :: iunit=10, dy=841 ,dx=681, iostatus
CHARACTER(LEN=200) :: filename,basedir
REAL, DIMENSION(:), ALLOCATABLE :: data_out

filename=TRIM('/short/w35/jtk561/tempo_BOM/tmin/20000102_tmin.flt')
! the *4 is needed, i think because the total files size divide by (dx*dy) gives
! me a factor of 4 exactly, this seems to be a machine related thing....
OPEN(iunit, file=filename, form='unformatted',status='old',action='read',access='direct',recl=dx*dy*4,iostat=iostatus)
IF (iostatus .NE. 0) THEN
   PRINT *, 'Problem opening file'
   STOP
END IF

ALLOCATE(data_out(dx*dy))
READ(unit=iunit,iostat=iostatus,rec=1) data_out
IF (iostatus .NE. 0) THEN
   PRINT *, 'Problem reading file'
   STOP
END IF

PRINT *, data_out

DEALLOCATE(data_out)

END PROGRAM test_read
