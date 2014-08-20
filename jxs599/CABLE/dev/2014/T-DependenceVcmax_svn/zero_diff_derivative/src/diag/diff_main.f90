
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

PROGRAM zero_diff_main 
   
   IMPLICIT NONE
   
   ! base of filenames (created by CABLE) - passed to PROGRAM 
   CHARACTER(LEN=300) :: filename1, filename2 

   ! dimx = typically # landpoints over which the var is specified per timestep 
   ! dimy = # timesteps
   INTEGER ::dimx, dimy 

   ! data arrays  
   REAL, DIMENSION(:,:), POINTER:: newdata, olddata

      IF( IARGC() > 0 ) THEN
         CALL GETARG(1, filename1)
         CALL GETARG(2, filename2)
      ELSE
         STOP 'This program requires arguments to be passed to it.'
      ENDIF

      ! read info about the spec. binary data which was created by the 
      ! host , i.e. how many points per timestep, how many timesteps.  
      CALL read_txt_file( TRIM(filename1), dimx, dimy )
      
      ALLOCATE( newdata(dimy,dimx) )
      ALLOCATE( olddata(dimy,dimx) )

      ! read the binary data andstore in 2nd arg 
      !CALL read_dat_file( TRIM(filename1), olddata, dimx, dimy  )
      CALL read_dat_file( TRIM(filename2), newdata, dimx,dimy )

      ! compute difference b/n old & new binary 
      !CALL comp_diff( olddata, newdata, dimx, dimy )

    END PROGRAM zero_diff_main 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

SUBROUTINE read_txt_file( Lfilename, dimx, dimy )

   CHARACTER(LEN=*), INTENT(IN) :: Lfilename
   INTEGER, INTENT(OUT) :: dimx, dimy
   
   INTEGER, PARAMETER :: gok=0
   INTEGER :: gopenstatus
   CHARACTER(LEN=99) :: trash 
   INTEGER :: i

      OPEN( UNIT=1,FILE=Lfilename//'.dat', STATUS="old",ACTION="read",      & 
            IOSTAT=gopenstatus )
         
         IF(gopenstatus==gok) THEN
            
            READ(1,*)
            READ (1,*), trash  
            READ (1,*), trash
            READ (1,*), trash

            READ (1,*), trash
            READ (1,*), dimx 
            READ (1,*), trash
            READ (1,*), dimy 
         
         ELSE

            WRITE (*,*), Lfilename//'.dat',' NOT found to read'
            STOP
      
         ENDIF

      CLOSE(1)
     
END SUBROUTINE read_txt_file 

!==========================================================================!

SUBROUTINE Read_dat_file( Lfilename, ar_data, dimx,dimy )

   INTEGER, INTENT(IN) :: dimx, dimy
   REAL, INTENT(OUT), DIMENSION(dimy,dimx) :: ar_data
   CHARACTER(LEN=*), INTENT(IN) :: Lfilename
   
   INTEGER, PARAMETER :: gok=0
   INTEGER :: gopenstatus
   INTEGER :: i

      OPEN(UNIT=2, FILE=Lfilename//'.bin', STATUS="unknown", ACTION="read", &
               IOSTAT=gopenstatus, FORM="unformatted" )
         
         IF(gopenstatus==gok) THEN
            
            DO i=1,dimy
               READ(2), ar_data(i,:) 
               PRINT *, 'data : ',ar_data(i,:)
            ENDDO
     
         ELSE
            WRITE (*,*), Lfilename//'.bin',' NOT found for read'
            STOP
     
         ENDIF
      
      CLOSE(2)

END SUBROUTINE read_dat_file 

!==========================================================================!
      
SUBROUTINE comp_diff( olddata, newdata, dimx, dimy )

   INTEGER, INTENT(IN) :: dimx, dimy
   REAL, INTENT(IN), DIMENSION(dimy, dimx) :: newdata, olddata
   REAL, DIMENSION(dimy,dimx) :: diff_data
   REAL :: sum_data

      diff_data = olddata - newdata
      sum_diff = SUM(diff_data)
      
      PRINT *, 'summed difference between old and new data :sum_diff'
      PRINT *, sum_diff
      PRINT *, ''

END SUBROUTINE comp_diff

!=======================================================================!






