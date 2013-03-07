
MODULE debug_read_mod

   IMPLICIT NONE

CONTAINS

!=============================================================================!

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
         
         ALLOCATE( newdata(dimy,dimx) )
         ALLOCATE( olddata(dimy,dimx) )
 
   END SUBROUTINE read_txt_file 

   !==========================================================================!
   !==========================================================================!

   SUBROUTINE Read_dat_file( Lfilename, ar_data, dimx, dimy )
   
      REAL, INTENT(OUT), DIMENSION(:,:) :: ar_data
      CHARACTER(LEN=*), INTENT(IN) :: Lfilename
      INTEGER, INTENT(IN) :: dimx, dimy
      
      INTEGER, PARAMETER :: gok=0
      INTEGER :: gopenstatus
      INTEGER :: i

         OPEN(UNIT=2, FILE=Lfilename//'.bin', STATUS="unknown", ACTION="read", &
                  IOSTAT=gopenstatus, FORM="unformatted" )
            
            IF(gopenstatus==gok) THEN
               
               DO i=1,dimy
                  READ(2), trash 
                  ar_data(j,i,:) = ar_Nvars( ( (j-1)*dimx )+1 : j*dimx )
               enddo
        
            ELSE
               WRITE (*,*), Lfilename//'.bin',' NOT found for read'
               STOP
        
            ENDIF
         
         CLOSE(2)

   end subroutine read_dat_file 

   !==========================================================================!
   !==========================================================================!

end module debug_read_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
