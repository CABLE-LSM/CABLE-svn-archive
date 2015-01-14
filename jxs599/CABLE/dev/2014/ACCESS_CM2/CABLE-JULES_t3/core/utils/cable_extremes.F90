

MODULE cable_extremes_module
   IMPLICIT NONE
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
 
   interface cable_extremes
      module procedure cable_extremes1, cable_extremes2
   end interface cable_extremes
   
CONTAINS

SUBROUTINE cable_extremes1(fname,field,mype)

   real, dimension(:,:) :: field
   character(len=*), dimension(:) :: fname
   
   integer, optional :: mype 
   integer :: i,j,k
   integer :: n,m,op
   real :: emax, emin, emean, emode
   real :: erange 
   real :: edbin 
   real, dimension(100) :: bin
   integer, dimension(100) :: ibin
   integer :: ib, ibmax, binmax, maxbin

   n = size(fname)
   m = size(field,2)

   ! for each field in fname(i)
   do i=1, n  
      
      emax =  MAXVAL( field(i,:) )
      emin =  MINVAL(field(i,:) )
      emean =  SUM(field(i,:) ) / ( m )
    
      erange = emax - emin
      edbin = erange / 100. ! for 100 bins

      bin(1) = emin          

      ! define bins per fname(i)
      do ib=2, 100
         bin(ib) = bin(ib-1) + edbin
      enddo

      ibin =0
      ! for each Element in field 
      do j=1, m  
      
         ! Assignn each Element to a bin 
         do ib=1, 99

            IF( field(i,j) >= bin(ib) .AND. &
                field(i,j) < bin(ib+1) ) THEN
               
               !if(ib==1 )print *, "jhan:field1 ", field(i,j)
               ibin(ib) = ibin(ib) + 1
            
            ENDIF   
          
         enddo ! DO LOOP over fill bins

      enddo ! DO LOOP over elements 
 
      binmax = 0 
      
      ! find max bin per field 
      do ib=1, 99

         IF( ibin(ib) > binmax ) THEN 
            binmax = ibin(ib) 
            maxbin = ib 
         ENDIF   
       
      enddo ! DO LOOP over bins
     
     print *, "jhan:bins1 ", bin 
     print *, "jhan:bins1 count", ibin 
     
      
      emode = bin(maxbin) 

         print *, ""
         print *, "CABLE_log: "
         print *, "   Field ", fname(i)
         print *, "   Min ", emin
         print *, "   Max ", emax
         print *, "   Mean ",emean
         print *, "   Mode ",emode
         print *, "End CABLE_log: "
         print *, ""
      
      enddo ! DO LOOP over fname(i)


END SUBROUTINE cable_extremes1


SUBROUTINE cable_extremes2(fname,field,mype)

   real, dimension(:,:,:) :: field
   character(len=*), dimension(:) :: fname
   
   integer, optional :: mype 
   integer :: i,j,k
   integer :: n,m,op
   real :: emax, emin, emean, emode
   real :: erange 
   real :: edbin 
   real, dimension(100) :: bin
   integer, dimension(100) :: ibin
   integer :: ib, ibmax, binmax, maxbin

   n = size(fname)
   m = size(field,2)
   op= size(field,3)

   ! for each field in fname(i)
   do i=1, n  
      
      emax =  MAXVAL( field(i,:,:) )
      emin =  MINVAL(field(i,:,:) )
      emean =  SUM(field(i,:,:) ) / ( m*op )
    
      erange = emax - emin
      edbin = erange / 100. ! for 100 bins

      bin(1) = emin          

      ! define bins per fname(i)
      do ib=2, 100
         bin(ib) = bin(ib-1) + edbin
      enddo

      ibin =0
      ! for each Element in field 
      do j=1, m  
         
         do k=1, op  
      
            ! Assignn each Element to a bin 
            do ib=1, 99

               IF( field(i,j,k) >= bin(ib) .AND. &
                   field(i,j,k) < bin(ib+1) ) THEN
                  !if(ib==1) print *, "jhan:field2 ", field(i,j,k)
                  ibin(ib) = ibin(ib) + 1
            
               ENDIF   
          
            enddo ! DO LOOP over fill bins
            
         enddo ! DO LOOP over elements 

      enddo ! DO LOOP over elements 
 
      binmax = 0 
      
      ! find max bin per field 
      do ib=1, 99

         IF( ibin(ib) > binmax ) THEN 
            binmax = ibin(ib) 
            maxbin = ib 
         ENDIF   
       
      enddo ! DO LOOP over bins
      
      emode = bin(maxbin) 

     print *, "jhan:bins2 ", bin 
     print *, "jhan:bins2 count", ibin 
     
         print *, ""
         print *, "CABLE_log: "
         print *, "   Field ", fname(i)
         print *, "   Min ", emin
         print *, "   Max ", emax
         print *, "   Mean ",emean
         print *, "   Mode ",emode
         print *, "End CABLE_log: "
         print *, ""
      
      enddo ! DO LOOP over fname(i)

END SUBROUTINE cable_extremes2

END MODULE cable_extremes_module


