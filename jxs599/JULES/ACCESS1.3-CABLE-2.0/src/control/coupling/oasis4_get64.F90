#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
       SUBROUTINE OASIS4_get64(data_64,row_length,rows,      &
     &                         o4info,ierr64,ft              &
     & ,halo_type_no_halo,gc_all_proc_group,nproc,mype)
!
! Description: This routine acts as a 64 bit interface to the 
!              get process. Currently it contains a special
!              gather and scatter operation in order to accommodate 
!              bit comparable running using the non-conservative 
!              OASIS4 coupler.
!              Most of this functionality will eventually be 
!              removed.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.4    Jan 2007  Original code. R. Hill, R. Barnes
!================================================================

! DEPENDS ON: OASIS4_grid32_mod
       USE OASIS4_grid32_mod
!
!  Description: Buffer routine between OASIS 32 bit operations
!               and UM 64 bit data. Get incoming data.
!
       IMPLICIT NONE

       integer  rows
       integer  row_length
       integer  ierr64, o4info, ft
       real data_64(row_length,rows,1)
       real data_64_save(row_length,rows,1)

      integer  rows_64
      integer  rows_global

      integer  i, j, k,n
      integer  scatter_pe
      integer  ft,halo_type_no_halo
      integer  gc_all_proc_group
      integer  info64
      integer  nproc
      integer  mype
      integer  start_row

      integer  ierr64
      integer  ierror
      integer  o4info64, max_o4_rows
      integer  info
      real ,allocatable,dimension(:,:) :: data_recv
      real ,allocatable,dimension(:,:,:) :: data_temp

      logical, parameter :: verbose = .true.
      character*(80) :: cmessage


      data_64_save = data_64

      ! See how many rows there are in the global domain
      rows_global = rows
      CALL GC_ISUM(1, nproc, info, rows_global)

      ALLOCATE (DATA_RECV(row_length,rows_global))

      IF (MYPE.EQ.NPROC-1) THEN
         max_o4_rows = rows_global-(2*(NPROC-1))
      ELSE
         max_o4_rows = 2
      END IF

      start_row = (MYPE*2)+1

      ALLOCATE (DATA_TEMP(1:row_length,1:max_o4_rows,1))
      DATA_TEMP = 0.0 

! DEPENDS ON: OASIS4_grid32_mod
       CALL GET32_OASIS4(DATA_TEMP,row_length,max_o4_rows,    &
     &                   o4info,ierr64)

      ! Gather the data we've just GOT. The difficulty
      ! here is that we can't use gather_field because
      ! that contains lots of assumptions about the
      ! size of each domain which are not true for
      ! this data!


      ! So here's what we do.....
      ! Assuming each PE other than O4CNTLPE
      ! only deals with 2 rows, we send 2 rows of data
      ! to PE o4CNTLPE.
      DO I = 0,O4CNTLPE-1

         CALL GC_SSYNC(NPROC,INFO)
         IF (MYPE.EQ.I) THEN
            CALL GC_RSEND(I,row_length*2,o4CNTLPE,INFO,    &
     &           DATA_RECV,data_temp)
         END IF
         CALL GC_SSYNC(NPROC,INFO)
         IF (MYPE.EQ.O4CNTLPE) THEN
            CALL GC_RRECV(I,row_length*2,I,INFO,           &     
     &                 DATA_RECV(1,(I*2)+1),data_temp)
         END IF
         CALL GC_SSYNC(NPROC,INFO)

      END DO

      ! Then we fill up the rest of the array with what we
      ! have on the o4CNTLPE PE.
      IF (MYPE.EQ.o4CNTLPE) THEN
         WRITE(6,*) "ASSIGNING DATA_RECV" 
         N = 0
         DO J = start_row,rows_global
            N = N + 1
            DO I = 1, row_length
               DATA_RECV(I,J) = data_temp(I,N,1)
            END DO
         END DO
      END IF

      ! Ok we should now have a global array which
      ! contains all our received regridded data.
      ! Now we can scatter this data according to
      ! the true MPP decomposition.
      SCATTER_PE = o4CNTLPE
! DEPENDS ON: SCATTER_FIELD
      CALL SCATTER_FIELD(data_64,DATA_RECV,            &          
     &  row_length,                                    & 
     &  rows,                                          & 
     &  row_length,                                    &
     &  rows_global,                                   &
     &  ft,halo_type_no_halo,                          &
     &  scatter_pe,GC_ALL_PROC_GROUP,info,cmessage)


      IF (MYPE.EQ.0) THEN
          ! This is the Antarctic handler.
          ! Restore uncorrupted data over the most
          ! southerly 10 rows - i.e. points where
          ! we don't expect to have incoming data
          ! from the ocean. (HadGAM resolution dependent). 
          DO J = 1, 10
            DO I = 1, row_length
            data_64(I,J,1)=data_64_save(I,J,1)
            END DO
          END DO
      END IF 

      !Now we're ready, deallocate our temporary arrays
      DEALLOCATE (DATA_TEMP)
      DEALLOCATE (DATA_RECV)

      END SUBROUTINE OASIS4_get64
#endif


