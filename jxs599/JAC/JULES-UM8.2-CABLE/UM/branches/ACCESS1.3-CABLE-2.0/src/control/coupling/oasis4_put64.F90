#if defined(OASIS4)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE OASIS4_put64(data_64,row_length,rows,ierr64  &
     & ,ft,halo_type_no_halo,gc_all_proc_group,nproc,mype)
!
! Description: This routine acts as a 64 bit interface to the 
!              put process. Currently it contains a gather
!              operation in order to accommodate bit comparable
!              running using the non-conservative OASIS4 coupler.
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
!               and UM 64 bit data.
! 
       IMPLICIT NONE

       integer  rows
       integer  row_length
       integer  ierr64, ft
       real data_64(row_length,rows,1)

       INTEGER I,J, N
      integer  rows_64
      integer  rows_32_s
      integer  rows_32_e
      integer  info64
      integer  halo_type_no_halo
      integer  gc_all_proc_group
      integer  nproc
      integer  mype


      INTEGER  GATHER_PE
      integer  ierr64, ft
      integer  ierror
      integer  info
      integer  rows_global
      character*(80) :: cmessage

      integer new_rows
      real ,allocatable,dimension(:,:)  :: data_send
      real ,allocatable,dimension(:,:,:) :: ds_chunk
       

      rows_global = rows
      CALL GC_ISUM(1, nproc, info, rows_global)

      ALLOCATE (DATA_SEND(row_length,rows_global))
 
      ! We now have to gather the data and redistribute it in
      ! accordance with what OASIS4 thinks our dimensions are
      GATHER_PE = NPROC-1

      ! 1st: Gather all onto 1 PE
! DEPENDS ON: GATHER_FIELD
      CALL GATHER_FIELD(data_64(1,1,1),data_send,                &            
     &  row_length,                                              &
     &  rows,                                                    &
     &  row_length,                                              &
     &  rows_global,                                             &
     &  ft,halo_type_no_halo,                                    &
     &  gather_pe,GC_ALL_PROC_GROUP,info,cmessage)

      ! 2nd: Broadcast the global array to everyone
      CALL GC_RBCAST(666,row_length*rows_global,GATHER_PE,       &
     &               nproc,info,DATA_SEND)

      IF (MYPE.EQ.GATHER_PE) THEN
         rows_32_s = (2*MYPE)+1
         rows_32_e = rows_global
      ELSE
         rows_32_s = (2*MYPE)+1
         rows_32_e = rows_32_s+1
      END IF

      NEW_ROWS = rows_32_e-rows_32_s+1
      ALLOCATE (DS_CHUNK(row_length,new_rows,1))

      N=0
      DO J=rows_32_s,rows_32_e
         N=N+1
         DO I = 1, row_length
            ds_chunk(I,N,1)=data_send(I,J)
         END DO
      END DO

! DEPENDS ON: OASIS4_grid32_mod
      CALL PUT32_OASIS4(ds_chunk,row_length,new_rows,ierr64)

      DEALLOCATE (DS_CHUNK)
      DEALLOCATE (DATA_SEND)

      END SUBROUTINE OASIS4_put64
#endif


