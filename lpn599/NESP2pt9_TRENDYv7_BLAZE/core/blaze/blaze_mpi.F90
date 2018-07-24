MODULE BLAZE_MPI

  USE MPI
  USE cable_mpicommon
  USE cable_def_types_mod, ONLY: ncp
  USE BLAZE

  ! Total number of restart parameters for BLAZE
  INTEGER, PARAMETER :: n_blaze_restart = 9

  ! Total number of output parameters for BLAZE
  INTEGER, PARAMETER :: n_blaze_output = 0

  ! Total number of restart parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_restart = 0

  ! Total number of input parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_input = 0

  ! Total number of output parameters for SIMFIRE
  INTEGER, PARAMETER :: n_simfire_output = 0

  
CONTAINS


SUBROUTINE master_blaze_types (comm, wland, mp, blaze_restart_ts, BLAZE)

  ! Send blaze restart data to workers  
  
  IMPLICIT NONE

  INTEGER              , INTENT(IN)  :: comm ! MPI communicator to talk to the workers
  TYPE(lpdecomp_t)     , INTENT(IN)  :: wland
  INTEGER              , INTENT(IN)  :: mp   ! Number of gridcells
  TYPE(TYPE_BLAZE)     , INTENT(IN)  :: BLAZE
  INTEGER, DIMENSION(:), INTENT(OUT) :: blaze_restart_ts,blaze_ts 

  ! MPI: temp arrays for marshalling all types into a struct
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blocks
  INTEGER(KIND=MPI_ADDRESS_KIND), ALLOCATABLE, DIMENSION(:) :: displs
  INTEGER, ALLOCATABLE, DIMENSION(:) :: types
  INTEGER :: ntyp ! number of worker's types

  INTEGER :: last2d, i

  ! MPI: block lenghts for hindexed representing all vectors
  INTEGER, ALLOCATABLE, DIMENSION(:) :: blen

  ! MPI: block lengths and strides for hvector representing matrices
  INTEGER :: r1len, r2len
  INTEGER(KIND=MPI_ADDRESS_KIND) :: r1stride, r2stride

  INTEGER :: tsize, totalrecv, totalsend, mp
  INTEGER(KIND=MPI_ADDRESS_KIND) :: text, tmplb

  INTEGER :: rank, off, cnt
  INTEGER :: bidx, midx, vidx, ierr


  ! Restart value handles (restart_blaze_ts()) 
  
  ntyp = n_blaze_restart

  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))

  mp = LEN(
  
  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx = 0

     ! ------------- 2D arrays -------------

     ! Annual (daily) rainfall (ncells,366)
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AnnRainf(off,1), displs(bidx), ierr) ! 1
     CALL MPI_Type_create_hvector (ndoy, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground life woody biomass 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr) ! 2
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground life grassy biomass 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr) ! 3
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground woody litter
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr) ! 4
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground grassy litter 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr) ! 5
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx
     
     ! ------------- 1D vectors -------------

     ! Integer days since last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
     blocks(bidx) = i1len

     ! Real(sp) Last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! current KBDI
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! DEADWOOD
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! ------------- Wrap up -------------

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_restart_ts(rank), ierr)
     CALL MPI_Type_commit (blaze_restart_ts(rank), ierr)

     CALL MPI_Type_size (blaze_restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (blaze_restart_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)

  !=============================================================================
  ! Standard desired Output (blaze_ts()) 
  !=============================================================================

!   INTEGER,  DIMENSION(:),  ALLOCATABLE :: DSLR,Flix
!   REAL,     DIMENSION(:),  ALLOCATABLE :: RAINF, KBDI, LR, U10,RH,TMAX,TMIN
!   REAL,     DIMENSION(:),  ALLOCATABLE :: F,FLI,ROS,Z,D,w,DFLI,AB
!   REAL,     DIMENSION(:,:),ALLOCATABLE :: TO, AGC_g, AGC_w

  IF( TRIM(BLAZE%OUTPUT) == "full") THEN
     ntyp = n_blaze_output_full
  ELSE
     ntyp = n_blaze_output
  ENDIF
  
  ALLOCATE (blocks(ntyp))
  ALLOCATE (displs(ntyp))
  ALLOCATE (types(ntyp))
  
  istride  = mp * extid ! short integer 
  r1stride = mp * extr1 ! single precision
  r2stride = mp * extr2 ! double precision

  DO rank = 1, wnp
     off = wland(rank)%patch0
     cnt = wland(rank)%npatch

     i1len = cnt * extid  
     r1len = cnt * extr1
     r2len = cnt * extr2

     bidx = 0

     ! ------------- 2D arrays -------------

     ! Annual (daily) rainfall (ncells,366)
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AnnRainf(off,1), displs(bidx), ierr) ! 1
     CALL MPI_Type_create_hvector (ndoy, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground life woody biomass 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLB_w(off,1), displs(bidx), ierr) ! 2
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground life grassy biomass 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLB_g(off,1), displs(bidx), ierr) ! 3
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground woody litter
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_w(off,1), displs(bidx), ierr) ! 4
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     ! Above ground grassy litter 
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%AGLit_g(off,1), displs(bidx), ierr) ! 5
     CALL MPI_Type_create_hvector (ncp, r1len, r1stride, MPI_BYTE, &
     &                             types(bidx), ierr)
     blocks(bidx) = 1

     last2d = bidx
     
     ! ------------- 1D vectors -------------

     ! Integer days since last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DSLR(off), displs(bidx), ierr)
     blocks(bidx) = i1len

     ! Real(sp) Last rainfall
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%LR(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! current KBDI
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%KBDI(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! DEADWOOD
     bidx = bidx + 1
     CALL MPI_Get_address (BLAZE%DEADWOOD(off), displs(bidx), ierr)
     blocks(bidx) = r1len

     ! ------------- Wrap up -------------

     types(last2d+1:bidx) = MPI_BYTE

     ! MPI: sanity check
     IF (bidx /= ntyp) THEN
        WRITE (*,*) 'invalid blaze ntyp constant, fix it!'
        CALL MPI_Abort (comm, 1, ierr)
     END IF

     CALL MPI_Type_create_struct (bidx, blocks, displs, types, blaze_restart_ts(rank), ierr)
     CALL MPI_Type_commit (blaze_restart_ts(rank), ierr)

     CALL MPI_Type_size (blaze_restart_ts(rank), tsize, ierr)
     CALL MPI_Type_get_extent (blaze_restart_ts(rank), tmplb, text, ierr)

     WRITE (*,*) 'restart results recv from worker, size, extent, lb: ', &
   &       rank,tsize,text,tmplb

     totalrecv = totalrecv + tsize

     ! free the partial types used for matrices
     DO i = 1, last2d
        CALL MPI_Type_free (types(i), ierr)
     END DO

  END DO

  WRITE (*,*) 'total size of restart fields received from all workers: ', totalrecv

  ! MPI: check whether total size of received data equals total
  ! data sent by all the workers
  totalsend = 0
  CALL MPI_Reduce (MPI_IN_PLACE, totalsend, 1, MPI_INTEGER, MPI_SUM, &
    &     0, comm, ierr)

  WRITE (*,*) 'total size of restart fields sent by all workers: ', totalsend

  IF (totalrecv /= totalsend) THEN
          WRITE (*,*) 'error: restart fields totalsend and totalrecv differ'
          CALL MPI_Abort (comm, 0, ierr)
  END IF

  DEALLOCATE(types)
  DEALLOCATE(displs)
  DEALLOCATE(blocks)


 END SUBROUTINE master_blaze_types



END MODULE BLAZE_MPI
