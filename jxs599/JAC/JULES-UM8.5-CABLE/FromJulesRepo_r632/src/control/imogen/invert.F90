#if !defined(UM_JULES)
!-----------------------------------------------------------------------
! This routine is designed to invert a tri-diagonal matrix
! and then test that all is OK.
!
! The matrix is of the form ``A u_new = B u_old + k''
! Matrix A and B are tri-diagonal
!
! The solution is du/dz - pu = lambda
!-----------------------------------------------------------------------
  SUBROUTINE INVERT(                                              &
    U_OLD,U_NEW,P,LAMBDA_OLD,LAMBDA_NEW,R1,R2,DZ_TOP,N_OLEVS      &
  )

    USE logging_mod, ONLY : log_fatal

    IMPLICIT NONE

    INTEGER ::                                                    &
      N_OLEVS  !IN Number of layers in main routine above

    REAL ::                                                       &
      P,                                                          &
               !IN P Value in mixed-boundary condition
      LAMBDA_OLD,                                                 &
               !IN Lambda value at old timestep
      LAMBDA_NEW,                                                 &
               !IN Lambda value at new timestep
      R1(N_OLEVS),                                                &
               !IN Mesh ratio 
      R2(N_OLEVS),                                                &
               !IN Mesh ratio 
      DZ_TOP,                                                     &
               !IN Top layer mesh size
      FACTOR,                                                     &
               !WORK Local dummy variable
      DUMMY1,                                                     &
               !WORK Local dummy variable
      DUMMY2   !WORK Local dummy variable

    REAL ::                                                       &
      A(N_OLEVS,N_OLEVS),                                         &
               !WORK Matrix A
      B(N_OLEVS,N_OLEVS),                                         &
               !WORK Matrix B
      K(N_OLEVS),                                                 &
               !WORK Matrix K
      U_OLD(1:N_OLEVS),                                           &
               !IN Old values of U
      U_NEW(1:N_OLEVS)
               !OUT New values of U
      
    REAL ::                                                       &
      A_L(N_OLEVS,N_OLEVS),                                       &
               !WORK A local working value of A - values ma
      F_L(N_OLEVS)
               !WORK Local working value of ``B u_old + k''

    INTEGER :: I,J     !WORK Looping parameters



! Set everything to zero in the first instance 
    A(:,:) = 0.0
    B(:,:) = 0.0
    A_L(:,:) = 0.0
    F_L(:) = 0.0
    K(:) = 0.0


! Evaluate the values of A,B,K for the particular use here.
    A(1,1) = (1.0 + R1(1)) + R1(1) * DZ_TOP * P
    A(1,2) = -R1(1)

    B(1,1) = (1.0 - R1(1)) - R1(1) * DZ_TOP * P
    B(1,2) = R1(1)

    K(1)   = -R1(1) * DZ_TOP * (LAMBDA_OLD + LAMBDA_NEW) 

    DO I = 2,N_OLEVS-1
      A(I,I-1) = -R1(I)
      A(I,I)   = 1.0 + R1(I) + R2(I)
      A(I,I+1) = -R2(I)

      B(I,I-1) = R1(I)
      B(I,I)   = 1.0 - R1(I) - R2(I)
      B(I,I+1) = R2(I)
    ENDDO

    A(N_OLEVS,N_OLEVS-1) = -R1(N_OLEVS)
    A(N_OLEVS,N_OLEVS)   = 1.0 + R1(N_OLEVS) + R2(N_OLEVS)

    B(N_OLEVS,N_OLEVS-1) = R1(N_OLEVS)
    B(N_OLEVS,N_OLEVS)   = 1.0 - R1(N_OLEVS) - R2(N_OLEVS)


! First evaluate F_L (which is initially `` B u_old + k'')
    F_L(1) = B(1,1) * U_OLD(1) + B(1,2) * U_OLD(2) + K(1)  

    DO I = 2,N_OLEVS-1 
      F_L(I) = B(I,I-1)*U_OLD(I-1) + B(I,I)*U_OLD(I)              &
             + B(I,I+1)*U_OLD(I+1) + K(I)     
    ENDDO
 
    F_L(N_OLEVS) = B(N_OLEVS,N_OLEVS-1) * U_OLD(N_OLEVS-1)        &
                 + B(N_OLEVS,N_OLEVS)*U_OLD(N_OLEVS) + K(N_OLEVS)


! Now set A_L = A     
    A_L(1,1) = A(1,1)
    A_L(1,2) = A(1,2)

    DO I = 2,N_OLEVS-1
      A_L(I,I-1) = A(I,I-1)
      A_L(I,I)   = A(I,I)
      A_L(I,I+1) = A(I,I+1)
    ENDDO

    A_L(N_OLEVS,N_OLEVS-1) = A(N_OLEVS,N_OLEVS-1)
    A_L(N_OLEVS,N_OLEVS)   = A(N_OLEVS,N_OLEVS)

! Now go through loops working back to find U_NEW(1)
    FACTOR = A_L(N_OLEVS-1,N_OLEVS)/A_L(N_OLEVS,N_OLEVS)
    A_L(N_OLEVS,N_OLEVS-1) = A_L(N_OLEVS,N_OLEVS-1) * FACTOR
    A_L(N_OLEVS,N_OLEVS) = A_L(N_OLEVS-1,N_OLEVS)
    F_L(N_OLEVS) = F_L(N_OLEVS) * FACTOR

    A_L(N_OLEVS-1,N_OLEVS-1) = A_L(N_OLEVS-1,N_OLEVS-1)           &
                             - A_L(N_OLEVS,N_OLEVS-1)
    A_L(N_OLEVS-1,N_OLEVS) = 0.0
    F_L(N_OLEVS-1) = F_L(N_OLEVS-1) - F_L(N_OLEVS)
      
    DO I = N_OLEVS-1,2,-1
      FACTOR = A_L(I-1,I)/A_L(I,I)       
      A_L(I,I-1) = A_L(I,I-1) * FACTOR
      A_L(I,I) = A_L(I-1,I)
      F_L(I) = F_L(I) * FACTOR

      A_L(I-1,I-1) = A_L(I-1,I-1) - A_L(I,I-1)       
      A_L(I-1,I) = 0.0
      F_L(I-1) = F_L(I-1) - F_L(I)
    ENDDO

! Now explicitly calculate the U_NEW values
    U_NEW(1) = F_L(1) / A_L(1,1)

    DO I = 2,N_OLEVS
      U_NEW(I) = (F_L(I) - (A_L(I,I-1)*U_NEW(I-1)))/A_L(I,I)
    ENDDO

! Now perform a check
    DO I = 1,N_OLEVS
      DUMMY1 = 0.0
      DUMMY2 = 0.0
      DO J = 1,N_OLEVS
        DUMMY1 = DUMMY1 + A(I,J)*U_NEW(J)
        DUMMY2 = DUMMY2 + B(I,J)*U_OLD(J)
      ENDDO
      DUMMY2 = DUMMY2 + K(I) 
        
      IF(ABS(DUMMY1-DUMMY2) >= 0.0001) THEN 
        CALL log_fatal("INVERT", "Invert unsuccessful")
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE INVERT
#endif
