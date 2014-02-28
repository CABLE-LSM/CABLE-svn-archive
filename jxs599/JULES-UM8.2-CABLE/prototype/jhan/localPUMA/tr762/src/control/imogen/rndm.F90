#if !defined(UM_JULES)
!**********************************************************************
!
! This routine calculates the random behaviour. Following the ITE
! model, the a random number generator is not called - instead,
! seed(1),..,seed(4) is updated.
!
!**********************************************************************

  SUBROUTINE RNDM(RANDOM_NUM,SEED)

    IMPLICIT NONE

    REAL :: RANDOM_NUM
              !OUT The random number.

    INTEGER ::                                                    &
      SEED(4),                                                    &
              !IN/OUT The seeding numbers
      I       !WORK Integer



    SEED(4) = 3 * SEED(4) + SEED(2)
    SEED(3) = 3 * SEED(3) + SEED(1)
    SEED(2) = 3 * SEED(2)
    SEED(1) = 3 * SEED(1)

    I = INT(SEED(1) / 1000.0)
    SEED(1) = SEED(1) - I * 1000
    SEED(2) = SEED(2) + I

    I = INT(SEED(2) / 100.0)
    SEED(2) = SEED(2) - 100 * I
    SEED(3) = SEED(3) + I

    I = INT(SEED(3) / 1000.0)
    SEED(3) = SEED(3) - I * 1000
    SEED(4) = SEED(4) + I

    I = INT(SEED(4) / 100.0)
    SEED(4) = SEED(4) - 100 * I

    RANDOM_NUM = (((FLOAT(SEED(1)) * 0.001 + FLOAT(SEED(2)))      &
                  * 0.01 + FLOAT(SEED(3))) * 0.001                &
                  + FLOAT(SEED(4))) * 0.01

    RETURN

  END SUBROUTINE RNDM
#endif
