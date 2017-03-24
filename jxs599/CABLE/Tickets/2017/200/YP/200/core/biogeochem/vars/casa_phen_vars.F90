MODULE phenvariable
  USE casadimension
  IMPLICIT NONE
  TYPE phen_variable
    INTEGER,   DIMENSION(:),  POINTER :: phase
    REAL(r_2), DIMENSION(:),  POINTER :: TKshed
    INTEGER,   DIMENSION(:,:),POINTER :: doyphase
    REAL, DIMENSION(:),  POINTER :: phen   ! fraction of max LAI
    REAL, DIMENSION(:),  POINTER :: aphen  ! annual leaf on sum
    INTEGER,   DIMENSION(:,:),POINTER :: phasespin
    INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_1
    INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_2
    INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_3
    INTEGER,   DIMENSION(:,:),POINTER :: doyphasespin_4

  END type phen_variable

CONTAINS

SUBROUTINE alloc_phenvariable(phen,arraysize)

  IMPLICIT NONE
  TYPE(phen_variable), INTENT(INOUT) :: phen
  INTEGER,             INTENT(IN) :: arraysize

  ALLOCATE(phen%Tkshed(mvtype))
  ALLOCATE(phen%phase(arraysize),         &
           phen%doyphase(arraysize,mphase))
  ALLOCATE(phen%phen(arraysize), &
       phen%aphen(arraysize), &
       phen%phasespin(arraysize,mdyear), &
       phen%doyphasespin_1(arraysize,mdyear), &
       phen%doyphasespin_2(arraysize,mdyear), &
       phen%doyphasespin_3(arraysize,mdyear), &
       phen%doyphasespin_4(arraysize,mdyear))
END SUBROUTINE alloc_phenvariable

End MODULE phenvariable

