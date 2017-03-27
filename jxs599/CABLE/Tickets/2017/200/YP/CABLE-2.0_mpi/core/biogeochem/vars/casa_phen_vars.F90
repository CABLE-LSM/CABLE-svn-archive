MODULE phenvariable
  USE casadimension
  IMPLICIT NONE
  TYPE phen_variable
    INTEGER,   DIMENSION(:),  POINTER :: phase        
    REAL(r_2), DIMENSION(:),  POINTER :: TKshed
    INTEGER,   DIMENSION(:,:),POINTER :: doyphase
  END type phen_variable

CONTAINS

SUBROUTINE alloc_phenvariable(phen,arraysize)
!SUBROUTINE alloc_phenvariable(phen,arraysize,mvt)
  IMPLICIT NONE
  TYPE(phen_variable), INTENT(INOUT) :: phen
  INTEGER,             INTENT(IN) :: arraysize
!  INTEGER,        INTENT(IN) :: mvt

  ALLOCATE(phen%Tkshed(mvtype))
!  ALLOCATE(phen%Tkshed(mvt))
  ALLOCATE(phen%phase(arraysize),         &
           phen%doyphase(arraysize,mphase))
END SUBROUTINE alloc_phenvariable

End MODULE phenvariable

