SUBROUTINE casa_pdummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Pplant(:,:) = casapool%Nplant(:,:) / casapool%ratioNPplant(:,:)

END SUBROUTINE casa_pdummy

