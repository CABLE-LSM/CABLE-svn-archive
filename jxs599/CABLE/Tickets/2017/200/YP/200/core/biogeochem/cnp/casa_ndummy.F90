SUBROUTINE casa_ndummy(casapool)
  IMPLICIT NONE
  TYPE (casa_pool),             INTENT(INOUT) :: casapool

  casapool%Nplant(:,:) = casapool%Cplant(:,:) * casapool%ratioNCplant(:,:)

END SUBROUTINE casa_ndummy


