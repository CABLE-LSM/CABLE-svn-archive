MODULE cbl_rhoch_module

IMPLICIT NONE

PUBLIC :: calc_rhoch
PRIVATE

CONTAINS

! this subroutine called from _init_radiation on cable_albedo.F90 pathway, explicit and implict
SUBROUTINE calc_rhoch( c1,rhoch, mp, nrb, taul, refl )

IMPLICIT NONE

!model dimensions
INTEGER, INTENT(IN) :: mp                     ! number of "tiles"
INTEGER, INTENT(IN) :: nrb                    ! # of rad. bands VIS,NIR(,LW)

REAL, INTENT(OUT) :: c1(mp,nrb)
REAL, INTENT(OUT) :: rhoch(mp,nrb)            ! REFLection black horiz leaves
REAL, INTENT(IN) :: taul(mp,nrb)              ! Leaf Transmisivity
REAL, INTENT(IN) :: refl(mp,nrb)              ! Leaf Reflectivity

c1(:,1) = SQRT(1.0 - taul(:,1) - refl(:,1))
c1(:,2) = SQRT(1.0 - taul(:,2) - refl(:,2))
c1(:,3) = 1.0

! Canopy C%REFLection black horiz leaves
! (eq. 6.19 in Goudriaan and van Laar, 1994):
rhoch = (1.0 - c1) / (1.0 + c1)

RETURN
END SUBROUTINE calc_rhoch

END MODULE cbl_rhoch_module
