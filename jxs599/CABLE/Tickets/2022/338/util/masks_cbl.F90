MODULE cbl_masks_mod

IMPLICIT NONE

PUBLIC L_tile_pts
PUBLIC fveg_mask
PUBLIC fsunlit_mask
PUBLIC fsunlit_veg_mask

!mask TRUE where tile fraction is greater than zero
LOGICAL, ALLOCATABLE :: L_tile_pts(:,:)

CONTAINS

SUBROUTINE fveg_mask( veg_mask, mp, lai_thresh, reducedLAIdue2snow )

IMPLICIT NONE
LOGICAL, INTENT(OUT),  ALLOCATABLE :: veg_mask(:)

INTEGER, INTENT(IN) :: mp
REAL,    INTENT(IN) :: lai_thresh
REAL,    INTENT(IN) :: reducedLAIdue2snow(mp)
!local vars
INTEGER :: i

IF ( .NOT. ALLOCATED(veg_mask)) ALLOCATE( veg_mask(mp) )

veg_mask(:) = .FALSE.
! Define vegetation mask:
DO i=1, mp
  IF ( reducedLAIdue2snow(i) > LAI_thresh ) veg_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fveg_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fsunlit_mask( sunlit_mask, mp, coszen_tols, coszen )

IMPLICIT NONE
LOGICAL, INTENT(OUT), ALLOCATABLE :: sunlit_mask(:)

INTEGER, INTENT(IN) :: mp
REAL,    INTENT(IN) :: coszen_tols
REAL,    INTENT(IN) :: coszen(mp)
!local vars
INTEGER :: i

IF ( .NOT. ALLOCATED(sunlit_mask)) ALLOCATE( sunlit_mask(mp) )

sunlit_mask = .FALSE.
! Define sunlit mask:
DO i=1, mp
  IF ( coszen(i) > coszen_tols ) sunlit_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fsunlit_mask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE fsunlit_veg_mask(  sunlit_veg_mask, veg_mask, sunlit_mask, mp )

IMPLICIT NONE
LOGICAL, INTENT(OUT), ALLOCATABLE :: sunlit_veg_mask(:)

INTEGER, INTENT(IN) :: mp
LOGICAL, INTENT(IN) :: veg_mask(mp)
LOGICAL, INTENT(IN) :: sunlit_mask(mp)
!local vars
INTEGER :: i

IF ( .NOT. ALLOCATED(sunlit_veg_mask)) ALLOCATE( sunlit_veg_mask(mp) )

sunlit_veg_mask = .FALSE.
! Define sunlit AND vegetation mask:
DO i=1, mp
  IF ( veg_mask(i) .AND.  sunlit_mask(i) ) sunlit_veg_mask(i) = .TRUE.
END DO

RETURN
END SUBROUTINE fsunlit_veg_mask

END MODULE cbl_masks_mod
