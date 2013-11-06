      MODULE decomp

      IMPLICIT NONE

      PUBLIC

      TYPE lpdecomp_t
          INTEGER :: landp0      ! starting land point index
          INTEGER :: nland       ! number of landpoints

          INTEGER :: patch0      ! starting patch index in global CABLE vars
          INTEGER :: npatch      ! sum of patches for all landpoints of this
                                 ! worker
      END TYPE

      CONTAINS

      SUBROUTINE load_info (decompf,wnp,wland)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: decompf
      INTEGER, INTENT(OUT) :: wnp
      TYPE(lpdecomp_t), POINTER, DIMENSION(:), INTENT(OUT) :: wland

      INTEGER :: n, landp0, nland, patch0, npatch
      INTEGER :: fd, rank

      fd = 10
      OPEN (fd, FILE = decompf )

      READ (fd, *) wnp

      ALLOCATE (wland(wnp))

      DO rank = 1, wnp
         READ (fd,*) n,landp0,nland,patch0,npatch
         wland(rank)%landp0 = landp0
         wland(rank)%nland = nland
         wland(rank)%patch0 = patch0
         wland(rank)%npatch = npatch
      END DO

      CLOSE (fd)

      RETURN

      END SUBROUTINE load_info

      SUBROUTINE global2local (wnp, wland, point, global, rank, local)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: wnp
      TYPE(lpdecomp_t), POINTER, DIMENSION(:), INTENT(IN) :: wland
      LOGICAL, INTENT(IN) :: point
      INTEGER, INTENT(IN) :: global
      INTEGER, INTENT(OUT) :: rank, local

      INTEGER :: i

      DO i = 1, wnp
         SELECT CASE (point)
         CASE(.true.)
                 IF (global .ge. wland(i)%landp0) THEN
                        rank = i
                 ELSE
                        EXIT
                 END IF
         CASE(.false.)
                 IF (global .ge. wland(i)%patch0) THEN
                        rank = i
                 ELSE
                        EXIT
                 END IF
         END SELECT
      END DO

      IF (point) THEN
              local = global - wland(rank)%landp0 + 1
      ELSE
              local = global - wland(rank)%patch0 + 1
      END IF

      RETURN

      END SUBROUTINE global2local

      SUBROUTINE local2global (wnp, wland, point, global, rank, local)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: wnp
      TYPE(lpdecomp_t), POINTER, DIMENSION(:), INTENT(IN) :: wland
      LOGICAL, INTENT(IN) :: point
      INTEGER, INTENT(OUT) :: global
      INTEGER, INTENT(IN) :: rank, local

      IF (point) THEN
              global = wland(rank)%landp0 + local - 1
      ELSE
              global = wland(rank)%patch0 + local - 1
      END IF

      RETURN

      END SUBROUTINE local2global

      END MODULE decomp

      PROGRAM qdecomp

      USE decomp

      IMPLICIT NONE

      CHARACTER(LEN=99) :: decompfile
      CHARACTER(LEN=8) :: grid
      CHARACTER(LEN=10) :: argstr

      TYPE(lpdecomp_t), POINTER, DIMENSION(:) :: wland

      INTEGER :: wnp, nargs, global, rank, local

      LOGICAL :: g2l, point

      nargs = iargc()

      IF (nargs .lt. 3 .or. nargs .gt. 4) THEN
          PRINT *, &
          'Usage: qdecomp file [point|patch] [global|rank local]'
          STOP
      END IF

      CALL getarg(1, decompfile)
      CALL getarg(2, grid)

      IF (grid == 'point') THEN
              point = .true.
      ELSE IF (grid == 'patch') THEN
              point = .false.
      ELSE
          PRINT *, &
          'Usage: qdecomp file [point|patch] [global|rank local]'
          STOP
      END IF

      IF (nargs .eq. 3) THEN
              g2l = .true.
              CALL getarg(3, argstr)
              READ (argstr,*) global
      ELSE
              g2l = .false.
              CALL getarg(3, argstr)
              READ (argstr,*) rank
              CALL getarg(4, argstr)
              READ (argstr,*) local
      END IF

      CALL load_info (decompfile, wnp, wland)

      IF (g2l) THEN
              CALL global2local (wnp, wland, point, global, rank, local)
              PRINT *,&
              'global ',grid,global,' is local ',local,' on rank ',rank
              
      ELSE
              CALL local2global (wnp, wland, point, global, rank, local)
              PRINT *,&
              'local ',grid,local,' on rank ',rank,' is global ',global
      END IF

      DEALLOCATE (wland)

      END PROGRAM qdecomp

