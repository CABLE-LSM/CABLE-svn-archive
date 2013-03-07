module arraydiag_m
   use ieee_arithmetic
   implicit none
   interface arraydiag
      module procedure arraydiag1, arraydiag2, arraydiag3, arraydiag4
   end interface

contains

   subroutine arraydiag1(str,a,unit,stats)
      character(len=*), intent(in) :: str
      real, dimension(:), intent(in) :: a
      integer, intent(in) :: unit
      logical, intent(in), optional :: stats
      integer :: i, nnan, nanloc
      logical, dimension(size(a)) :: amask
      real :: mean, var

      amask = ieee_is_finite(a)

      if ( present(stats)) then
         if (stats) then
            mean = sum(a,amask)/count(amask)
            var = sum((a-mean)**2,amask)/count(amask)
            write(unit,*) trim(str), " Mean:", mean, "SD:",  sqrt(var)
         end if
      end if
      write(unit,*) trim(str), " Count:", count(amask)
      write(unit,*) trim(str), " Range:", minval(a,amask), maxval(a,amask)
      write(unit,*) trim(str), " Minloc:", minloc(a,amask)
      write(unit,*) trim(str), " Maxloc:", maxloc(a,amask)

      if (any(ieee_is_nan(a))) then
         write(unit,*) trim(str), " Array has", count(ieee_is_nan(a)), "NaN values"
         write(unit,*) trim(str), " First NaN at", maxloc(merge(1,0,ieee_is_nan(a)))
      end if
      if (.not. all(amask)) then
         write(unit,*) trim(str), " Array has", count(.not. amask), "non-finite values"
         write(unit,*) trim(str), " First non-finite at", maxloc(merge(0,1,amask))
      end if

   end subroutine arraydiag1

   subroutine arraydiag2(str,a,unit,stats)
      character(len=*), intent(in) :: str
      real, dimension(:,:), intent(in) :: a
      integer, intent(in) :: unit
      logical, intent(in), optional :: stats
      integer :: i, nnan, nanloc
      logical, dimension(size(a,1),size(a,2)) :: amask
      real :: mean, var

      amask = ieee_is_finite(a)

      if ( present(stats)) then
         if (stats) then
            mean = sum(a,amask)/count(amask)
            var = sum((a-mean)**2,amask)/count(amask)
            write(unit,*) trim(str), " Mean:", mean, "SD:",  sqrt(var)
         end if
      end if
      write(unit,*) trim(str), " Range:", minval(a,amask), maxval(a,amask)
      write(unit,*) trim(str), " Minloc:", minloc(a,amask)
      write(unit,*) trim(str), " Maxloc:", maxloc(a,amask)

      if (any(ieee_is_nan(a))) then
         write(unit,*) trim(str), " Array has", count(ieee_is_nan(a)), "NaN values"
         write(unit,*) trim(str), " First NaN at", maxloc(merge(1,0,ieee_is_nan(a)))
      end if
      if (.not. all(amask)) then
         write(unit,*) trim(str), " Array has", count(.not. amask), "non-finite values"
         write(unit,*) trim(str), " First non-finite at", maxloc(merge(0,1,amask))
      end if

   end subroutine arraydiag2

   subroutine arraydiag3(str,a,unit,stats)
      character(len=*), intent(in) :: str
      real, dimension(:,:,:), intent(in) :: a
      integer, intent(in) :: unit
      logical, intent(in), optional :: stats
      integer :: i, nnan, nanloc
      logical, dimension(size(a,1),size(a,2),size(a,3)) :: amask
      real :: mean, var

      amask = ieee_is_finite(a)

      if ( present(stats)) then
         if (stats) then
            mean = sum(a,amask)/count(amask)
            var = sum((a-mean)**2,amask)/count(amask)
            write(unit,*) trim(str), " Mean:", mean, "SD:",  sqrt(var)
         end if
      end if
      write(unit,*) trim(str), " Range:", minval(a,amask), maxval(a,amask)
      write(unit,*) trim(str), " Minloc:", minloc(a,amask)
      write(unit,*) trim(str), " Maxloc:", maxloc(a,amask)

      if (any(ieee_is_nan(a))) then
         write(unit,*) trim(str), " Array has", count(ieee_is_nan(a)), "NaN values"
         write(unit,*) trim(str), " First NaN at", maxloc(merge(1,0,ieee_is_nan(a)))
      end if
      if (.not. all(amask)) then
         write(unit,*) trim(str), " Array has", count(.not. amask), "non-finite values"
         write(unit,*) trim(str), " First non-finite at", maxloc(merge(0,1,amask))
      end if


   end subroutine arraydiag3

   subroutine arraydiag4(str,a,unit,stats)
      character(len=*), intent(in) :: str
      real, dimension(:,:,:,:), intent(in) :: a
      integer, intent(in) :: unit
      logical, intent(in), optional :: stats
      integer :: i, nnan, nanloc
      logical, dimension(size(a,1),size(a,2),size(a,3),size(a,4)) :: amask
      real :: mean, var

      amask = ieee_is_finite(a)

      if ( present(stats)) then
         if (stats) then
            mean = sum(a,amask)/count(amask)
            var = sum((a-mean)**2,amask)/count(amask)
            write(unit,*) trim(str), " Mean:", mean, "SD:",  sqrt(var)
         end if
      end if
      write(unit,*) trim(str), " Range:", minval(a,amask), maxval(a,amask)
      write(unit,*) trim(str), " Minloc:", minloc(a,amask)
      write(unit,*) trim(str), " Maxloc:", maxloc(a,amask)

      if (any(ieee_is_nan(a))) then
         write(unit,*) trim(str), " Array has", count(ieee_is_nan(a)), "NaN values"
         write(unit,*) trim(str), " First NaN at", maxloc(merge(1,0,ieee_is_nan(a)))
      end if
      if (.not. all(amask)) then
         write(unit,*) trim(str), " Array has", count(.not. amask), "non-finite values"
         write(unit,*) trim(str), " First non-finite at", maxloc(merge(0,1,amask))
      end if


   end subroutine arraydiag4

end module arraydiag_m
