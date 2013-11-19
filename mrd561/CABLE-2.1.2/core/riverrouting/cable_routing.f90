module cable_routing

  use routing_constants, only :: rrlat,rrlon,rrlatlon
  use cable_types

  implicit none



contains

 
!-----------------------------------------------------------------------------
  integer function dirc2latindex(in_dirc)
    implicit none
    integer, intent(in) :: in_dirc

    dirc2latindex = 0
    if ((in_dirc .le. 1) .or.  (in_dirc .eq. 9)) dirc2latindex  = rrlat
    if ((in_dirc .ge. 3) .and. (in_dirc .le. 5)) dirc2latindex = -rrlat
    dirc2latindex = max(1,dirc2latindex)
    dirc2latindex = min(dirc2latindex,rrlatlon)

  end function dirc2latindex
!-----------------------------------------------------------------------------
  

!-----------------------------------------------------------------------------
  integer function dirc2lonindex(in_dirc)
    implicit none
    integer, intent(in) :: in_dirc

    dirc2lonindex = 0
    if ((in_dirc .ge. 1) .and. (in_dirc .le. 3)) dirc2lonindex  = 1
    if ((in_dirc .ge. 5) .and. (in_dirc .le. 7)) dirc2lonindex = -1
    dirc2lonindex = max(1,dirc2lonindex)
    dirc2lonindex = min(dirc2lonindex,rrlatlon) 

  end function dirc2lonindex
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  subroutine find_downstream_index(dstrm_index,river_dirc)

    integer, dimension(:), intent(out) :: dstrm_index
    integer, dimension(:), intent(in ) :: river_dirc
    
    integer :: i,j,n,ir,jr
    dwnstrm_index(:) = 0
    do j=1,rrlat
      do i=1,rrlon
        n = (j-1)*rrlon + i
        if (rdirc(n) /= 0) then
          ir = i + dirc2latindex(abs(river_dirc(n)))
          jr = j + dirc2lonindex(abs(river_dirc(n)))
          if (ir .lt. 1     ) ir = ir + rtmlon
          if (ir .gt. rtmlon) ir = ir - rtmlon
          if (jr .lt. 1 .or. jr .gt. rrlat .or. ir .lt. 1 .or. ir .gt. rrlon) then
            stop
          endif
          nr = (jr-1)*rrlon + ir
          if (n == nr) then
            stop
          endif
          dstrm_index(n) = nr
        endif
      enddo
    enddo

  end subroutine find_downstream_index
!-----------------------------------------------------------------------------




!-----------------------------------------------------------------------------
  subroutine find_stream_length(dstrm_index)
  end subroutine find_stream_length
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  subroutine map_cable2river()

  end subroutine map_cable2river
!-----------------------------------------------------------------------------
  
!-----------------------------------------------------------------------------
  subroutine step_riverflow()

  end subroutine step_riverflow
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
  subroutine map_river2cable()

  end subroutine map_river2cable
!-----------------------------------------------------------------------------

 
!-----------------------------------------------------------------------------
  subroutine mass_balance()

  end surbourinte mass_balance
!-----------------------------------------------------------------------------
o







end module cable_routing