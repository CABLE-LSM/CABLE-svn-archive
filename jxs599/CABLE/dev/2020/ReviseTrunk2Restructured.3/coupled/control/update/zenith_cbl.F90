module zenith_mod_cbl

REAL, ALLOCATABLE :: latitude(:,:) ! CABLE needs lat for zenith calc per dt   

contains

subroutine latitudeFix_cbl( latitude )
USE theta_field_sizes,        ONLY: t_i_length, t_j_length
USE model_grid_mod, ONLY : grid_lat => latitude
implicit none
real, allocatable :: latitude(:,:)
logical, save :: first_call = .TRUE.

  IF ( first_call) THEN 
    ALLOCATE( latitude(t_i_length,t_j_length) )
    write(6,*) "    "
    write(6,*) "To hyper-accurately test CABLE-JULES vs CABLE-CABLE we need to "
    write(6,*) "use the same zenith angle. It seems that JULES does not        "
    write(6,*) "remember *latitude* past the first timestep anyway, making it  "
    write(6,*) "less intrusive just to include CABLE's calc of zenith angle    "
    write(6,*) "    "
  
    latitude = grid_lat
    first_call = .FALSE.
  END IF

End subroutine latitudeFix_cbl

subroutine zenith_cbl( cos_zenith )
USE datetime_mod,       ONLY: l_360, l_leap
USE datetime_utils_mod, ONLY: day_of_year
USE model_time_mod,     ONLY: current_time
USE theta_field_sizes,  ONLY: row_length => t_i_length,   & 
                              rows       => t_j_length

USE cbl_sinbet_mod, ONLY : sinbet

implicit none
REAL, INTENT(OUT) :: cos_zenith(row_length, rows) ! Cosine of zenith angle
integer ::  curr_day_number
real :: current_hour
integer :: i, j

curr_day_number = day_of_year( current_time%year, current_time%month,        &
                               current_time%day, l_360, l_leap )
current_hour = (current_time%time/3600.) + .25

  ! Elemental function, in this case over (row_length,rows)
  do i=1,  row_length
    do j=1,  rows
      cos_zenith(i,j) = sinbet( real(curr_day_number),latitude(i,j),         &
                                real(current_hour) ) 
      write(6,*) "lat ij_cbl(i,j) ", latitude(i,j)
    end do
  end do

End subroutine zenith_cbl

Endmodule zenith_mod_cbl
