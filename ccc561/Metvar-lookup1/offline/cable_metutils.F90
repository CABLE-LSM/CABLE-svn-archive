!==============================================================================
! This source code is part of the
! Australian Community Atmosphere Biosphere Land Exchange (CABLE) model.
! This work is licensed under the CSIRO Open Source Software License
! Agreement (variation of the BSD / MIT License).
!
! You may not use this file except in compliance with this License.
! A copy of the License (CSIRO_BSD_MIT_License_v2.0_CABLE.txt) is located
! in each directory containing CABLE code.
!
! ==============================================================================
! Purpose: lookup tables for possible names for the met. forcing variables in
!          netcdf files
!
! Contact: jhan.srbinovsky@csiro.au
!
! History: First version by Claire Carouge (09/2021)
!
! ==============================================================================
!
module Cable_MetUtils
  USE netcdf
  
  implicit none
    
  PUBLIC :: metrain_names

  ! Rain variable
  CHARACTER(LEN=8),DIMENSION(5) :: metrain_names=(/                          &
                        'prcp    ', 'pr      ', 'RAIN    ', 'Rainf   ',        &
                        'Precip  '                                             &
                        /)
    
CONTAINS
  !==============================================================================
  !
  ! Name: find_metvarid
  !
  ! Purpose: Find the var ID in the met. file. It will use a lookup table of 
  !   possible variable names to find the variable in the netcdf file. 
  !
  ! CALLed from: open_met_file
  !
  !==============================================================================
  SUBROUTINE find_metvarid(file_id, possible_names, varid, status)
    
    IMPLICIT NONE

    ! Arguments
    ! ID of the netcdf with the met. forcing
    INTEGER, INTENT(IN) :: file_id  
    ! Lookup table for names of the variable in netcdf files
    CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: possible_names

    ! Output: variable id in netcdf file
    INTEGER, INTENT(OUT) :: varid 
    ! Output: error status if variable not found.
    INTEGER, INTENT(OUT) :: status
    
    ! Local variables
    INTEGER :: i  ! Loop index

    ! Initialise varid
    varid = -10
    DO i = 1, size(possible_names)
      status = NF90_INQ_VARID(file_id,trim(possible_names(i)),varid)
      ! Stop the loop if variable found
      IF (status == NF90_NOERR) exit
    END DO

  END SUBROUTINE find_metvarid

end module Cable_MetUtils