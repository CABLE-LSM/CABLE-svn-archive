MODULE other_constants
  USE define_dimensions, ONLY : i_d, r_1, nrb
  REAL, PARAMETER, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/) ! Gaussian integ. weights
!jhan:make spatially expl in Mk3L
  REAL, PARAMETER, DIMENSION(nrb) :: refl = (/ 0.07, 0.425, 0.00 /) ! YP nov2009
  REAL, PARAMETER, DIMENSION(nrb) :: taul = (/ 0.07, 0.425, 0.00/)  ! leaf transmittance
   !--- jhan: can make these trigger of #defines/namelist
   real, parameter :: RAD_THRESH = 0.01 
   real, parameter :: LAI_THRESH = 0.01 
END MODULE other_constants
