SUBROUTINE avgsoil(veg,soil,casamet)
! Get avg soil moisture, avg soil temperature
! need to estimate the land cell mean soil temperature and moisture weighted by the area fraction
! of each tile within the land cell

  IMPLICIT NONE
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (soil_parameter_type),   INTENT(INOUT) :: soil ! soil parameters
  TYPE (casa_met),              INTENT(INOUT) :: casamet

  INTEGER                     :: ns,nland
  logical :: Ticket146 = .false.
  casamet%tsoilavg   = 0.0
  casamet%moistavg   = 0.0
  casamet%btran      = 0.0

  DO ns = 1, ms
  DO nland=1,mp
    casamet%tsoilavg(nland)  = casamet%tsoilavg(nland)+veg%froot(nland,ns)  &
                             * casamet%tsoil(nland,ns)
    casamet%moistavg(nland)  = casamet%moistavg(nland)+ veg%froot(nland,ns) &
                           * min(soil%sfc(nland),casamet%moist(nland,ns))
    casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
            * (min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland)) &
            /(soil%sfc(nland)-soil%swilt(nland))

  if( .NOT. Ticket146) &
    ! Ticket#121
    casamet%btran(nland)     = casamet%btran(nland)+ veg%froot(nland,ns)  &
            * (max(min(soil%sfc(nland),casamet%moist(nland,ns))-soil%swilt(nland),0.0)) &
            /(soil%sfc(nland)-soil%swilt(nland))

  ENDDO
  ENDDO

END SUBROUTINE avgsoil

