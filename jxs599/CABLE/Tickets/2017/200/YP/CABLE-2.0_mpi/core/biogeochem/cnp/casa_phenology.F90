SUBROUTINE phenology(iday,veg,phen)
  IMPLICIT NONE
  INTEGER,              INTENT(IN)    :: iday
  TYPE (veg_parameter_type),    INTENT(INOUT) :: veg  ! vegetation parameters
  TYPE (phen_variable), INTENT(INOUT) :: phen

  ! local variables (temprary)
  INTEGER :: np
  INTEGER, DIMENSION(mp)  :: days,days1to2, days2to3, days3to4, days4to1 

!  PRINT *, 'Within SUBROUTINE phenology, mp = ', mp 
  DO np=1,mp
    days1to2(np) = phen%doyphase(np,2) - phen%doyphase(np,1)
    days2to3(np) = phen%doyphase(np,3) - phen%doyphase(np,2)
    days3to4(np) = phen%doyphase(np,4) - phen%doyphase(np,3)
    days4to1(np) = phen%doyphase(np,1) - phen%doyphase(np,4)
    IF(days1to2(np) < 0) days1to2(np) = days1to2(np) +365
    IF(days2to3(np) < 0) days2to3(np) = days2to3(np) +365
    IF(days3to4(np) < 0) days3to4(np) = days3to4(np) +365
    IF(days4to1(np) < 0) days4to1(np) = days4to1(np) +365 
  ENDDO
  ! compute leaf phenology
  DO np=1,mp
    SELECT CASE(phen%phase(np))
      CASE(0)
        days(np) = iday - phen%doyphase(np,4)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days4to1(np)) phen%phase(np) =1
      CASE(1)
        days(np) = iday - phen%doyphase(np,1)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days1to2(np)) phen%phase(np) =2
      CASE(2)
        days(np) = iday - phen%doyphase(np,2)
        IF(days(np) <0) days(np) = days(np) + 365
        IF(days(np) > days2to3(np)) phen%phase(np) =3
      CASE(3)
        days(np) = iday - phen%doyphase(np,3)
        IF(days(np) < 0) days(np) = days(np) +365
        IF(days(np) > days3to4(np)) phen%phase(np) =0
    END SELECT
  ENDDO

  WHERE(veg%iveg==1 .or. veg%iveg ==2 )
       phen%phase = 2
  ENDWHERE

END SUBROUTINE phenology


