module cable_iday_mod
   implicit none
   public
   integer, pointer, save :: iday_number   

contains

subroutine iday_kick(iday)
   integer,target :: iday
   iday_number => iday  
end subroutine iday_kick

end module cable_iday_mod
