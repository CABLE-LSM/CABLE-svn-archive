
module test_module 
   use generic_pars 
   implicit none
   private
   public :: test_sr
      interface test_sr
         module procedure test_sr1, test_sr2
      end interface test_sr

      contains

      subroutine test_sr1
         implicit none
         return
      end subroutine test_sr1

      subroutine test_sr2
         implicit none
         return
      end subroutine test_sr2

end module test_module 


program test_driver
   use test_module 
   implicit none
      call test_sr
   stop
end program test_driver

