!!this code illustrates the dynamics of using arrays in another function which have been declared in a module outside of the function  - pointers need to be used

module func_mod 
   use var_mod 
   implicit none
      contains
      subroutine sr1
         implicit none
           call falloc
            ar = 2.1
            print *,ar
         return
      end subroutine sr1
end module func_mod


