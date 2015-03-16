

MODULE cable_farray_module
   IMPLICIT NONE
   INTEGER, PARAMETER :: gok=0
   INTEGER :: galloctest=1
 
   interface cable_farray 
      module procedure cable_farray1, cable_farray2
   end interface cable_farray 
   
   integer, parameter ::                     & 
      farray_nmax = 50
    
   character(len=30), dimension(farray_nmax) :: &
      farray_names  
       
   real, dimension(:,:), allocatable :: & 
      farray_fields

   real, dimension(:,:,:), allocatable :: & 
      farray_fields2

CONTAINS
 
!==========================================================================!

SUBROUTINE cable_farray1( mp, CheckNames, CheckFields, &
                  n1,f1, n2,f2, n3,f3, n4,f4, n5,f5, n6,f6, n7,f7, & 
                  n8,f8, n9,f9, n10,f10, n11,f11, n12,f12, n13,f13, &
                  n14,f14, n15,f15, n16,f16, n17,f17, n18,f18, & 
                  n19,f19, n20,f20, n21,f21, n22,f22, n23,f23, & 
                  n24,f24, n25,f25, n26,f26, n27,f27, n28,f28, & 
                  n29,f29, n30,f30, n31,f31, n32,f32, n33,f33, & 
                  n34,f34, n35,f35, n36,f36, n37,f37, n38,f38, & 
                  n39,f39, n40,f40, n41,f41, n42,f42, n43,f43, & 
                  n44,f44, n45,f45, n46,f46, n47,f47, n48,f48, & 
                  n49,f49, n50,f50 & 
               )

   integer :: mp
    
   character(len=*), optional :: & 
                        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, & 
                        n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
                        n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
                        n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
                        n41, n42, n43, n44, n45, n46, n47, n48, n49, n50  

   real, dimension(:), optional :: & 
                        f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, & 
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,& 
                        f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,& 
                        f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,& 
                        f41, f42, f43, f44, f45, f46, f47, f48, f49, f50
   
   character(len=30), dimension(:), allocatable :: CheckNames 
   real, dimension(:,:), allocatable :: CheckFields

   integer :: i, k

   allocate( farray_fields( farray_nmax, mp) ) 
   i = 0

   if( present (n1) .AND. present(f1) ) then
      CALL fill_farray( n1, f1, i )   
   else
      print *, "CABLE_log: cable_farray missing dummy args"
      return
   endif
   
   if( present (n2) .AND. present(f2) ) CALL fill_farray( n2, f2, i )   
   if( present (n3) .AND. present(f3) ) CALL fill_farray( n3, f3, i )   
   if( present (n4) .AND. present(f4) ) CALL fill_farray( n4, f4, i )   
   if( present (n5) .AND. present(f5) ) CALL fill_farray( n5, f5, i )   
   if( present (n6) .AND. present(f6) ) CALL fill_farray( n6, f6, i )   
   if( present (n7) .AND. present(f7) ) CALL fill_farray( n7, f7, i )   
   if( present (n8) .AND. present(f8) ) CALL fill_farray( n8, f8, i )   
   if( present (n9) .AND. present(f9) ) CALL fill_farray( n9, f9, i )   
   if( present (n10) .AND. present(f10) ) CALL fill_farray( n10, f10, i )   
   if( present (n11) .AND. present(f11) ) CALL fill_farray( n11, f11, i )   
   if( present (n12) .AND. present(f12) ) CALL fill_farray( n12, f12, i )   
   if( present (n13) .AND. present(f13) ) CALL fill_farray( n13, f13, i )   
   if( present (n14) .AND. present(f14) ) CALL fill_farray( n14, f14, i )   
   if( present (n15) .AND. present(f15) ) CALL fill_farray( n15, f15, i )   
   if( present (n16) .AND. present(f16) ) CALL fill_farray( n16, f16, i )   
   if( present (n17) .AND. present(f17) ) CALL fill_farray( n17, f17, i )   
   if( present (n18) .AND. present(f18) ) CALL fill_farray( n18, f18, i )   
   if( present (n19) .AND. present(f19) ) CALL fill_farray( n19, f19, i )   
   if( present (n20) .AND. present(f20) ) CALL fill_farray( n20, f20, i )   
   if( present (n21) .AND. present(f21) ) CALL fill_farray( n21, f21, i )   
   if( present (n22) .AND. present(f22) ) CALL fill_farray( n22, f22, i )   
   if( present (n23) .AND. present(f23) ) CALL fill_farray( n23, f23, i )   
   if( present (n24) .AND. present(f24) ) CALL fill_farray( n24, f24, i )   
   if( present (n25) .AND. present(f25) ) CALL fill_farray( n25, f25, i )   
   if( present (n26) .AND. present(f26) ) CALL fill_farray( n26, f26, i )   
   if( present (n27) .AND. present(f27) ) CALL fill_farray( n27, f27, i )   
   if( present (n29) .AND. present(f29) ) CALL fill_farray( n29, f29, i )   
   if( present (n30) .AND. present(f30) ) CALL fill_farray( n30, f30, i )   
   if( present (n31) .AND. present(f31) ) CALL fill_farray( n31, f31, i )   
   if( present (n32) .AND. present(f32) ) CALL fill_farray( n32, f32, i )   
   if( present (n33) .AND. present(f33) ) CALL fill_farray( n33, f33, i )   
   if( present (n34) .AND. present(f34) ) CALL fill_farray( n34, f34, i )   
   if( present (n35) .AND. present(f35) ) CALL fill_farray( n35, f35, i )   
   if( present (n36) .AND. present(f36) ) CALL fill_farray( n36, f36, i )   
   if( present (n37) .AND. present(f37) ) CALL fill_farray( n37, f37, i )   
   if( present (n38) .AND. present(f38) ) CALL fill_farray( n38, f38, i )   
   if( present (n39) .AND. present(f39) ) CALL fill_farray( n39, f39, i )   
   if( present (n40) .AND. present(f40) ) CALL fill_farray( n40, f40, i )   
   if( present (n41) .AND. present(f41) ) CALL fill_farray( n41, f41, i )   
   if( present (n42) .AND. present(f42) ) CALL fill_farray( n42, f42, i )   
   if( present (n50) .AND. present(f50) ) CALL fill_farray( n50, f50, i )   
   
   allocate( CheckNames(i) )
   allocate( CheckFields(i,mp) )

   do k=1,i
      CheckNames(k) = farray_names(k) 
      CheckFields(k,:) = farray_fields(k,:)
   enddo        

   deallocate( farray_fields )
    
END SUBROUTINE cable_farray1



SUBROUTINE cable_farray2( mp, np, CheckNames, CheckFields, &
                  n1,f1, n2,f2, n3,f3, n4,f4, n5,f5, n6,f6, n7,f7, & 
                  n8,f8, n9,f9, n10,f10, n11,f11, n12,f12, n13,f13, &
                  n14,f14, n15,f15, n16,f16, n17,f17, n18,f18, & 
                  n19,f19, n20,f20, n21,f21, n22,f22, n23,f23, & 
                  n24,f24, n25,f25, n26,f26, n27,f27, n28,f28, & 
                  n29,f29, n30,f30, n31,f31, n32,f32, n33,f33, & 
                  n34,f34, n35,f35, n36,f36, n37,f37, n38,f38, & 
                  n39,f39, n40,f40, n41,f41, n42,f42, n43,f43, & 
                  n44,f44, n45,f45, n46,f46, n47,f47, n48,f48, & 
                  n49,f49, n50,f50 & 
               )

   integer :: mp, np
    
   character(len=*), optional :: & 
                        n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, & 
                        n11, n12, n13, n14, n15, n16, n17, n18, n19, n20, &
                        n21, n22, n23, n24, n25, n26, n27, n28, n29, n30, &
                        n31, n32, n33, n34, n35, n36, n37, n38, n39, n40, &
                        n41, n42, n43, n44, n45, n46, n47, n48, n49, n50  

   real, dimension(:,:), optional :: & 
                        f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, & 
                        f11, f12, f13, f14, f15, f16, f17, f18, f19, f20,& 
                        f21, f22, f23, f24, f25, f26, f27, f28, f29, f30,& 
                        f31, f32, f33, f34, f35, f36, f37, f38, f39, f40,& 
                        f41, f42, f43, f44, f45, f46, f47, f48, f49, f50
   
   character(len=30), dimension(:), allocatable :: CheckNames 
   real, dimension(:,:,:), allocatable :: CheckFields

   integer :: i, k

   allocate( farray_fields2( farray_nmax, mp, np) ) 
   i = 0

   if( present (n1) .AND. present(f1) ) then
      CALL fill_farray2( n1, f1, i )   
   else
      print *, "CABLE_log: cable_farray missing dummy args"
      return
   endif
   
   if( present (n2) .AND. present(f2) ) CALL fill_farray2( n2, f2, i )   
   if( present (n3) .AND. present(f3) ) CALL fill_farray2( n3, f3, i )   
   if( present (n4) .AND. present(f4) ) CALL fill_farray2( n4, f4, i )   
   if( present (n5) .AND. present(f5) ) CALL fill_farray2( n5, f5, i )   
   if( present (n6) .AND. present(f6) ) CALL fill_farray2( n6, f6, i )   
   if( present (n7) .AND. present(f7) ) CALL fill_farray2( n7, f7, i )   
   if( present (n8) .AND. present(f8) ) CALL fill_farray2( n8, f8, i )   
   if( present (n9) .AND. present(f9) ) CALL fill_farray2( n9, f9, i )   
   if( present (n10) .AND. present(f10) ) CALL fill_farray2( n10, f10, i )   
   if( present (n11) .AND. present(f11) ) CALL fill_farray2( n11, f11, i )   
   if( present (n12) .AND. present(f12) ) CALL fill_farray2( n12, f12, i )   
   if( present (n13) .AND. present(f13) ) CALL fill_farray2( n13, f13, i )   
   if( present (n14) .AND. present(f14) ) CALL fill_farray2( n14, f14, i )   
   if( present (n15) .AND. present(f15) ) CALL fill_farray2( n15, f15, i )   
   if( present (n16) .AND. present(f16) ) CALL fill_farray2( n16, f16, i )   
   if( present (n17) .AND. present(f17) ) CALL fill_farray2( n17, f17, i )   
   if( present (n18) .AND. present(f18) ) CALL fill_farray2( n18, f18, i )   
   if( present (n19) .AND. present(f19) ) CALL fill_farray2( n19, f19, i )   
   if( present (n20) .AND. present(f20) ) CALL fill_farray2( n20, f20, i )   
   if( present (n21) .AND. present(f21) ) CALL fill_farray2( n21, f21, i )   
   if( present (n22) .AND. present(f22) ) CALL fill_farray2( n22, f22, i )   
   if( present (n23) .AND. present(f23) ) CALL fill_farray2( n23, f23, i )   
   if( present (n24) .AND. present(f24) ) CALL fill_farray2( n24, f24, i )   
   if( present (n25) .AND. present(f25) ) CALL fill_farray2( n25, f25, i )   
   if( present (n26) .AND. present(f26) ) CALL fill_farray2( n26, f26, i )   
   if( present (n27) .AND. present(f27) ) CALL fill_farray2( n27, f27, i )   
   if( present (n29) .AND. present(f29) ) CALL fill_farray2( n29, f29, i )   
   if( present (n30) .AND. present(f30) ) CALL fill_farray2( n30, f30, i )   
   if( present (n31) .AND. present(f31) ) CALL fill_farray2( n31, f31, i )   
   if( present (n32) .AND. present(f32) ) CALL fill_farray2( n32, f32, i )   
   if( present (n33) .AND. present(f33) ) CALL fill_farray2( n33, f33, i )   
   if( present (n34) .AND. present(f34) ) CALL fill_farray2( n34, f34, i )   
   if( present (n35) .AND. present(f35) ) CALL fill_farray2( n35, f35, i )   
   if( present (n36) .AND. present(f36) ) CALL fill_farray2( n36, f36, i )   
   if( present (n37) .AND. present(f37) ) CALL fill_farray2( n37, f37, i )   
   if( present (n38) .AND. present(f38) ) CALL fill_farray2( n38, f38, i )   
   if( present (n39) .AND. present(f39) ) CALL fill_farray2( n39, f39, i )   
   if( present (n40) .AND. present(f40) ) CALL fill_farray2( n40, f40, i )   
   if( present (n41) .AND. present(f41) ) CALL fill_farray2( n41, f41, i )   
   if( present (n42) .AND. present(f42) ) CALL fill_farray2( n42, f42, i )   
   if( present (n50) .AND. present(f50) ) CALL fill_farray2( n50, f50, i )   
   
   allocate( CheckNames(i) )
   allocate( CheckFields(i,mp, np) )

   do k=1,i
      CheckNames(k) = farray_names(k) 
      CheckFields(k,:,:) = farray_fields2(k,:,:)
   enddo        

   deallocate( farray_fields2 )
    
END SUBROUTINE cable_farray2

subroutine fill_farray( n, f, i )   
   
   character(len=*) :: n
   real, dimension(:) :: f
   integer :: i
      
      i=i+1
      farray_names(i) = n
      farray_fields(i,:) = f
   
end subroutine fill_farray 


subroutine fill_farray2( n, f, i )   
   
   character(len=*) :: n
   real, dimension(:,:) :: f
   integer :: i
      
      i=i+1
      farray_names(i) = n
      farray_fields2(i,:,:) = f
   
end subroutine fill_farray2 


!==========================================================================!
 
END MODULE cable_farray_module
