module mic_variable
  use mic_constant
  IMPLICIT NONE
  SAVE

  TYPE mic_parameter
  real(r_2), dimension(:),      allocatable  :: xav,xak,xdesorp,xbeta,xdiffsoc,rootbetax
  real(r_2), dimension(:,:),    allocatable  :: K1,K2,K3,J1,J2,J3
  real(r_2), dimension(:,:),    allocatable  :: V1,V2,V3,W1,W2,W3
  real(r_2), dimension(:,:),    allocatable  :: desorp
  real(r_2), dimension(:,:),    allocatable  :: Q1,Q2,fm,fs
  real(r_2), dimension(:,:),    allocatable  :: mgeR1,mgeR2,mgeR3,mgeK1,mgeK2,mgeK3
  real(r_2), dimension(:,:),    allocatable  :: tvmicR,tvmicK,betamicR,betamicK
  real(r_2), dimension(:,:),    allocatable  :: fmetave
  real(r_2), dimension(:,:,:),  allocatable  :: cn_r
  real(r_2), dimension(:,:),    allocatable  :: fr2p,fk2p,fr2c,fk2c,fr2a,fk2a
  real(r_2), dimension(:),      allocatable  :: xcnleaf,xcnroot,xcnwood,fligleaf,fligroot,fligwood
  real(r_2), dimension(:),      allocatable  :: diffsocx
  ! the following are alrealy available in CABLE
  integer,   dimension(:),      allocatable  :: pft,region,siteid
  real(r_2), dimension(:,:),    allocatable  :: sdepth,fracroot
  real(r_2), dimension(:,:),    allocatable  :: clay
  real(r_2), dimension(:,:),    allocatable  :: csoilobs,bulkd  
  END TYPE mic_parameter
  
  TYPE mic_input
  real(r_2), dimension(:,:),    allocatable  :: tavg,wavg
  real(r_2), dimension(:),      allocatable  :: dleaf,dwood,droot
  real(r_2), dimension(:,:),    allocatable  :: cinputm
  real(r_2), dimension(:,:),    allocatable  :: cinputs
  real(r_2), dimensioN(:),      allocatable  :: fcnpp

  END TYPE mic_input
 
  TYPE mic_output
  real(r_2), dimension(:,:),    allocatable  :: rsoil   
  END TYPE mic_output
  
  TYPE mic_cpool
  real(r_2), dimension(:,:,:),  allocatable  :: cpool
  END TYPE mic_cpool
 
  TYPE mic_npool
  real(r_2), dimension(:,:),    allocatable  :: mineralN
  END TYPE mic_npool 
  
 
 CONTAINS

  SUBROUTINE mic_allocate_parameter(mp,ms,micparam)
   IMPLICIT NONE
   TYPE(mic_parameter), INTENT(INOUT)  :: micparam
   integer  mp,ms

    ! tunable model parameters
    allocate(micparam%xav(mp),      &
	         micparam%xak(mp),      &
             micparam%xdesorp(mp),  &
             micparam%xbeta(mp),    &
             micparam%xdiffsoc(mp), &
             micparam%rootbetax(mp))
    allocate(micparam%K1(mp,ms),  &
             micparam%K2(mp,ms),  & 
             micparam%K3(mp,ms),  & 
             micparam%J1(mp,ms),  & 
             micparam%J2(mp,ms),  & 
             micparam%J3(mp,ms),  & 
             micparam%V1(mp,ms),  & 
             micparam%V2(mp,ms),  & 
             micparam%V3(mp,ms),  & 
             micparam%W1(mp,ms),  & 
             micparam%W2(mp,ms),  & 
             micparam%W3(mp,ms),  & 
             micparam%desorp(mp,ms), &
             micparam%Q1(mp,ms),     &
             micparam%Q2(mp,ms),     &
             micparam%fm(mp,ms),     &
             micparam%fs(mp,ms),     &
             micparam%mgeR1(mp,ms),  & 
             micparam%mgeR2(mp,ms),  & 
             micparam%mgeR3(mp,ms),  & 
             micparam%mgeK1(mp,ms),  & 
             micparam%mgeK2(mp,ms),  & 
             micparam%mgeK3(mp,ms),  & 
             micparam%fmetave(mp,ms),&
             micparam%tvmicR(mp,ms), &
             micparam%tvmicK(mp,ms), &
             micparam%betamicR(mp,ms),     &
             micparam%betamicK(mp,ms),     &
             micparam%cn_r(mp,ms,mcpool),  &
             micparam%fr2p(mp,ms),   & 
             micparam%fk2p(mp,ms),   & 
             micparam%fr2c(mp,ms),   & 
             micparam%fk2c(mp,ms),   &
             micparam%fr2a(mp,ms),   & 
             micparam%fk2a(mp,ms))

    allocate(micparam%xcnleaf(mp),   &
             micparam%xcnroot(mp),   &
             micparam%xcnwood(mp),   &
             micparam%fligleaf(mp),  &
             micparam%fligroot(mp),  &
             micparam%fligwood(mp),  &
             micparam%diffsocx(mp))

    allocate(micparam%pft(mp),       &
             micparam%region(mp),    &
             micparam%siteid(mp))

    allocate(micparam%sdepth(mp,ms),   &
             micparam%fracroot(mp,ms), &
             micparam%clay(mp,ms),     &
             micparam%csoilobs(mp,ms), &
             micparam%bulkd(mp,ms)) 
   
  END SUBROUTINE mic_allocate_parameter
  
  SUBROUTINE mic_allocate_input(mp,ms,micinput)
   IMPLICIT NONE
   integer mp,ms
   TYPE(mic_input), INTENT(INOUT)  :: micinput

    allocate(micinput%tavg(mp,ms),    &
             micinput%wavg(mp,ms),    &
             micinput%fcnpp(mp),      &
             micinput%dleaf(mp),      &
             micinput%dwood(mp),      &
             micinput%droot(mp),      &
             micinput%cinputm(mp,ms), &
             micinput%cinputs(mp,ms) )
   
  END SUBROUTINE mic_allocate_input
  
  SUBROUTINE mic_allocate_output(mp,ms,micoutput)
   IMPLICIT NONE
   TYPE(mic_output), INTENT(INOUT)  :: micoutput
   integer  mp,ms
   
   allocate(micoutput%rsoil(mp,ms))

  END SUBROUTINE mic_allocate_output  
  
  SUBROUTINE mic_allocate_cpool(mp,ms,miccpool)
   IMPLICIT NONE
   integer mp,ms
   TYPE(mic_cpool), INTENT(INOUT)  :: miccpool

   allocate(miccpool%cpool(mp,ms,mcpool))
   
  END SUBROUTINE mic_allocate_cpool 

  SUBROUTINE mic_allocate_npool(mp,ms,micnpool)
   IMPLICIT NONE
   integer mp,ms
   TYPE(mic_npool), INTENT(INOUT)  :: micnpool

   ALLOCATE(micnpool%mineralN(mp,ms))
   
  END SUBROUTINE mic_allocate_npool 
  
end module mic_variable