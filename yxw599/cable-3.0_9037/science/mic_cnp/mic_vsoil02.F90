module mic_constant
  IMPLICIT NONE
  integer,  parameter  :: r_2 = SELECTED_REAL_KIND(12, 60)
  integer,  parameter  :: diag=0       ! =1 for printout 0 no prinout
  integer,  parameter  :: outp=1       ! output site
  integer,  parameter  :: msite=182    ! number of sites
  integer,  parameter  :: mp=1000      ! number of site the model runs for
  integer,  parameter  :: ms=15        ! soil layers
  integer,  parameter  :: mcpool=7     ! number of C pools
  integer,  parameter  :: nfvar=36     ! number of data input variables
  real,     parameter  :: delt= 1.0    ! one hour
  integer,  parameter  :: mpft=15      ! number of PFTs
  real(r_2),PARAMETER  :: diffsoc  =(1.0/24.0)* 2.74e-3  !cm2/hour   
                                        ! m2/hour  ! see Table 1,Camino-Serrano et al. (2018)
  
end module mic_constant

module mic_variable
  use mic_constant
  IMPLICIT NONE
  SAVE

  TYPE mic_parameter
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
  real(r_2), dimension(:,:),    allocatable  :: tavg,wavg,tair
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
             micinput%tair(mp,365),   &
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

 real*8 function functn(nx,xparam16)
   use mic_constant
   use mic_variable
   implicit none
    TYPE(mic_parameter)  :: micparam
    TYPE(mic_input)      :: micinput
    TYPE(mic_cpool)      :: miccpool
    TYPE(mic_npool)      :: micnpool
    TYPE(mic_output)     :: micoutput

    !local variables
    real*8,    dimension(16)           :: xparam16
    integer    nx
    real*8,    dimension(nx)           :: xopt
    real*8     totcost
    integer isoc14,kinetics,pftopt,jopt,nyeqpool
    real(r_2), dimension(mp,ms,mcpool)  :: cpooleq
    real(r_2), dimension(ms)            :: zse  ! value set in "vmic_param_constant"
    data zse/0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1/  !in m, layer thickness 

      open(1,file='params1.txt')
      open(91,file='modobs.txt')
      open(92,file='modobs2.txt')
      read(1,*) 
      read(1,*) isoc14,kinetics,pftopt,jopt

      xopt =xparam16(1:nx)
      if(jopt==1) read(1,*) xopt(1:nx)

      totcost = 0.0
      nyeqpool= 434

      call mic_allocate_parameter(mp,ms,micparam)
      call mic_allocate_input(mp,ms,micinput)
      call mic_allocate_output(mp,ms,micoutput)
      call mic_allocate_cpool(mp,ms,miccpool)
      call mic_allocate_npool(mp,ms,micnpool)
  
      print *, 'getdat'  
      call getdata(micinput,micparam,micnpool)
  
      print *, 'vmicsoil'

      call vmicsoil(nx,xopt,kinetics,pftopt,nyeqpool, &
                micparam,micinput,miccpool,micnpool,micoutput,zse,cpooleq)

      call calcost(nx,isoc14,pftopt,xopt,micparam,zse,cpooleq,totcost) 
      print *, 'finished'
      close(1)
      close(91)
      close(92)
  
      functn = totcost

END function functn

  subroutine getdata(micinput,micparam,micnpool)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool
    real(r_2),  dimension(182,nfvar)     :: forcdata
    real(r_2),  dimension(212,365)       :: tairx
    !local variables
    INTEGER m,np,ns,nd,ndx,npx

      open(10,file='input/ICP_profile50.txt')          ! open input file
101   format(a80)
      read(10,*)                                 ! Variable name, jump the first line
      read(10,*)
      read(10,*)
      m = 0
      do np=1,50
         m = m + 1
         read(10,*,end=15) forcdata(m,1:nfvar)
      enddo
15    close(10)
      open(11,file='input/China_profile134.txt')
      read(11,*)                                ! Variable name, jump the first line
      read(11,*)
      read(11,*)
      do np=1,132
         m = m + 1
         read(11,*,end=25) forcdata(m,1:nfvar)
      enddo
25    close(11)

   
      micparam%csoilobs(:,:) = -999.0
      do np=1, mp
         npx = np
!         if(np>182) npx =1
   
         micparam%pft(np)    = int(forcdata(npx,3))
         micparam%region(np) = int(forcdata(npx,1))
         micparam%siteid(np) = int(forcdata(npx,2))

         do ns=1,ms
            micinput%tavg(np,ns)     = forcdata(npx,4)  ! average temperature in deg C
            micinput%wavg(np,ns)     = forcdata(npx,5)  ! average soil water content unit??
            micparam%clay(np,ns)     = forcdata(npx,6)  ! clay content 
            micnpool%mineralN(np,ns) = forcdata(npx,7)*0.001 ! mineral N: "0.001" mg N /kg soil --> g N /kg soil
          enddo !"ns"
          if(micparam%region(np)==1) then                  !ICP sites
             micparam%bulkd(np,1)       = forcdata(npx,28)  ! 0-10cm
             micparam%bulkd(np,2)       = forcdata(npx,30)  ! 10-20cm
             micparam%bulkd(np,3)       = forcdata(npx,32)  ! 20-40cm
             micparam%bulkd(np,4)       = forcdata(npx,34)  ! 40-80cm
             micparam%bulkd(np,5)       = micparam%bulkd(np,4)      ! 80-100cm

             micparam%csoilobs(np,1)    = forcdata(npx,27)
             micparam%csoilobs(np,2)    = forcdata(npx,29)
             micparam%csoilobs(np,3)    = forcdata(npx,31)  
             micparam%csoilobs(np,4)    = forcdata(npx,33)
             micparam%csoilobs(np,5)    = micparam%csoilobs(npx,4) 
          endif
    
          if(micparam%region(np)==3) then                 !China-C sites
            micparam%bulkd(np,1)       = 1000.0*forcdata(npx,27) ! 0-10cm ! original data have unit of Mg/m3
            micparam%bulkd(np,2)       = 1000.0*forcdata(npx,28) !10-20cm
            micparam%bulkd(np,3)       = 1000.0*forcdata(npx,29) !20-30cm
            micparam%bulkd(np,4)       = 1000.0*forcdata(npx,30) !30-50cm
            micparam%bulkd(np,5)       = 1000.0*forcdata(npx,31) !50-100cm

            micparam%csoilobs(np,1)    = forcdata(npx,32)
            micparam%csoilobs(np,2)    = forcdata(npx,33)
            micparam%csoilobs(np,3)    = forcdata(npx,34)
            micparam%csoilobs(np,4)    = forcdata(npx,35)
            micparam%csoilobs(np,5)    = forcdata(npx,36)
          endif
          ! make sure "*delt" is not repeated in the model called by rk4
          micinput%fcnpp(np)      = (forcdata(npx,14)+forcdata(npx,15)+forcdata(npx,16)+forcdata(npx,17))*delt
          micinput%Dleaf(np)      = (forcdata(npx,14)+forcdata(npx,15))/(24.0*365.0)*delt    !gc/m2/delt
          micinput%Droot(np)      = forcdata(npx,16)/(24.0*365.0)*delt     !gc/m2/delt
          micinput%Dwood(np)      = forcdata(npx,17)/(24.0*365.0)*delt     !gc/m2/delt

          micparam%xcnleaf(np)    = forcdata(npx,18)
          micparam%xcnroot(np)    = forcdata(npx,19)
          micparam%xcnwood(np)    = forcdata(npx,20)
          micparam%fligleaf(np)   = forcdata(npx,21)
          micparam%fligroot(np)   = forcdata(npx,22)
          micparam%fligwood(np)   = forcdata(npx,23)

      enddo    ! "np=1,mp"

      ! read in daily temperature data
      open(12,file='input/tair_C14_ICP_CN.txt')
      read(12,*)
      do nd=1,365
         read(12,*) ndx,(tairx(ns,nd),ns=1,212)
      enddo
      close(12)  
      do np=1,mp
         npx = np
!         if(np>182) npx = 1
         do nd=1,365
             micinput%tair(np,nd)=tairx(micparam%siteid(npx),nd)
         enddo
      enddo 

   end subroutine getdata
   
SUBROUTINE vmic_param_constant(nx,xopt,zse,micparam) 
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout) :: micparam
    integer nx
    real*8,    dimension(nx)           :: xopt
!local variables   
    real(r_2), dimension(mpft,ms)      :: froot
    real(r_2), dimension(mpft)         :: rootbeta
    data rootbeta/0.962,0.962,0.961,0.976,0.962,0.966,0.943,0.962,0.966,0.943,0.966,0.966,0.966,0.972,0.962/
    real(r_2), dimension(ms)           :: zse
    real(r_2), dimension(mpft)         :: totroot
    real(r_2) xdiffsoc,rootbetax
    integer ipft,np,ns
    real(r_2)  depths1,depths2
    real(r_2) Kor, Kok, Q1, Q2, fm, fs

      xdiffsoc=xopt(5); rootbetax=xopt(7)
  
      Kor = 4.0; Kok = 4.0; Q1= Kor; Q2  = Kok; fm = 0.05; fs= 0.05
      micparam%Q1(:,:) = Q1; micparam%Q2(:,:) =Q2; micparam%fm(:,:)=fm; micparam%fs(:,:)=fs
      
      depths1=0.0;depths2=0.0
      do ns=1,ms
         depths2 = depths2 + zse(ns)
         do ipft=1,mpft
             froot(ipft,ns) = (1.0/rootbetax) *( exp(-rootbetax*depths1)-exp(-rootbetax*depths2))
         enddo
         depths1=depths2
      enddo
      
      do ipft=1,mpft
         totroot(ipft) =sum(froot(ipft,1:ms))
      enddo

      !normalizing
      do ns=1,ms
         do ipft=1,mpft
            froot(ipft,ns) = froot(ipft,ns)/totroot(ipft)
         enddo
      enddo

      ! calculate mp by ms all parameter values
      do np=1, mp
         if(micparam%pft(np)<1 .or. micparam%pft(np)>15) print *, 'error: PFT incorrect', np, micparam%pft(np)
         do ns=1,ms
            micparam%sdepth(np,ns)   = zse(ns)
            micparam%fracroot(np,ns) = froot(micparam%pft(np),ns)
          enddo !"ns"
          micparam%diffsocx(np) = xdiffsoc * diffsoc  !"diffsoc" from mic_constant
      enddo    ! "np=1,mp"
  
      if(diag==1) then
         print *, micparam%fracroot(outp,:) 
         print *, micparam%sdepth(outp,:)
         print *, micparam%diffsocx(outp)
      endif
END SUBROUTINE vmic_param_constant  

subroutine vmic_param_time(nx,xopt,micparam,micinput,micnpool)
    ! time-dependent model parameters, called every time step if the forcing, such air temperature
    ! varies every time step
    ! otherwise only called at the start the integration	
    use mic_constant 
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    integer   nx
    real*8,   dimension(nx)              :: xopt
    ! local variables
    real(r_2)   xav, xak, xdesorp, xbeta

      xav = xopt(1); xak= xopt(2); xdesorp = xopt(3); xbeta=xopt(4)
      ! compute fractions
      call bgc_fractions(micparam,micinput)
      ! compute microbial growth efficiency
      call mget(micparam,micinput,micnpool)
      ! compute microbial turnover rates
      call turnovert(xbeta,micparam,micinput)
      call Desorpt(xdesorp,micparam) 
      call Vmaxt(xav,micparam,micinput)
      call Kmt(xak,micparam,micinput)

  
end subroutine vmic_param_time


subroutine vmic_init(nx,xopt,kinetics,pftopt,micparam,micinput,miccpool,micnpool)
    use mic_constant
    use mic_variable
    implicit none
    integer nx
    real*8, dimension(nx)                :: xopt 
    integer  kinetics,pftopt
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_cpool),     INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool

    ! local variables
    ! for numerical solution
    real(r_2),    parameter            :: tol = 1.0E-04
    real(r_2),    parameter            :: tolx = 0.0001
    real(r_2),    parameter            :: tolf = 0.000001
    integer,     parameter             :: ntrial = 100
    integer np,ns,ip
    real(r_2),    dimension(mcpool)     :: cpooldef,xpool0,y

!	  print *, 'calling vmic_init'

      cpooldef(1) = 16.5*0.1;     cpooldef(2) = 16.5*0.1
      cpooldef(3) = 16.5*0.025;   cpooldef(4) = 16.5*0.025
      cpooldef(5) = 16.5*0.1125;  cpooldef(6) = 16.5*0.375;  cpooldef(7) = 16.5*0.2625

      do ip=1,mcpool
         miccpool%cpool(:,:,ip) = cpooldef(ip)
      enddo
  
      call vmic_param_time(nx,xopt,micparam,micinput,micnpool)
  
      do np=1,mp
         if(micparam%pft(np)==pftopt) then
            do ns=1,ms

               ! initial pool sizes
               do ip=1,mcpool
                  xpool0(ip) = miccpool%cpool(np,ns,ip)
               enddo

               call mnewt(ntrial,np,ns,kinetics,micparam,micinput,xpool0,tolx,tolf)
               call vmic_c(np,ns,kinetics,micparam,micinput,xpool0,y)
               if(maxval(xpool0(1:mcpool))>1.0e4.or.minval(xpool0(1:mcpool))<0.0) then
                  xpool0 = cpooldef
               endif
               do ip=1,mcpool
                  miccpool%cpool(:,:,ip) = xpool0(ip)
               enddo
            enddo
         endif
       enddo
  

end subroutine vmic_init

subroutine vmicsoil(nx,xopt,kinetics,pftopt,nyeqpool, &
                    micparam,micinput,miccpool,micnpool,micoutput,zse,cpooleq)
    use mic_constant 
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(INout)   :: micparam
    TYPE(mic_input),     INTENT(INout)   :: micinput
    TYPE(mic_cpool),     INTENT(INOUT)   :: miccpool
    TYPE(mic_npool),     INTENT(INOUT)   :: micnpool
    TYPE(mic_output),    INTENT(INout)   :: micoutput

    integer nx,isoc14,kinetics,pftopt,nyeqpool
    real*8, dimension(nx)             :: xopt

    ! local variables
    ! for numerical solution
    real(r_2),    parameter            :: tol = 1.0E-04
    real(r_2),    parameter            :: tolx = 0.0001
    real(r_2),    parameter            :: tolf = 0.000001
    integer,      parameter            :: ntrial = 100

    ! variables exchange with vmic_parameters
    real(r_2) xav, xak, xdesorp,xavexp

    ! variables exchange with vmic_cost
    real(r_2), dimension(mp,ms,mcpool)  :: cpooleq
    real(r_2), dimension(ms)            :: zse  

    ! local variables
    real(r_2),    dimension(mcpool)    :: xpool0,xpool1
    real(r_2),    dimension(ms)        :: ypooli,ypoole,fluxsoc

    integer       ndelt,n1,n2,i,year,ip,np,ns,ny
    real(r_2)     timex,delty,fluxdocsx,diffsocxx

    integer j


      call vmic_param_constant(nx,xopt,zse,micparam) 
!	  call vmic_param_time(nx,xopt,micparam,micinput,micnpool)
!      call vmic_init(miccpool)
      call vmic_init(nx,xopt,kinetics,pftopt,micparam,micinput,miccpool,micnpool)

      ndelt   = int(24*365/delt) ! number of time step per year in "delt"

      do year=1,nyeqpool+66
	     print *, 'run for year and % to completion', year, real(year)/real(nyeqpool+66)
         ny = year-nyeqpool         
         do i=1,365    ! 365 days in a year

           ! using the daily air temperature 
           do np=1,mp
              micinput%tavg(np,1:ms) = micinput%tair(np,i)        ! using daily air temperature
           enddo
           ! only called if forcing varies every time step
           call vmic_param_time(nx,xopt,micparam,micinput,micnpool)

            do np=1,mp
            if(micparam%pft(np)==pftopt) then

               ! do MIMICS for each soil layer
               do ns=1,ms

                  do ip=1,mcpool
                     xpool0(ip) = miccpool%cpool(np,ns,ip)
                  enddo

                  timex=real(i*delt)
                  delty = real(ndelt)/(365.0*delt)  ! time step in rk4 in "delt"
                  call rk4modelx(timex,delty,np,ns,kinetics,micparam,micinput,xpool0,xpool1)              
                  ! the following used to avoid poolsize<0.0
                  do ip=1,mcpool
                     xpool0(ip) = max(1.0e-8,xpool1(ip))
                  enddo

                  xpool1=xpool0 
                  do ip=1,mcpool
                     miccpool%cpool(np,ns,ip) = xpool1(ip)
                  enddo
           
               enddo ! "ns"
  
               if(diag==1) then  
                  print *, 'year day site np1', year, i, outp,micparam%diffsocx(outp)
                  do ns=1,ms
                     print *, ns, miccpool%cpool(outp,ns,:) 
                  enddo  
               endif
  
               do ip=1,mcpool
                  do ns=1,ms
                     ypooli(ns) = miccpool%cpool(np,ns,ip)      ! in mg c/cm3
                  enddo  !"ns"

                  fluxsoc(:) = 0.0  ! This flux is added in "modelx"
                  diffsocxx= micparam%diffsocx(np)
  
                !  print *, 'before bioturb'
                !  print *, year,i,ip, int(delty),ms,zse(1),delt,diffsocxx,fluxsoc(1)
                !  print *, ip,ypooli(:)
  
                  ! do bioturbation for "ndelt" timesteps (to calculate equilirium pool only)
  
                  call bioturb(int(delty),ms,zse,delt,diffsocxx,fluxsoc,ypooli,ypoole)  ! only do every 24*delt

                !  print *, 'after bioturb'
                !  print *, ip, ypoole(:)

                  do ns=1,ms
                     miccpool%cpool(np,ns,ip) = ypoole(ns)
                  enddo
               enddo !"ip"

            endif   !pft(np) = pftopt
            enddo !"np"


            ! print *, 'fcost: ',year, i, micinput%tavg(1,1)
            ! do ns=1,ms
            !   print *, miccpool%cpool(1,ns,:)
            ! enddo

         enddo  ! i (day)  
      enddo !year
  
  
      cpooleq(:,:,:) = miccpool%cpool(:,:,:)
  
    end subroutine vmicsoil 

    subroutine calcost(nx,isoc14,pftopt,xopt,micparam,zse,cpooleq,totcost)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(IN)    :: micparam
    integer nx,isoc14,pftopt
    real*8  totcost
    real*8, dimension(nx)              :: xopt
    real(r_2), dimension(ms)           :: zse
    real(r_2), dimension(mp,ms,mcpool) :: cpooleq
    ! cost function
    real(r_2), dimension(mp)           :: xcost
    real(r_2), dimension(mp,ms)        :: xobs, xmod


    integer   np,ns,i,ip,msobs
    real(r_2)  xbdz


    ! units: carbon flux: mg C cm-3 delt-1

      ! write out the observed and modelled SOC (gC/kg soil) for each site
      xcost(1:mp) = 0.0; xobs(:,:)=0.0; xmod(:,:)=0.0
      do np=1,mp
         if(micparam%pft(np)==pftopt) then

            xbdz=0.0
            xobs(np,1) = micparam%csoilobs(np,1)
            xmod(np,1) = 1000.0 * sum(cpooleq(np,1,3:7))/micparam%bulkd(np,1)
            xobs(np,2) = micparam%csoilobs(np,2)
            xmod(np,2) = 1000.0 * sum(cpooleq(np,2,3:7))/micparam%bulkd(np,2)
            if(micparam%region(np) ==1) then ! ICP sites
               xobs(np,3) = micparam%csoilobs(np,3)
               xmod(np,3) = sum(cpooleq(np,3,3:7))*zse(3) + sum(cpooleq(np,4,3:7))*zse(4) 
               xbdz        = micparam%bulkd(np,3)*(zse(3)+zse(4))
               xmod(np,3) = 1000.0* xmod(np,3)/xbdz

               xobs(np,4) = micparam%csoilobs(np,4)
               xmod(np,4) = sum(cpooleq(np,5,3:7))*zse(5)+sum(cpooleq(np,6,3:7))*zse(6)   &
                          +sum(cpooleq(np,7,3:7))*zse(7)+sum(cpooleq(np,8,3:7))*zse(8)                         
               xbdz       = micparam%bulkd(np,4)*(zse(5)+zse(6)+zse(7)+zse(8))
               xmod(np,4) = 1000.0* xmod(np,4)/xbdz
            endif
            if(micparam%region(np) ==3) then ! China-C site
               xobs(np,3) = micparam%csoilobs(np,3)
               xmod(np,3) = 1000.0*sum(cpooleq(np,3,3:7))/micparam%bulkd(np,3)

               xobs(np,4) = micparam%csoilobs(np,4)
               xmod(np,4) = sum(cpooleq(np,4,3:7))*zse(4)+sum(cpooleq(np,5,3:7))*zse(5)
               xbdz       = micparam%bulkd(np,4)*(zse(4)+zse(5))
               xmod(np,4) = 1000.0* xmod(np,4)/xbdz

               xobs(np,5) = micparam%csoilobs(np,5)
               xmod(np,5) = sum(cpooleq(np,6,3:7))*zse(6)+sum(cpooleq(np,7,3:7))*zse(7)  &
                           +sum(cpooleq(np,8,3:7))*zse(8)+sum(cpooleq(np,9,3:7))*zse(9)  &
                           +sum(cpooleq(np,10,3:7))*zse(10)

               xbdz       = micparam%bulkd(np,5)*(zse(6)+zse(7)+zse(8)+zse(9)+zse(10))
               xmod(np,5) = 1000.0* xmod(np,5)/xbdz                                

            endif
            xcost(np) = 0.0
            msobs = 5

            do i=1,msobs
               if(xmod(np,i) <0.0.or.xmod(np,i) >1.0e3) then
                  print *, 'abnormal value of model simulation site=', np
                  print *, 'parameter values = ',  xopt(:)
                  print *,  xobs(np,i),xmod(np,i),xobs(np,i)-xmod(np,i),xcost(np)
                  do ns=1,ms
                     print *, ' modelled pool size= ', ns,cpooleq(np,ns,:)
                  enddo
                  xmod(np,i) = -1000.0
               endif

               if(xobs(np,i) > 0.0) then
                  xcost(np) = xcost(np) + (xobs(np,i)-xmod(np,i))**2
               endif
               write(91,901) micparam%siteid(np),micparam%pft(np),msobs,i,xobs(np,i),xmod(np,i),xcost(np)
            enddo

            do ns=1,ms
                write(92,921) micparam%siteid(np),micparam%pft(np), ns, (1000.0*cpooleq(np,ns,ip),ip=1,mcpool)
            enddo
        endif  !pft(np)=pftopt

      enddo  !"np"
      totcost = sum(xcost(1:mp))

901   format(4(i3,2x),5(f10.4,2x))
921   format(3(i3,2x),10(f12.4,1x))
    end subroutine calcost

    subroutine rk4modelx(timex,delty,np,ns,kinetics,micparam,micinput,xpool0,xpool1)
    use mic_constant
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(IN)  :: micparam
    TYPE(mic_input),     INTENT(IN)  :: micinput
    integer      np,ns, kinetics
    real(r_2)    timex,delty,h
    real(r_2),   dimension(mcpool),intent(inout)     :: xpool0,xpool1
    real(r_2),   dimension(mcpool)                   :: y1,y2,y3,y4,dy1dt,dy2dt,dy3dt,dy4dt

     h=delty
     y1(:) = xpool0(:)
 
     call vmic_c(np,ns,kinetics,micparam,micinput,y1,dy1dt)
     y2(:) = y1(:) + 0.5 * h * dy1dt(:)
     call vmic_c(np,ns,kinetics,micparam,micinput,y2,dy2dt)
     y3(:) = y1(:) + 0.5 * h * dy2dt(:)
     call vmic_c(np,ns,kinetics,micparam,micinput,y3,dy3dt)
     y4(:) = y1(:) +       h * dy3dt(:)
     call vmic_c(np,ns,kinetics,micparam,micinput,y4,dy4dt)

     xpool1(:) = xpool0(:) + (dy1dt(:)/6.0 + dy2dt(:)/3.0 + dy3dt(:)/3.0 + dy4dt(:)/6.0) * h

    end subroutine rk4modelx

    subroutine Kmt(xak,micparam,micinput)
      ! unit: mg Mic C/cm3
      use mic_constant
      use mic_variable
      implicit none
      real(r_2), parameter            :: sk =0.017
      real(r_2), parameter            :: skx=0.027
      real(r_2), parameter            :: ak = 10.0
      real(r_2), parameter            :: bk = 3.19
      real(r_2), parameter            :: xk1 =8.0
      real(r_2), parameter            :: xk2 =2.0
      real(r_2), parameter            :: xk3 =4.0
      real(r_2), parameter            :: xj1 =2.0
      real(r_2), parameter            :: xj2 =4.0
      real(r_2), parameter            :: xj3 =6.0
      real(r_2)                          xak
     TYPE(mic_parameter), INTENT(INOUT)   :: micparam
     TYPE(mic_input),     INTENT(IN)      :: micinput

  
      ! local variable
      real(r_2), dimension(mp,ms)     :: xkclay,km,kmx
      integer np,ns

      do np=1,mp
      do ns=1,ms
         xkclay(np,ns) = 1.0/ exp(-2.0*sqrt(micparam%clay(np,ns)))
         km(np,ns) =  xak * ak * exp(sk * micinput%tavg(np,ns) + bk)
         micparam%K1(np,ns) =  km(np,ns)/xk1
         micparam%K3(np,ns) =  km(np,ns) * xkclay(np,ns)/xk3
         micparam%J1(np,ns) =  km(np,ns)/xj1
         micparam%J3(np,ns) =  km(np,ns) * xkclay(np,ns)/xj3

         kmx(np,ns) =  xak * ak * exp(skx * micinput%tavg(np,ns) + bk)       
         micparam%K2(np,ns) =  kmx(np,ns)/xk2
         micparam%J2(np,ns) =  kmx(np,ns)/xj2
       enddo
       enddo

       if(diag==1) then   
          print *, 'Kmt',micparam%clay(outp,1),micinput%tavg(outp,1),km(outp,1),kmx(outp,1)
          print *, micparam%K1(outp,1)
          print *, micparam%K2(outp,1)
          print *, micparam%K3(outp,1)
          print *, micparam%J1(outp,1)
          print *, micparam%J2(outp,1)
          print *, micparam%J3(outp,1)   
       endif
   
    end subroutine Kmt

    subroutine Vmaxt(xav,micparam,micinput)
      ! mg Cs per mg mic C per hour
      use mic_constant
      use mic_variable
      implicit none
      real(r_2), parameter :: sv = 0.063
      real(r_2), parameter :: av = 8.0e-6
      real(r_2), parameter :: bv = 5.47
      real(r_2), parameter :: xv1= 10.0
      real(r_2), parameter :: xv2= 2.0
      real(r_2), parameter :: xv3= 10.0
      real(r_2), parameter :: xw1= 3.0
      real(r_2), parameter :: xw2= 3.0
      real(r_2), parameter :: xw3= 2.0
      real(r_2)  xav
      TYPE(mic_parameter), INTENT(INOUT)  :: micparam
      TYPE(mic_input),     INTENT(IN)     :: micinput
  

      ! local variables
      real(r_2),dimension(mp,ms) :: vmax
      integer np,ns

       do np=1,mp
          do ns=1,ms 
  
           vmax(np,ns) =  xav * av * exp(sv*micinput%tavg(np,ns) + bv) * delt
           micparam%V1(np,ns)   =  xv1 * vmax(np,ns) 
           micparam%V2(np,ns)   =  xv2 * vmax(np,ns) 
           micparam%V3(np,ns)   =  xv3 * vmax(np,ns) 
      
           micparam%W1(np,ns)   =  xw1 * vmax(np,ns) 
           micparam%W2(np,ns)   =  xw2 * vmax(np,ns)  
           micparam%W3(np,ns)   =  xw3 * vmax(np,ns) 
          enddo
        enddo
         
        if(diag==1) then 
           print *, 'Vmaxt',micinput%tavg(outp,1),vmax(outp,1)
           print *, micparam%V1(outp,1)
           print *, micparam%V2(outp,1)
           print *, micparam%V3(outp,1)
           print *, micparam%W1(outp,1)
           print *, micparam%W2(outp,1)
           print *, micparam%W3(outp,1)
        endif

    end subroutine Vmaxt

    subroutine Desorpt(xdesorp,micparam)
      use mic_constant
      use mic_variable
      implicit none
      real(r_2)              xdesorp
      TYPE(mic_parameter), INTENT(INOUT)    :: micparam 
      integer np,ns 

      do np=1,mp
      do ns=1,ms 
         micparam%desorp(np,ns) = xdesorp * (1.5e-5) * exp(-1.5*micparam%clay(np,ns)) 
      enddo
      enddo

      if(diag==1) then
         print *, 'Desorpt'
         print *, micparam%desorp(outp,:)
      endif
  
    end subroutine Desorpt

  subroutine mget(micparam,micinput,micnpool)
     use mic_constant
     use mic_variable
     implicit none
     real(r_2),    parameter  :: cuemax    = 0.80
     real(r_2),    parameter  :: cue_coef1 = 0.66
     real(r_2),    parameter  :: cue_coef2 = 1.23
 
     real(r_2),    parameter  :: epislon1 = 0.5
     real(r_2),    parameter  :: epislon2 = 0.25
     real(r_2),    parameter  :: epislon3 = 0.7
     real(r_2),    parameter  :: epislon4 = 0.35
     TYPE(mic_parameter), INTENT(INOUT)  :: micparam 
     TYPE(mic_input),     INTENT(IN)     :: micinput
     TYPE(mic_npool),     INTENT(IN)     :: micnpool 

     ! local variables
     integer np,ns

       do np=1,mp
       do ns=1,ms 
          ! variable mge 

          micparam%mgeR1(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,1)/micparam%cn_r(np,ns,3)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))

          micparam%mgeR2(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,2)/micparam%cn_r(np,ns,3)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))

          micparam%mgeR3(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,7)/micparam%cn_r(np,ns,3)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))

          micparam%mgeK1(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,1)/micparam%cn_r(np,ns,4)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))

          micparam%mgeK2(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,2)/micparam%cn_r(np,ns,4)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))

          micparam%mgeK3(np,ns) = cuemax*min(1.0,(micparam%cn_r(np,ns,7)/micparam%cn_r(np,ns,4)) &
                                  **(cue_coef1*(micnpool%mineralN(np,ns)-cue_coef2)))
         ! fixed mge
          micparam%mgeR1(np,ns) = epislon1
          micparam%mgeR2(np,ns) = epislon2
          micparam%mgeR3(np,ns) = epislon1
          micparam%mgeK1(np,ns) = epislon4
          micparam%mgeK2(np,ns) = epislon3
          micparam%mgeK3(np,ns) = epislon3
       enddo
       enddo

  end subroutine mget

  subroutine turnovert(xbeta,micparam,micinput)
      use mic_constant
      use mic_variable
      implicit none
      TYPE(mic_parameter), INTENT(INOUT)   :: micparam
      TYPE(mic_input),     INTENT(IN)      :: micinput  

      integer nx
      real(r_2)  xbeta
      real(r_2), parameter       :: xtv = 100.0
      real(r_2),  parameter      :: betamic  = 1.0   
  
      ! local variable
      integer np,ns
      real(r_2), dimension(mp)    :: tvref
 
        do np=1,mp
           tvref(np) = sqrt(micinput%fcnpp(np)/xtv)
           tvref(np) = max(0.6,min(1.3,tvref(np)))          ! 0.8-1.2 based on Wieder et al., 2015
           do ns=1,ms
              micparam%tvmicR(np,ns)   = 0.00052 * tvref(np) * exp(0.3 * micparam%fmetave(np,ns)) * delt
              micparam%tvmicK(np,ns)   = 0.00024 * tvref(np) * exp(0.1 * micparam%fmetave(np,ns)) * delt
              micparam%betamicR(np,ns) = betamic * xbeta
              micparam%betamicK(np,ns) = betamic * xbeta
           enddo
        enddo

        if(diag==1) then
          print *, 'turnovert'
          print *, micparam%tvmicR(outp,:) 
          print *, micparam%tvmicR(outp,:) 
        endif

  end subroutine turnovert


    subroutine bgc_fractions(micparam,micinput)
    use mic_constant
    use mic_variable
    implicit none
    real(r_2), parameter                :: fmicsom1=0.432
    real(r_2), parameter                :: fmicsom2=0.098
    real(r_2), parameter                :: fmicsom3=10.56
    real(r_2), parameter                :: fmicsom4=29.78
    real(r_2), parameter                :: fmicsom5=2.61
    TYPE(mic_parameter), INTENT(INOUT)  :: micparam  
    TYPE(mic_input),     INTENT(INOUT)  :: micinput
    !local variables
    integer np,ns
    real(r_2), dimension(mp)            :: fmetleaf,fmetroot,fmetwood
    real(r_2), dimension(mp,ms)         :: dleafx,drootx,dwoodx
    real(r_2), dimension(mp,ms)         :: cinputm,cinputs
    real(r_2), dimension(mp,ms,2)       :: cninp
    
    
      do np=1,mp
  
          fmetleaf(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligleaf(np) * micparam%xcnleaf(np)))
          fmetroot(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligroot(np) * micparam%xcnroot(np)))
          fmetwood(np) = max(0.0, 1.0 * (0.85 - 0.013 * micparam%fligwood(np) * micparam%xcnwood(np)))

          ! Initial C:N ratio of each C pool
          do ns=1,ms

             ! **this is a temporary solution, to be modified after N cycle is included
             micparam%cn_r(np,ns,1)=max( 5.0,0.5*(micparam%xcnleaf(np)+micparam%xcnroot(np)))
             micparam%cn_r(np,ns,2)=max(10.0,0.5*(micparam%xcnleaf(np)+micparam%xcnwood(np)))
             micparam%cn_r(np,ns,3)=7.4
             micparam%cn_r(np,ns,4)=13.4
             micparam%cn_r(np,ns,5)=12.0
             micparam%cn_r(np,ns,6)=16.0
             micparam%cn_r(np,ns,7)=10.0
 
 
       ! here zse in m, litter input in g/m2/delt, *0.001 to mgc/cm3/delt and "zse" in m.
             if(ns==1) then
                dleafx(np,ns) = 0.001* micinput%dleaf(np)/micparam%sdepth(np,1)                      ! mgc/cm3/delt
                drootx(np,ns) = 0.001* micparam%fracroot(np,1) * micinput%droot(np)/micparam%sdepth(np,1)     ! mgc/cm3/delt
                dwoodx(np,ns) = 0.001* micinput%dwood(np)/micparam%sdepth(np,1)                      ! mgc/cm3/delt
             else
                dleafx(np,ns) = 0.0
                drootx(np,ns) = 0.001 * micparam%fracroot(np,ns) * micinput%droot(np)/micparam%sdepth(np,ns)  ! mgc/cm3/delt
                dwoodx(np,ns) = 0.0
             endif

          !! calculate soil texture and litter quaility dependent parameter values
          ! C input to metabolic litter 
             micinput%cinputm(np,ns) = dleafx(np,ns)*fmetleaf(np)        &
                                     + drootx(np,ns)*fmetroot(np)        &
                                     + dwoodx(np,ns)*fmetwood(np)         
          ! C input to structural litter
             micinput%cinputs(np,ns) = dleafx(np,ns)*(1.0-fmetleaf(np))  &
                                     + drootx(np,ns)*(1.0-fmetroot(np))  &
                                     + dwoodx(np,ns)*(1.0-fmetwood(np)) 
             if((dleafx(np,ns)+drootx(np,ns)+dwoodx(np,ns))>0.0) then 
          ! C:N input of litter input to the metabolic pool 
                cninp(np,ns,1) = micinput%cinputm(np,ns)                          &
                               /(dleafx(np,ns)*fmetleaf(np)/micparam%xcnleaf(np)  &
                               +drootx(np,ns)*fmetroot(np)/micparam%xcnroot(np)   &
                               +dwoodx(np,ns)*fmetwood(np)/micparam%xcnwood(np))
          ! C:N input of litter input to the structural pool
                cninp(np,ns,2) = micinput%cinputs(np,ns)                               &
                               /(dleafx(np,ns)*(1.0-fmetleaf(np))/micparam%xcnleaf(np) &
                               +drootx(np,ns)*(1.0-fmetroot(np))/micparam%xcnroot(np)  &
                               +dwoodx(np,ns)*(1.0-fmetwood(np))/micparam%xcnwood(np))

                micparam%fmetave(np,ns) = (dleafx(np,ns)*fmetleaf(np) + drootx(np,ns)*fmetroot(np) + dwoodx(np,ns)*fmetwood(np))  &
                                        /(dleafx(np,ns) + drootx(np,ns) + dwoodx(np,ns)+ 1.0e-10)
             else
               if(ns==1) then
                  cninp(np,ns,1)          = micparam%xcnleaf(np)
                  cninp(np,ns,2)          = micparam%xcnleaf(np)
                  micparam%fmetave(np,ns) = fmetleaf(np)
               else
                  cninp(np,ns,1)          = micparam%xcnroot(np)
                  cninp(np,ns,2)          = micparam%xcnroot(np)
                  micparam%fmetave(np,ns) = fmetroot(np)
               endif
             endif

             micparam%cn_r(np,ns,1) = cninp(np,ns,1); micparam%cn_r(np,ns,2)=cninp(np,ns,2)

             micparam%fr2p(np,ns) = fmicsom1 * 0.30 * exp(1.3*micparam%clay(np,ns)) *1.0                   ! 3.0
             micparam%fk2p(np,ns) = fmicsom2 * 0.20 * exp(0.8*micparam%clay(np,ns)) *1.0                   ! 3.0
             micparam%fr2c(np,ns) = min(1.0-micparam%fr2p(np,ns), fmicsom3 * 0.10 * exp(-fmicsom5 * micparam%fmetave(np,ns))*1.0 )    ! 9.0   to invoid a negative value of fr2a  ZHC
             micparam%fk2c(np,ns) = min(1.0-micparam%fk2p(np,ns), fmicsom4 * 0.30 * exp(-fmicsom5 * micparam%fmetave(np,ns))*1.0)     ! 9.0   to invoid a negative value of fk2a ZHC
             micparam%fr2a(np,ns) = 1.00 - micparam%fr2p(np,ns) - micparam%fr2c(np,ns)
             micparam%fk2a(np,ns) = 1.00 - micparam%fk2p(np,ns) - micparam%fk2c(np,ns)
          enddo   !"ns"
      enddo       !"np"

      if(diag==1) then
         print *,'bgc_fraction parameters'
         print *, micinput%cinputm(outp,:)
         print *, micinput%cinputs(outp,:)
         print *, micparam%fmetave(outp,:)
         print *, micparam%cn_r(outp,:,1) 
         print *, micparam%cn_r(outp,:,2)
         print *, micparam%fr2p(outp,:) 
         print *, micparam%fk2p(outp,:) 
         print *, micparam%fr2c(outp,:)
         print *, micparam%fk2c(outp,:)
         print *, micparam%fr2a(outp,:) 
         print *, micparam%fk2a(outp,:)
      endif
   
   
   end subroutine bgc_fractions

   subroutine bioturb(ndelt,ms,zse,delt,diffsocxx,fluxsoc,xpooli,xpoole)
   ! multi-layered soil BGC including DOC and bioturbation using microbially-based BGC modeling
   ! step 1: litter-C and SOC bioturbation treated as a diffusion process
   ! step 2: advection of DOC along with water flux
   ! solve dc/dt=Dd2c/dx +F(z) where c is total SOC concentration in each soil layer
   ! bioturbation diffusion rate 
   ! boundary conditions: at the top     -Ddc/dx = F0+F(1)  at x=0
   !                      at the bottom: dC/dx=0            at x=h
   ! using the fully implicit method together with Thomas method
   ! unit for pool:                 mgc/cm3      (=kg C/m3)
   !      for flux:                 mgc/cm3/delt (=kg c/m3/delt): g/m2/delt = 0.1 mg/cm2/delt
   !      for length:               cm
   !      for diffsion coefficient: cm2/delt
   use mic_constant,  ONLY : r_2
   implicit none
   integer                        ndelt,ms
   real(r_2), dimension(ms)    :: zse
   real(r_2)                      delt,diffsocxx
   real(r_2), dimension(ms)    :: xpooli,xpoole,xpool,fluxsoc 
   ! local variables
   integer                        i,j
   real(r_2)                      deltD,tot0, tot1, totflux
   real(r_2), dimension(ms)    :: xzse
   real(r_2), dimension(ms+1)  :: sdepthx
   real(r_2)                      coeffA, coeffB
   real(r_2), dimension(ms)    :: at,bt,ct,rt
  
     ! calculate the mid-point of each layer
     sdepthx(1) = 0.0          ! depth of a layer from the top (x_0.5=0.0 eg soil surface)
     do j=2,ms+1
        sdepthx(j) = sdepthx(j-1) + zse(j-1)*100.0     ! depth of the bottom of each layer (eg x_j+0.5)
                                                       !*100 to convert from m to cm
     enddo

     do j=1,ms
        xzse(j) = 0.5 * (sdepthx(j) + sdepthx(j+1))    ! depth of midpoint of a layer j  (x_j)
     enddo

     deltD = diffsocxx * delt
     xpool = xpooli
     tot0 = 0.0
     do j=1,ms
        tot0 = tot0 + xpool(j) * zse(j)*100.0         !*100 convert m to cm
     enddo
  
     do i=1,ndelt
        do j=1,ms
           if(j==1) then
              coeffB = 1.0/(sdepthx(2)-sdepthx(1))
              coeffA = deltD*coeffB/(xzse(2)-xzse(1))
              ! Crank-Nicholson
              at(1) = 0.0
              bt(1) = 1.0 + 0.5 * coeffA
              ct(1) =     - 0.5 * coeffA
              rt(1) = (1.0-0.5*coeffA) * xpool(1) + 0.5 * coeffA * xpool(2) &
                    +  fluxsoc(1) * delt 
           endif
           if(j>1.and.j<ms) then
             coeffA = deltD/((xzse(j+1)-xzse(j))*(sdepthx(j+1)-sdepthx(j)))
             coeffB = (xzse(j+1)-xzse(j))/(xzse(j)-xzse(j-1))
             ! Crank-Nicholson
             at(j) =    -0.5 * coeffA * coeffB
             bt(j) = 1.0+0.5 * coeffA *(1.0+coeffB)
             ct(j) =    -0.5 * coeffA
             rt(j) = 0.5 * coeffA * coeffB * xpool(j-1)        &
                     +(1.0-0.5* coeffA*(1.0+coeffB))*xpool(j)  &
                     + 0.5* coeffA * xpool(j+1)                &
                     + fluxsoc(j) *delt
           endif
           if(j==ms) then
               coeffA = deltD/((xzse(ms)-xzse(ms-1))*(sdepthx(ms+1) - sdepthx(ms)))
             ! Crank-Nicholson
               at(ms) = -0.5 * coeffA
               bt(ms) = 1.0 + 0.5 * coeffA
               ct(ms) = 0.0
               rt(ms) = 0.5* coeffA  * xpool(ms-1) + (1.0-0.5 * coeffA) * xpool(ms) &
                      + fluxsoc(ms) * delt
           endif
        enddo
        call tridag(at,bt,ct,rt,xpool,ms)
     enddo
     xpoole = xpool
     
     tot1 = 0.0
     totflux=0.0
     do j=1,ms
        tot1 = tot1 + xpool(j) * zse(j) *100.0
        totflux = totflux + fluxsoc(j) * zse(j) *100.0
     enddo
  
end subroutine bioturb

   subroutine tridag(at,bt,ct,rt,u,ms)
   ! solving the triadigonal matrix (numerical recipes, p43)
   ! linear equation: A* u(i-1) + B *u(i) + C * u(i+1) = R, 
   ! where i is soil layer, u(i-1), u(i) and u(i+1) are at time step t
   ! NOTE: bt(1) should not be equal to 0.0, otherwise rewrite the equation
    use mic_constant,  ONLY : r_2
    implicit none
    integer, parameter    :: nmax=500
    integer ms
    real(r_2), dimension(ms)    :: at,bt,ct,rt,u
    integer j
    real(r_2) bet
    real(r_2), dimension(nmax) :: gam
     
      bet  = bt(1)
      u(1) = rt(1)/bet
      do j=2,ms
         gam(j) = ct(j-1)/bet
         bet = bt(j)-at(j)* gam(j)
         if(bet ==0) then
            print *, 'triag failed'
            stop
         endif
         u(j) = (rt(j) - at(j) * u(j-1))/bet
      enddo
      do j=ms-1,1,-1
         u(j) = u(j) -gam(j+1) * u(j+1)
      enddo
    end subroutine tridag


    subroutine advecdoc(deltx,zse,fluxsoilwx,fluxdocsx,vsoilwx,ypool)
    ! to be modified using an implicit solver to ensure mass conservation
    !
    use mic_constant
    implicit none
    real(r_2)                          deltx
    real(r_2), dimension(ms)        :: zse
    real(r_2), dimension(ms)        :: fluxsoilwx,vsoilwx,ypool
    real(r_2), dimension(ms)        :: dypool,ypool1
    real(r_2)                          fluxdocsx,totdoc0,totdoc1,fluxdocbot 
    integer ns,iter
     
     ypool1= ypool
     fluxdocbot = 0.0
     do iter=1,100
      do ns=1,ms
        if(ns==1) then
           dypool(1)  = (fluxdocsx - fluxsoilwx(1)*ypool1(1)/vsoilwx(1))*deltx*0.01/zse(1)
        else
           dypool(ns) = (fluxsoilwx(ns-1)*ypool1(ns-1)/vsoilwx(ns-1) &
                        -fluxsoilwx(ns)  *ypool1(ns)/vsoilwx(ns))       *deltx*0.01/zse(ns)
        endif
        if(ns==ms) then
           fluxdocbot = fluxdocbot + fluxsoilwx(ns)  *ypool1(ns)/vsoilwx(ns) *deltx* 0.01
        endif
      enddo
      ypool1 = max(0.0,ypool1+dypool)
     enddo
     ! check mass conservation
     totdoc0=0.0; totdoc1=0.0
     do ns=1,ms
        totdoc0 = totdoc0 + ypool(ns)  *zse(ns)
        totdoc1 = totdoc1 + ypool1(ns) *zse(ns)
     enddo
    ! print *, 'mass cons DOC', totdoc0,totdoc1,(totdoc1-totdoc0)-(fluxdocsx - fluxdocbot)*deltx
    
     ypool = ypool1
     
    end subroutine advecdoc


    subroutine vmic_c(np,ns,kinetics,micparam,micinput,xpool,y)
    ! MIMICS as modified by Zhang et al. (2019, GCB).
    ! Seven pools: metabolic litter (1), Structural litter, microbe-R (3), microbe-K(4),
    !              Physical protected (5), chemically-protected (6), active (7)
    ! for each layer
    ! input: aboveground litter and belowground in the surface layer
    !        belowground input in other layers
    ! kinetics: Michaelis-Mennten
    ! unit:
    ! all carbon pools : mg C/cm3
    ! time step :        one hour
    !
     use mic_constant
     use mic_variable
     implicit none
     real(r_2),  parameter                          ::  epislon1 = 0.5
     real(r_2),  parameter                          ::  epislon2 = 0.25
     real(r_2),  parameter                          ::  epislon3 = 0.7
     real(r_2),  parameter                          ::  epislon4 = 0.35
     real(r_2),  parameter                          ::  betamic  = 1.0   
     TYPE(mic_parameter), INTENT(IN)     :: micparam
     TYPE(mic_input),     INTENT(IN)     :: micinput

     real(r_2),  dimension(mcpool),  INTENT(IN)        :: xpool 
     real(r_2),  dimension(mcpool),  INTENT(INOUT)     :: y   !=dxpool/dt     ! local variables
     ! local variables
     integer     np,ns,kinetics  

     real(r_2)  betamicR,betamicK,                 &
                cinputmx,cinputsx,fmx,fsx,         &
                fr2px,fr2cx,fr2ax,                 &
                fk2px,fk2cx,fk2ax,                 &
                mgeRx1,mgeRx2,mgeRx3,              &
                mgeKx1,mgeKx2,mgeKx3,              &
                tvmicRx,tvmicKx,                   &
                tavgx,clayx,                       &
                desorpx,                           &
                V1x,V2x,V3x,W1x,W2x,W3x,           &
                J1x,J2x,J3x,K1x,K2x,K3x,           &
                Q1x,Q2x

     real(r_2) cfluxm2r, cfluxm2k, cfluxs2r, cfluxs2k, cfluxr,   cfluxk
     real(r_2) cfluxr2p, cfluxk2p, cfluxp2a, cfluxr2c, cfluxk2c
     real(r_2) cfluxc2a, cfluxr2a, cfluxk2a, cfluxa2r, cfluxa2k
 
      cinputmx = micinput%cinputm(np,ns);   cinputsx = micinput%cinputs(np,ns)
      tavgx    = micinput%tavg(np,ns);      clayx    = micparam%clay(np,ns)   

      fmx      = micparam%fm(np,ns);        fsx      = micparam%fs(np,ns)  
      fr2px    = micparam%fr2p(np,ns);      fr2cx    = micparam%fr2c(np,ns)
      fr2ax    = micparam%fr2a(np,ns);      fk2px    = micparam%fk2p(np,ns)
      fk2cx    = micparam%fk2c(np,ns);      fk2ax    = micparam%fk2a(np,ns)
      mgeRx1   = micparam%mgeR1(np,ns);     mgeRx2   = micparam%mgeR2(np,ns);   mgeRx3 = micparam%mgeR3(np,ns)
      mgeKx1   = micparam%mgeK1(np,ns);     mgeKx2   = micparam%mgeK2(np,ns);   mgeKx3 = micparam%mgeK3(np,ns)
      tvmicRx  = micparam%tvmicR(np,ns);    tvmicKx  = micparam%tvmicK(np,ns) 
      desorpx  = micparam%desorp(np,ns) 
      V1x      = micparam%V1(np,ns);        V2x      = micparam%V2(np,ns);      V3x    = micparam%V3(np,ns)
      W1x      = micparam%W1(np,ns);        W2x      = micparam%W2(np,ns);      W3x    = micparam%W3(np,ns)
      K1x      = micparam%K1(np,ns);        K2x      = micparam%K2(np,ns);      K3x    = micparam%K3(np,ns)
      J1x      = micparam%J1(np,ns);        J2x      = micparam%J2(np,ns);      J3x    = micparam%J3(np,ns)
      Q1x      = micparam%Q1(np,ns);        Q2x      = micparam%Q2(np,ns)            
      betamicR = micparam%betamicR(np,ns);  betamicK = micparam%betamicK(np,ns)

      ! carbon fluxes
      if(kinetics==1) then
        ! forward Michaelis-Menten
        cfluxm2r = xpool(3) * V1x * xpool(1)/(K1x + xpool(1))
        cfluxs2r = xpool(3) * V2x * xpool(2)/(K2x + xpool(2))
        cfluxa2r = xpool(3) * V3x * xpool(7)/(K3x + xpool(7))

        cfluxm2k = xpool(4) * W1x * xpool(1)/(J1x + xpool(1))
        cfluxs2k = xpool(4) * W2x * xpool(2)/(J2x + xpool(2))
        cfluxa2k = xpool(4) * W3x * xpool(7)/(J3x + xpool(7))

        cfluxr   = tvmicRx * xpool(3) ** betamicR
        cfluxk   = tvmicKx * xpool(4) ** betamicK

        cfluxr2p = fr2px * cfluxr
        cfluxk2p = fk2px * cfluxk

        cfluxr2c = fr2cx   * cfluxr 
        cfluxk2c = fk2cx   * cfluxk

        cfluxp2a = desorpx * xpool(5)
        cfluxr2a = fr2ax   * cfluxr
        cfluxk2a = fk2ax   * cfluxk
        cfluxc2a = xpool(3)* V2x * xpool(6)/(Q1x*K2x + xpool(6))   &
                 + xpool(4)* W2x * xpool(6)/(Q2x*J2x + xpool(6))
      endif
      if(kinetics ==2 )then 
        !=======================================================
        ! reverse Michaelis-Menten
        cfluxm2r = xpool(1) * V1x * xpool(3)/(K1x + xpool(3))
        cfluxs2r = xpool(2) * V2x * xpool(3)/(K2x + xpool(3))
        cfluxa2r = xpool(7) * V3x * xpool(3)/(K3x + xpool(3))

        cfluxm2k = xpool(1) * W1x * xpool(4)/(J1x + xpool(4))
        cfluxs2k = xpool(2) * W2x * xpool(4)/(J2x + xpool(4))
        cfluxa2k = xpool(7) * W3x * xpool(4)/(J3x + xpool(4))

        cfluxr   = tvmicRx * xpool(3) ** betamicR
        cfluxk   = tvmicKx * xpool(4) ** betamicK


        cfluxr2p = fr2px   * cfluxr
        cfluxk2p = fk2px   * cfluxk

        cfluxr2c = fr2cx * cfluxr 
        cfluxk2c = fk2cx * cfluxk

        cfluxp2a = desorpx * xpool(5)
        cfluxr2a = fr2ax * cfluxr
        cfluxk2a = fk2ax * cfluxk
        cfluxc2a = xpool(6) * V2x * xpool(3)/(Q1x*K2x + xpool(3))   &
                 + xpool(6) * W2x * xpool(4)/(Q2x*J2x + xpool(4))
      endif

      !===================================================
      ! 
      ! metabolic litter  [=Im*(1-fm)-A1-A5]
      y(1) = cinputmx * (1.0-fmx) - cfluxm2r - cfluxm2k

      ! structural litter [=Is*(1-fs)-A2-A6]
      y(2) = cinputsx * (1.0-fsx) - cfluxs2r - cfluxs2k

      !microbe R          [mge1*A1+mge2*A2+mge3*A3-A4]
      y(3) = mgeRx1 * cfluxm2r + mgeRx2 * cfluxs2r + mgeRx3 * cfluxa2r - cfluxr

      !microbe K          [mge3*A5+mge4*A6+mge3*A7-A8]
      y(4) = mgeKx1 * cfluxm2k + mgeKx2 * cfluxs2k + mgeKx2 * cfluxa2k - cfluxk


      !physically protected SOM: [Lm*fm+fpr*A4+fpk*A8-A9]
      y(5) = cinputmx * fmx + cfluxr2p + cfluxk2p - cfluxp2a 

      ! chemically protected SOM: [Is*fs+fcr*A4+fck*A8-A10]
      y(6) = cinputsx * fsx + cfluxr2c + cfluxk2c - cfluxc2a 

      !active SOM: [far*A4+fak*A8+A9+A10-A3-A7]
      y(7) = cfluxr2a + cfluxk2a + cfluxp2a + cfluxc2a - cfluxa2r - cfluxa2k
  
   end subroutine vmic_c 
   
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ for Newton_ralphson method@@@@@@@@@@@@@@@@@
!  these are not used 
    subroutine mnewt(ntrial,np,ns,kinetics,micparam,micinput,xpool0,tolx,tolf)
    use mic_constant,  ONLY : r_2,mp,ms,mcpool
    use mic_variable
    implicit none
    TYPE(mic_parameter), INTENT(IN)     :: micparam
    TYPE(mic_input),     INTENT(IN)     :: micinput
    integer      ntrial,np,ns,kinetics
    real(r_2)    tolf,tolx
!    real(r_2),   dimension(mp,ms,nxp),intent(inout)       :: xparam
    real(r_2),   dimension(mcpool),   intent(inout)       :: xpool0
    real(r_2),   dimension(mcpool)                        :: y,p
    real(r_2),   dimension(mcpool,mcpool)                 :: dydx
    integer i,k
    real(r_2)    d, errf,errx
    integer,    dimension(mcpool)                          :: index

       do k=1,ntrial
          call vmic_c(np,ns,kinetics,micparam,micinput,xpool0,y)
          call modeljacx(np,ns,kinetics,micparam,micinput,xpool0,dydx)
          errf = 0.0
101   format(i5,2x,10(f10.4,1x))
          do i=1,mcpool
             errf = errf + abs(y(i))
          enddo
          if(errf <=tolf) go to 999
          do i=1,mcpool
             p(i) = -y(i)
          enddo
          call ludcmp(dydx,mcpool,index,d)
          call ludksb(dydx,mcpool,index,p)
          errx = 0.0
          do i=1,mcpool
             errx = errx + abs(p(i))
             xpool0(i) = max(xpool0(i) + p(i), 1.0e-3) ! C pool cannot be negative values or 0 ZHC
          enddo
          if(errx <=tolx) go to 999
       enddo
999    continue
!       print *,'number of trials=', k
    end subroutine mnewt

    SUBROUTINE ludcmp(a,np,indx,d)
      use mic_constant,  ONLY : r_2
      implicit none 
      integer,      parameter     :: nmax=500
      real(r_2),    parameter     :: tiny =1.0e-20
      INTEGER np
      integer,      dimension(np)     :: indx
      real(r_2),    dimension(nmax)   :: vv
      real(r_2),    dimension(np,np)  :: a
      real(r_2)     d
      integer  i,imax,j,k
      real(r_2)     aamax,dum,sum

      d=1.
      do 12 i=1,np
        aamax=0.
        do 11 j=1,np
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
            print *, 'not paused: singular matrix in ludcmp!'
            aamax = 1.0e-4
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,np
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,np
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,np
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.np)then
          dum=1./a(j,j)
          do 18 i=j+1,np
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue

    END SUBROUTINE ludcmp

    SUBROUTINE ludksb(a,np,indx,b)
    use mic_constant,  ONLY : r_2
    IMPLICIT NONE
      INTEGER np,indx(np)
      REAL(r_2) a(np,np),b(np)
      INTEGER i,ii,j,ll
      REAL(r_2) sum

      ii=0
      do 12 i=1,np
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=np,1,-1
        sum=b(i)
        do 13 j=i+1,np
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue

   END SUBROUTINE ludksb     

    subroutine modeljacx(np,ns,kinetics,micparam,micinput,xpool0,dydx)
      use mic_constant,  ONLY : r_2,mp,ms,mcpool
      use mic_variable
      implicit none
      real, parameter          :: eps = 1.0e-4
      integer np,ns,kinetics
      TYPE(mic_parameter), INTENT(IN)     :: micparam
      TYPE(mic_input),     INTENT(IN)     :: micinput
      real(r_2), dimension(mcpool,mcpool),intent(out)   :: dydx
      real(r_2), dimension(mcpool),intent(inout)        :: xpool0
      real(r_2), dimension(mcpool)                      :: y0,y
      ! local variables
      real(r_2) h, temp
      integer i,j
     
       call vmic_c(np,ns,kinetics,micparam,micinput,xpool0,y0)
       do j=1,mcpool
          temp = xpool0(j)
          h    = eps * abs(temp)
          if(h<=0.0) h=eps
          xpool0(j) = temp + h
          h = xpool0(j) -temp
          call vmic_c(np,ns,kinetics,micparam,micinput,xpool0,y)  

          xpool0(j) = temp
          do i=1,mcpool
             dydx(i,j) = (y(i) - y0(i))/h ! the difference in changes of C pool caused by manupulated changes in C pools.
          enddo
       enddo

    end subroutine modeljacx

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
