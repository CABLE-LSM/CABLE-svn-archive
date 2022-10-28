! this file includes codes for initilization, input and output related to mic_cnp.F90
! Strategy for implementing MIMICS into CABLE
!
!(1) allocate all the variables at the start
!(2) read the parameter lookup table
!    make the following parameters with (mp) dimensions and add to the "allocate"
!    xav,xak,xdesorp,xbeta,xdiffsoc and rootbetax
!(3) assign all the parameter values to 1:mp
!(4) each timestep: get the litter input from the previous day, and current time 
!    of soil temperature, moisture
!(5) output soil respiration and update pool sizes
!(6) add "mic" as a switch (=1,2 and 3 for C, CN and CNP cycles)
!(7) call "vmic" from "biogeochem" in "casa_inout.F90"
! strategy based on casa-cnp in CABLE
! (8) check the unit of "deltvmic" for consistency
!
!===========================================================================================
MODULE vmic_inout_mod
  USE cable_def_types_mod, ONLY : mp,ms,r_2,mvtype,mstype
  IMPLICIT NONE

subroutine mic_init(micparam,micinput,micoutput,miccpool,micnpool)
   use mic_constant
   use mic_variable
   implicit none
   TYPE(mic_parameter)  :: micparam
   TYPE(mic_input)      :: micinput
   TYPE(mic_cpool)      :: miccpool
   TYPE(mic_npool)      :: micnpool
   TYPE(mic_output)     :: micoutput
   
      call mic_allocate_parameter(micparam)
      call mic_allocate_input(micinput)
      call mic_allocate_output(micoutput)
      call mic_allocate_cpool(miccpool)
      call mic_allocate_npool(micnpool)
	  
end subroutine mic_init	 


subroutine mic_parameter(veg,soil,casabiome,zse)
! read the parameter lookup table and assigm them to (1:mp)
   use mic_constant
   use mic_variable
   implicit none
   TYPE(mic_parameter)                    :: micparam
   TYPE(mic_input)                        :: micinput
   TYPE (veg_parameter_type),  INTENT(IN) :: veg
   TYPE (soil_parameter_type), INTENT(IN) :: soil
   TYPE (casa_biome),          INTENT(IN) :: casabiome
   real(r_2), dimension(ms)               :: zse
   real(r_2), dimension(mpft)             :: xav_pft,      &
                                             xak_pft,      &
                                             xavexp_pft,   &
                                             xbeta_pft,    &
                                             xdiffsoc_pft, &
                                             xdesorpt_pft, &
                                             rootbetax_pft
								
! parameter estimates based on ICP-China sites for 4 PFTs based optimizinfg 5 parameters
! xav      xak       xavexp    xbeta    xdiffsoc   
! 14.524    15.091     1.008     1.255    18.329
!  5.475    13.271     1.011     1.446     4.636
!  6.125     9.869     1.324     1.582     1.051
!  4.757    19.013     1.271     1.955     0.263
!
!  7.72     14.31      1.15      1.56      6.07

   data      xav_pft/ 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72, 7.72/
   data      xak_pft/14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31,14.31/
   data   xavexp_pft/ 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15, 1.15/
   data    xbeta_pft/ 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56, 1.56/
   data xdiffsoc_pft/ 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07, 6.07/
   data rootbeta_pft/ 0.96, 0.96, 0.96, 0.98, 0.96, 0.97, 0.94, 0.96, 0.97, 0.94, 0.97, 0.97, 0.97, 0.97, 0.96, 0.96, 0.96/   

   ! local variables
   integer ivt,np,ns
   
   xdesorpt_pft(:) = 1.0
   do np=1,mp
      ivt   = veg%iveg(np))
	  isoil = soil%isoilm(np)
      xav(np)      = xav_pft(ivt)
      xak(np)      = xak_pft(ivt)
      xavexp(np)   = xavexp_pft(ivt)	  
      xbeta(np)    = xbeta_pft(ivt)
      xdiffsoc(np) = xdiffsoc_pft(ivt)
      rootbeta(np) = rootbeta_pft(ivt)
      xdesorpt(np) = xdesorpt_pft(ivt)
      ! soil parameters	  
      micparam%bulkd(np) = soil%rhosoil(np)
      micparam%clay(np)  = soil%clay(np)	
      ! vegetation parameters
       
      micparam%pft(np)     =  veg%iveg(np))	   
      micparam%xcnleaf(np) = 1.0/casabiome%ratioNCplantmax(ivt,1)
      micparam%xcnroot(np) = 1.0/casabiome%ratioNCplantmax(ivt,3)	  
      micparam%xcnwood(np) = 1.0/casabiome%ratioNCplantmax(ivt,2)
      micparam%fligleaf(np)=  casabiome%fracligninplant(ivt,1)
      micparam%fligroot(np)=  casabiome%fracligninplant(ivt,3)
      micparam%fligwood(np)=  casabiome%fracligninplant(ivt,2)
    enddo  
    do ns=1,ms
       zse(ns) = soil%zse(ns)
    enddo
	
end subroutine mic_parameter

subroutine mic_input(dleaf,dwood,droot,tsoil,moist,nsoilmin)
   use mic_constant
   use mic_variable
   implicit none
   TYPE(mic_parameter), INTENT(IN)          :: micparam
   TYPE(mic_input),     INTENT(OUT)         :: micinput
   TYPE(mic_npool),     INTENT(IN)          :: micnpool   
   real(r_2),  dimension(mp)                :: dleaf,dwood,droot,fcnpp,nsoilmin
   real(r_2),  dimension(mp,ms)             :: tsoil,moist
   real(r_2),  dimension(mpft)              :: fcnpp_pft
   data fcnpp_pft/500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0/
 
   ! local variables
   integer np,ivt
  
     do np=1,mp
	    !
        ivt =micparam%pft(np)
	    micinput%fcnpp(np) = fcnpp_pft(ivt)  ! gc/m2/year for calculating microbial turnover rate (constant!!)
		micinput%dleaf(np) = dleaf(np)       ! gc/m2/delt
		micinput%dwood(np) = dwood(np)       ! gc/m2/delt
        micinput%droot(np) = droot(np)       ! gc/m2/delt
		micinput%tavg(np,:)= tsoil(np,:)
		micinput%wavg(np,:)= moist(np,:)
		micnpool%mineralN(np,1:ms) = nsoilmin(np)   ! temporary solution
	 enddo	
end subroutine mic_input

END MODULE vmic_inout_mod