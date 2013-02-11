
outfile_func <- function( outfilename, outtype, sitename ) {        
   if(outtype=='pdf')
	   pdf(file=outfilename,paper='a4r',width=11,height=8) #filetype 
   if(outtype=='ps')
		postscript(file=outfilename,paper='special',width=11,height=8)
   if(outtype=='png')
		png(file=outfilename,width=2000,height=1414,pointsize=24)
   if(outtype=='jpg')
		jpeg(file=outfilename,width=4000,height=2828,pointsize=24)
}

# summing functions
fsumTWOd <- function( avcol, data, a, b, i, j ) {
   sub =  sum(data[(a+(j-1)*365):(b+(j-1)*365),i]) 
   val =  avcol + sub
   list( cummx=val, ittx=sub )
}

fsumONEd <- function( avcol, data, a, b, j ) {
   sub =  sum(data[(a+(j-1)*365):(b+(j-1)*365)]) 
   val =  avcol + sub
   list( cummx=val, ittx=sub )
}

func_ochk <- function( onc, flux ) {
   if (flux == 'Rnet') {
      if (is.character(onc$var$Rnet$name)) check = TRUE
      else check = FALSE
   }
   else  check = TRUE
}



#quantify goodness of fit 
compute_chisq <- function( n, check, O, E, scatter ) {
   dof = n -1 #degrees of freedom
   if( check ) { 
      sub=0.0
      for(i in 1:n) {
         DeltaO = ( O[i]-E[i] ) 
         itt = ( DeltaO * DeltaO ) / ( scatter[i] * scatter[i] ) # check this against qc stuff 
         if ( scatter[i] == 0) sub = sub
         else sub = sub + itt
#         print( paste( 'chi chi_i, obs, exp, sig ',sub, itt, O[i],E[i],scatter[i]   ) )
      } 
      val = (sub) / dof   #in principal this should be dof but as sample size is only small this 
                        #warps result from 1 (for n) for good fit
   }
   else val <- 100000.0
   chisq <-as.numeric( val )
   return(chisq)
}

##############################################################################################
### generic function to compute statistical scatter in data arrays passed, and upper/lower ###
### bounds on the data which can be plotted if so desired.                                 ###  
##############################################################################################
compute_scatter <-function( n, nyears, obs, obs_i ) { 
	sigma=c() 
	plus_sigma=c() 
	minus_sigma=c() 
   dof = nyears -1 
  	for(l in 1:n){ 
      scatter_sub=0.0
	   for(k in 1:nyears){ 
         scatter_t1 =  obs_i[k,l] -  obs[l]
         scatter_t2 = scatter_t1 * scatter_t1 
         scatter_sub = scatter_sub + scatter_t2 
      } 
      sigma[l] = sqrt( scatter_sub / nyears )
     # print( paste('sigma[l]  ',sigma[l], obs_i[k,l], obs[l] ) )
   }
   minus_sigma <-  obs - sigma
   plus_sigma<-  obs + sigma
   list( sigma=sigma, plus_sigma=plus_sigma, minus_sigma=minus_sigma )
}











#df=11 here # iwasn't able to make this dynamic
chisq_pdf_func<- function ( chi_x, nsteps, nvals ) {
   df =nvals-1 #=11 for seasonal #i.e. 12 - 1 
   chisq = array( 0,c(1,nsteps) ) #chi sq. pdf - dec.array (1*100) & init. to 0
   fGamma <- function(t) { z=df/2 ;t^(z-1) * exp(-t) } #t is the dummy argument of integration
   Gamma <- integrate( fGamma, lower = 0, upper=Inf )
   termG <- as.numeric(Gamma[1])
   term1 = 2^(df/2)
   for( i in 1:nsteps) {
      term2 = chi_x[i]^( (df/2) -1 )
      term3 = exp( -chi_x[i] /2 )
      chisq[i] = ( 1/ (term1*termG) ) * term2 * term3
   #   print( paste('integral ', chi_x[i], chisq[i], term1, termG, term2, term3 ) )
   }
   return(chisq)
}
























