
comp_data <- function( rdata ) {
   n_pfts_orig <- 5
   pseudo_len <- length(rdata$pseudo ) 

   #for test purposes
   Lnewpseudo_len <- pseudo_len 
   #Lnewpseudo_len <- 17

   #---dec. array
   Lnewpseudo <- array( 0, c(Lnewpseudo_len) )

   #because 4th dim (t) = 1 therefore only 3 dimensions here
   Lnewfield <- array( 0, c( length(rdata$lon), length(rdata$lat), Lnewpseudo_len ) )

   if(Lnewpseudo_len == pseudo_len) {
       for(i in 1:Lnewpseudo_len) {
         Lnewpseudo[i] = i 
         Lnewfield[,,i] = rdata$field[,,i] 
      }
   }

   if(Lnewpseudo_len > pseudo_len) {
      #print( paste('dim field', dim(rdata$field) ) )
      for(i in 1:Lnewpseudo_len) {
         Lnewpseudo[i] = i 
         #in all cases below alloc. of last tiles is common
         
         #---rewrite ancil. w 1st and last matching orig. and all 
         #---tiles added in between w no effect
         #-----------------------------------------
         #---1st 5 tiles repeated in postions 1:5 
         if(i <= n_pfts_orig) Lnewfield[,,i] = rdata$field[,,i] 
         #---tiles6:13 set to 0 (i.e. they have no effect)
         if(i > n_pfts_orig & i <= (n_pfts_orig+Lnewpseudo_len - pseudo_len) )  Lnewfield[,,i] = 0.0 

         #---rewrite ancil. w last matching orig. 1st 5tiles repeated twice - BUT
         #---half the frac. (net result should be same) 
         #-----------------------------------------
         #---1st 5 tiles repeated in postions 1:5 
         #if(i <= 5) Lnewfield[,,i] = 0.5 * rdata$field[,,i] 
         #---1st 5 tiles repeated in postions 6:10
         #if(i > 5 & i <= 10 ) Lnewfield[,,i] = 0.5 * rdata$field[,,i-5] 
         #---tiles in postions 11:13 set to 0 (i.e. they have no effect)
         #if(i > 10 & i <= 13 )  Lnewfield[,,i] = 0.0 
    
         #---tiles13:17 set to original 5:9(non-veg tiles)
         if(i > (n_pfts_orig+Lnewpseudo_len - pseudo_len) ) Lnewfield[,,i] = rdata$field[,,i-(Lnewpseudo_len - pseudo_len)] 
      }
   }

#   print( paste('zero base   ', rdata$field[0] ) )
#   print( paste('one base   ', rdata$field[1] ) )
   list(  newpseudo=Lnewpseudo, newfield=Lnewfield ) 
}


