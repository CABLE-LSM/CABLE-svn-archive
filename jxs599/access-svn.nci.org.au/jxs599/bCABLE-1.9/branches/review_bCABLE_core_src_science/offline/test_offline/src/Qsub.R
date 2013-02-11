
source( '../sites_main.nml' )

qs = c(2)
qj = c(length(sitenames))

for(j in 1:length(sitenames)) {

   qs[1] = paste('  sitenames = c( \'', sitenames[j], '\') ', sep='' )
   qs[2] = paste('  obsfiles = c( \'', obsfiles[j], '\') ', sep='' )

   write( qs, paste( 'qsmain_', j, '.nml', sep='' ) )
   
   qj[j] =j
   write( qj, 'qsj.j' )

}

#   system(paste('mv qsmain* ../')) 
#   system(paste('mv qsj.j ../')) 

