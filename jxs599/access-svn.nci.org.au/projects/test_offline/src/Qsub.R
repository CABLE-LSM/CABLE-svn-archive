
source( '../sites_main.txt' )

qs = c(2)
qj = c(length(sitenames))

for(j in 1:length(sitenames)) {

   qs[1] = paste('  sitenames = c( \'', sitenames[j], '\') ', sep='' )
   qs[2] = paste('  obsfiles = c( \'', obsfiles[j], '\') ', sep='' )

   write( qs, paste( 'qsmain_', j, '.txt', sep='' ) )
   
   qj[j] =j
   write( qj, 'qsj.j' )

}


