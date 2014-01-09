#FFLAGS = -O
#FC = gfortran 
FC = ifort 
FLINK = -c 
#NCI = -I /usr/local/include/  # -fbackslash
#NCMOD = /usr/local/include/NETCDF.mod # -fbackslash
FFLAGS = -O0 # -Wall 
#AUX= FLOCAL/
#NCLIBS = -L/usr/local/lib/ -lnetcdf	## for mac
#FPGPLOTLIBS = -L/sw/lib/pgplot/ -lpgplot -L/sw/lib/ -lpng  ## for mac
#NCLIBS = -L/usr/lib/ -lnetcdf
#FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng  
#PROG = ./runtest

#   rm -f netcdf.mod

## this is how i have always written makw=efiles
test : driver.f90
	$(FC) $(FFLAGS)	$(FLINK) 	driver.f90	
	$(FC) $(FFLAGS)	$(FLINK) 	subr.f90	
	$(FC) $(FFLAGS) -o test driver.o subr.o $(NCLIBS) 
### or just
#test : test.f90
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test1.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test2.f90	
#	$(F) $(FFLAGS) -o test test.f90 test1.o test2.o $(NCLIBS) 

## alternatively makefiles can be written in this format as cable is done
#cable_make: $(PROG)
#	$(PROG)
#
#$(PROG): test.o
#	$(FC) $(FFLAGS) -o $(PROG) test.o test1.o test2.o $(NCLIBS) 
#
#test.o : test.f90 test1.o test2.o
#	$(FC) $(FFLAGS) -c test.f90 
#
#test1.o : test1.f90 
#	$(FC) $(FFLAGS) -c test1.f90 
#
#test2.o : test2.f90 
#	$(FC) $(FFLAGS) -c test2.f90 
#
#
#









clean:
	rm -f *.o *~ *.mod






