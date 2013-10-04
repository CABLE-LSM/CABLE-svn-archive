FFLAGS = -O0 
FC = ifort
FLINK = -c 

FFLAGS = -O0 # -Wall 
AUX= FLOCAL/



#   rm -f netcdf.mod

## this is how i have always written makw=efiles
#test : test.f90
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test1.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test2.f90	
#	$(F) $(FFLAGS) -c test.f90 
#	$(F) $(FFLAGS) -o test test.o test1.o test2.o $(NCLIBS) 
### or just
#test : test.f90
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test1.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)test2.f90	
#	$(F) $(FFLAGS) -o test test.f90 test1.o test2.o $(NCLIBS) 

## alternatively makefiles can be written in this format as cable is done


program: program.F90
#	$(FC) $(FFLAGS) -c vars.f90 
	$(FC) $(FFLAGS) -o program program.F90 # vars.o 


test: test.f90
	$(FC) $(FFLAGS) -o test test.f90 

	
write_filter : write_filter.f90
	$(F)	$(FLINK) 	$(MYSUBROUTINES)filter_modules.f90	
	$(F)	$(FLINK) 	$(MYSUBROUTINES)write_filter_coords.f90	
#	$(F) $(FFLAGS)  write_filter  write_filter.f90 filter_modules.o write_filter_coords.o $(FPGPLOTLIBS) 

clean:
	rm -f *.o *~ *.mod






