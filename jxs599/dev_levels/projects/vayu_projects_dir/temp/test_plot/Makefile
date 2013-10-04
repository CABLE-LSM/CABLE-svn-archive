#FC = gfortran 
FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
FLINK = -c 
FFLAGS = -O3 -fpp 
AUX= subr/
FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  
#NCI = -I /usr/local/include/  # -fbackslash
#NCMOD = /usr/local/include/NETCDF.mod # -fbackslash
#NCLIBS = -L/usr/local/lib/ -lnetcdf	## for mac

#NCLIBS = -L/usr/lib/ -lnetcdf
#FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng  
#PROGMAKE = transform 
#PROGEXE = ./transform 

test_driver: test_driver.f90
	$(FC) $(FFLAGS) -c module_parameters.f90	
	$(FC) $(FFLAGS) -c test_procedure.f90 $(FPGPLOTLIBS)	
	$(FC) $(FFLAGS) -o test_driver test_driver.f90	test_procedure.o	module_parameters.o $(FPGPLOTLIBS)






#make_cable: $(PROGMAKE)
#	$(PROGMAKE)


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
#   $(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_parameters.f90	
#test1.o : test1.f90 
#	$(FC) $(FFLAGS) -c test1.f90 
#
#test2.o : test1.f90 
#	$(FC) $(FFLAGS) -c test2.f90 
#
#
#exe_cable: $(PROGEXE)
#	$(PROGEXE)
#$(PROGEXE): test.o
#	$(FC) $(FFLAGS) -o $(PROG) test.o test1.o test2.o $(NCLIBS) 
#


##
## these are examples of working makefiles
##   
#$(PROGMAKE): transform.o
#	$(FC) $(FFLAGS) -o $(PROGMAKE) transform.o four1.o module_parameters.o 
#
#transform.o: transform.f90 four1.o module_parameters.o 
#	$(FC) $(FFLAGS) $(FLINK) -c transform.f90
#
#four1.o : four1.f90 module_parameters.o
#	$(FC) $(FFLAGS) $(FLINK) four1.f90 
#
#module_parameters.o : $(AUX)module_parameters.f90
#	$(FC) $(FFLAGS) $(FLINK) $(AUX)module_parameters.f90
# 
#

#new_fesc: new_fesc.f90
##	ln -s $(NCMOD)	netcdf.mod
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_parameters.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_vars.f90	
#		$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_txtread.f90	
#		$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_derivatives.f90	
#		$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_closest.f90	
#		$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_intercepts.f90	
#		$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_interp.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)plot_routines.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_data_setup.f90	
#	$(F) $(FFLAGS)	$(FLINK) 	$(AUX)module_spec_analysis.f90	
#	$(F) $(FFLAGS) -c new_fesc.f90 
##	$(F) $(FFLAGS) -o new_fesc new_fesc.o module_parameters.o module_vars.o module_derivatives.o module_closest.o module_intercepts.o module_interp.o module_data_setup.o module_spec_analysis.o  module_txtread.o plot_routines.o  $(NCLIBS) $(FPGPLOTLIBS) 
#	$(F) $(FFLAGS) -o new_fesc new_fesc.o module_parameters.o module_vars.o module_derivatives.o module_closest.o module_intercepts.o module_interp.o module_data_setup.o module_spec_analysis.o  module_txtread.o plot_routines.o  $(NCLIBS) $(FPGPLOTLIBS) 
#
#	
#	
#
#
#program_GammaBH: program_GammaBH.f90
#	ln -s $(NCMOD)	netcdf.mod
##	$(F)	$(FLINK) 	$(LocalSubroutines)module_txtread.f90	
##	$(F) $(FFLAGS)   program_GammaBH  program_GammaBH.f90 module_txtread.o $(FPGPLOTLIBS) 
#	$(F) $(FFLAGS) -c program_GammaBH.f90 
#	$(F) $(FFLAGS) -o program_GammaBH  program_GammaBH.o $(NCLIBS) 
#	
#	
#	
#	
#	 
#plot: plot.f90
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)genericread_module.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)plot_modules.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)plot_routines.f90	
##	$(F) $(FFLAGS)  plot plot.f90 plot_modules.o plot_routines.o $(FPGPLOTLIBS) 
#	$(F) $(FFLAGS)  plot plot.f90 plot_modules.o plot_routines.o genericread_module.o $(FPGPLOTLIBS) 
#
##plot: plot.f90
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)readcode_vars.f90	
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)plot_cases_module.f90	
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)plot_intmod.f90	
#  # $(F) $(FFLAGS)  plot plot.f90 readcode_vars.o  genericread_module.o plot_cases_module.o plot_modules.o plot_intmod.o  $(FPGPLOTLIBS) 
#
#match: match.f90
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)match_modules.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)genericread_modules.f90	
#	$(F) $(FFLAGS) match match.f90 match_modules.o genericread_modules.o $(FPGPLOTLIBS)
# 
##match: match.f90
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)genericread_modules.f90	
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)match_modules.f90	
##	$(F) $(FFLAGS) match match.f90 match_modules.o $(FPGPLOTLIBS)
# 
#test: test.f90
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)read_modules.f90	
#	$(F) $(FFLAGS) test test.f90 genericread_modules.o read_modules.o $(FPGPLOTLIBS) 
#
#
#readfile: readfile.f90
#	$(F) $(FFLAGS) readfile readfile.f90 read_modules.o $(FPGPLOTLIBS) 
#
# 
#program_filter:	program_filter.f90
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)filter_modules.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)write_filter_coords.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)write_filter.f90	
##	$(F)	$(FLINK) 	$(MYSUBROUTINES)read_FUNC.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)plot_filter.f90 
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)contour.f90	
#	$(F)	$(FFLAGS)	program_filter program_filter.f90 write_filter_coords.o write_filter.o filter_modules.o plot_filter.o contour.o $(FPGPLOTLIBS)  
#
#	
#write_filter : write_filter.f90
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)filter_modules.f90	
#	$(F)	$(FLINK) 	$(MYSUBROUTINES)write_filter_coords.f90	
##	$(F) $(FFLAGS)  write_filter  write_filter.f90 filter_modules.o write_filter_coords.o $(FPGPLOTLIBS) 
#
clean:
	rm -f *.o *~ *.mod






