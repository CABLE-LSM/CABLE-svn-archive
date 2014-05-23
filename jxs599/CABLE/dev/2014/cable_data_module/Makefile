#program specific structure
AUX= subr/

#compiler independent
FLINK = -c 

#compiler specific
FC = /usr/local/bin/gfortran 
FFLAGS = -O0 -x f95-cpp-input  #-fdefault-real-8 -fdefault-double-8 -fdefault-integer-8     

#FC = ifort
#FFLAGS = -O0 -fpp -r8 -i8 

##FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng
#cherax?
#FPGPLOTLIBS = -L/cs/home/csdar/srb001/Public/pgplot -lpgplot -L/usr/lib -lpng -L/usr/lib -lz -L/usr/X11R6/lib/ -lX11  
# xe
#FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  

OBJ= cable_data_module.o 
cable_data_module.o: cable_data_module.F90
	$(FC) $(FFLAGS) -c cable_data_module.F90 
clean:
	rm -f *.o *~ *.mod






