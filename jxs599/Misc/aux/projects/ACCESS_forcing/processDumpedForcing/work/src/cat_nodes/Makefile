#program specific structure
AUX= subr/

#compiler independent
FLINK = -c 

#compiler specific
FC = gfortran 
FFLAGS = -O0 -x f95-cpp-input  -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8     

#FC = ifort
#FFLAGS = -O0 -fpp -r8 -i8 

##FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng
#cherax?
#FPGPLOTLIBS = -L/cs/home/csdar/srb001/Public/pgplot -lpgplot -L/usr/lib -lpng -L/usr/lib -lz -L/usr/X11R6/lib/ -lX11  
# xe
#FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  

OBJ= debug_common.o debug_read.o debug_write.o cat_nodes.o

diag_cat_nodes_r8: diag_cat_nodes.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_write.f90	
	$(FC) $(FFLAGS) -c $(AUX)cat_nodes.f90	
	$(FC) $(FFLAGS) -c diag_cat_nodes.f90 
	$(FC) $(FFLAGS) -o diag_cat_nodes_r8 diag_cat_nodes.o $(OBJ) 
clean:
	rm -f *.o *~ *.mod






