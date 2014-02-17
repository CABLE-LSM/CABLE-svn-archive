
#FC = gfortran 
# xe
#FC = /apps/intel-fc/11.0.074/bin/intel64/ifort
#cherax
FC = ifort
FLINK = -c 
CFLAGS = -O0 -fpp
FFLAGS = -O0 -fpp -r8 -i8 
AUX= subr/
#cherax
#FPGPLOTLIBS = -L/usr/lib/ -lpgplot -L/usr/lib/ -lpng
#FPGPLOTLIBS = -L/cs/home/csdar/srb001/Public/pgplot -lpgplot -L/usr/lib -lpng -L/usr/lib -lz -L/usr/X11R6/lib/ -lX11  
# xe
#FPGPLOTLIBS = -L/home/599/jxs599/Public/pgplot/ -lpgplot -L/usr/lib64 -lpng -L/usr/lib64 -lz -L/usr/X11R6/lib64/ -lX11  

OBJ= debug_common.o debug_read.o debug_write.o cat_nodes.o

diag_cat_nodes_r8: diag_cat_nodes.f90
	$(FC) $(FFLAGS) -c $(AUX)debug_common.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_read.f90	
	$(FC) $(FFLAGS) -c $(AUX)debug_write.f90	
	$(FC) $(FFLAGS) -c $(AUX)cat_nodes.f90	
	$(FC) $(FFLAGS) -o diag_cat_nodes_r8 diag_cat_nodes.f90 $(OBJ) $(FPGPLOTLIBS)

clean:
	rm -f *.o *~ *.mod






