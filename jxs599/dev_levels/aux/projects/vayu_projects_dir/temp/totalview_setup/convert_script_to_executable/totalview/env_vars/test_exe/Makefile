FC = ifort
FFLAGS = -O0 # -Wall 
AUX= FLOCAL/
PROG = um_run_tv


$(PROG): driver.f90
	$(FC) $(FFLAGS) -o $(PROG) driver.f90 

clean:
	rm -f *.o *~ *.mod






