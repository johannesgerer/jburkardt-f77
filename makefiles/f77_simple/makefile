#
#  On Chili and Kachina, this is where the FORTRAN77 compiler is:
#
F77 = /usr/bin/g77

f77_simple.out: f77_simple
	./f77_simple > f77_simple.out
f77_simple: f77_simple.o midpoint.o f.o
	$(F77) -o f77_simple f77_simple.o midpoint.o f.o
f77_simple.o: f77_simple.f
	$(F77) -c f77_simple.f
midpoint.o: midpoint.f
	$(F77) -c midpoint.f
f.o: f.f
	$(F77) -c f.f
clean:
	/bin/rm f77_simple.out f77_simple f77_simple.o midpoint.o f.o



