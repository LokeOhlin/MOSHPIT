##-------------------------------------------------------------------##
##     Makefile for compiling the dustywave example program          ##
##     Daniel Price and Guillaume Laibe                              ##
##     daniel.price@monash.edu, guillaume.laibe@monash.edu           ##
##-------------------------------------------------------------------##

F90  = gfortran
F2PY = f2py
F90FLAGS = "-O3 -Wall"

# define the implicit rule to make a .o file from a .f90 file
%.o : %.f90
	$(F90C) -c $(F90FLAGS) $< -o $@

SRC= cubicsolve.f90 dustywaves.f90
OBJ=$(SRC:.f90=.o)

dustywave: $(SRC)
	$(F2PY) --f90exec=$(F90) --f90flags=$(F90FLAGS) -c -m dustywaves $(SRC)
clean:
	rm *.o *.mod dustywave

cleanall:
	rm *.o *.mod dustywave dustywave*.out
