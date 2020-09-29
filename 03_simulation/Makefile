
f90 = gfortran

flags=-ff2c -ftree-vectorize -fimplicit-none -fbackslash -O3 -Wall

ompflags=-fopenmp -ftree-vectorize -fimplicit-none -fbackslash -O3 -Wall

all: MD.IS-SPA.multi.sym.f
	$(f90) $(ompflags) -o MD.IS-SPA.multi.sym.e MD.IS-SPA.multi.sym.f


