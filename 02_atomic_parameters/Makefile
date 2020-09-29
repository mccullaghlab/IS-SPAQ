
f90 = gfortran

flags=-ff2c -ftree-vectorize -fimplicit-none -fbackslash -O3 -Wall

ompflags=-fopenmp -ftree-vectorize -fimplicit-none -fbackslash -O3 -Wall

all: extract_LJparm.f SPA.poisson.omp.rc12.sym.f SPA.poisson2.omp.rc12.sym.f SPA.pol2.omp.rc12.sym.f
	$(f90) $(flags) -o extract_LJparm.e extract_LJparm.f
	$(f90) $(ompflags) -o SPA.poisson.omp.rc12.sym.e SPA.poisson.omp.rc12.sym.f
	$(f90) $(ompflags) -o SPA.poisson2.omp.rc12.sym.e SPA.poisson2.omp.rc12.sym.f
	$(f90) $(ompflags) -o SPA.pol2.omp.rc12.sym.e SPA.pol2.omp.rc12.sym.f


