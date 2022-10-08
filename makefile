# This is a commentary line in a makefile
# Start of the makefile

exe1: main.o user.o grids.o cfdlib.o visual.o rans_models.o
	gfortran -o exe1 main.o user.o grids.o visual.o cfdlib.o rans_models.o

user.mod: user.o user.f90
	gfortran -c user.f90

user.o: user.f90
	gfortran -c user.f90

rans_models.mod: rans_models.o rans_models.f90 user.mod
	gfortran -c rans_models.f90 

rans_models.o: rans_models.f90 user.mod
	gfortran -c rans_models.f90 

cfdlib.mod: cfdlib.o cfdlib.f90 user.mod rans_models.mod
	gfortran -c cfdlib.f90

cfdlib.o: cfdlib.f90 user.mod rans_models.mod
	gfortran -c cfdlib.f90

grids.mod: grids.o grids.f90 user.mod 
	gfortran -c grids.f90

grids.o: grids.f90 user.mod
	gfortran -c grids.f90

visual.mod: visual.o visual.f90 user.mod 
	gfortran -c visual.f90

visual.o: visual.f90 user.mod
	gfortran -c visual.f90

main.o: user.mod main.f90 grids.mod cfdlib.mod visual.mod rans_models.mod
	gfortran -c main.f90
clean:
	rm user.mod user.o main.o grids.o exe1 ex1_.dat grids.mod visual.o visual.mod cfdlib.o cfdlib.mod
# End of the makefile
