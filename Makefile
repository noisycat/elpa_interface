##############################
################ cray - edison
CC=mpiicc -cc=icc
CXX=mpiicpc -cxx=icpc
FC=mpiifort -f90=ifort
#CXXFLAGS=-O0 -g -fopenmp -DDEBUG -DTEST4
CXXFLAGS=-O0 -g -fopenmp -DTEST4 -DLOCKING_TIMING  -DWITH_MPI
CFLAGS=$(CXXFLAGS) 
FCFLAGS=$(CXXFLAGS) 

ELPA_LIB_DIR=/opt/apps/intel/16.0.3/elpa/2016.05.003/lib
ELPA_INC_DIR=/opt/apps/intel/16.0.3/elpa/2016.05.003/include/elpa_openmp-2016.05.003

ELPA_LIB=-Wl,-rpath=$(ELPA_LIB_DIR) -L$(ELPA_LIB_DIR) -lelpa_openmp # $(IPM_GNU)
ELPA_INC=$(ELPA_INC_DIR)

##############################
##############################
################ workstations
#CXX=mpic++
#FC=mpif90
#CXXFLAGS=-O0 -g -openmp -DDEBUG -DTEST4
#FCFLAGS=$(CXXFLAGS) 
#ELPA_LIB=-Wl,-rpath=/workspace/elpa/lib -L/workspace/elpa/lib -lelpa  -mkl=cluster
#ELPA_INC=-I/workspace/elpa/include/elpa-2014.06.001 -I/workspace/elpa/include/elpa-2014.06.001/modules -mod /workspace/elpa/include/elpa-2014.06.001/modules
# main target
#all : main tests
all: elpa2_test_real_c_version

debug: 
	$(MAKE) all CXXFLAGS="$(CXXFLAGS) -DPRINT_GATHER_DEBUG -DPRINT_GATHER_DEBUG_N=34"

main : main.o  solve_provided.o test_input.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(ELPA_LIB)

cleanobj:
	rm -f *.o

cleanexe:
	rm -f main

cleantests:
	rm -f numroc_fortran_test test_real2 

cleanoutput:
	rm -f EVs_real2_out.txt EVs_c_out.txt output*.txt *.out *.txt


clean: cleanobj cleanexe cleantests cleanoutput

# implicit rules
%.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@ -I$(ELPA_INC) -I./ 

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(ELPA_INC)

%.o : %.F90
	$(FC) $(FCFLAGS) -I$(ELPA_INC) -module $(ELPA_INC)/modules -I$(ELPA_INC)/modules -c $< -o $@ 

# C++ has so many strange pitfalls. But we live and learn
main.o : main.cpp elpa_interface.hpp test_functions.hpp 

#linking their new one

utility.o: utility.F90
	$(FC) $(FCFLAGS) -c $< -o $@ 

elpa2_test_real_c_version: elpa2_test_real_c_version.c utility.o
	$(CC) $(CFLAGS) $^ -o $@ -I. -I$(ELPA_INC) $(ELPA_LIB) -Wl,-rpath=/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64 -L/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/ -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lifcoremt -lmkl_intel_lp64


numroc_fortran_test: numroc_fortran_test.o
	$(FC) $(FCFLAGS) $^ -o $@ $(ELPA_LIB)

test_real2 : test_real2.F90
	$(FC) $(FCFLAGS) $< -o $@ -I$(ELPA_INC) $(ELPA_LIB)

tests: numroc_fortran_test 

