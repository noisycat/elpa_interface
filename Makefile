##############################
################ cray - edison
CXX=CC
FC=ftn
#CXXFLAGS=-O0 -g -fopenmp -DDEBUG -DTEST4
CXXFLAGS=-O3 -g -fopenmp -DTEST4 -DLOCKING_TIMING -DPRINT_GATHER_DEBUG
FCFLAGS=$(CXXFLAGS) 
ELPA_LIB=-Wl,-rpath=$(HOME)/gcc/elpa/lib -L$(HOME)/gcc/elpa/lib -lelpa_openmp $(IPM_GNU)
ELPA_INC=-I$(HOME)/elpa-2014.06.001 -I$(HOME)/elpa-2014.06.001/modules -J$(HOME)/elpa-2014.06.001/modules
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
all : main tests

main : main.o  solve_provided.o test_input.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(ELPA_LIB)

cleanobj:
	rm -f *.o

cleanexe:
	rm -f main

cleantests:
	rm -f numroc_fortran_test test_real2 test_hybrid

cleanoutput:
	rm -f EVs_real2_out.txt EVs_c_out.txt output*.txt


clean: cleanobj cleanexe cleantests cleanoutput

# implicit rules
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%.o : %.F90
	$(FC) $(FCFLAGS) -c $< -o $@ $(ELPA_INC)

# C++ has so many strange pitfalls. But we live and learn
main.o : main.cpp elpa_interface.hpp test_functions.hpp 

numroc_fortran_test: numroc_fortran_test.o
	$(FC) $(FCFLAGS) $^ -o $@ $(ELPA_LIB)

test_hybrid : test.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(ELPA_INC) $(ELPA_LIB)

test_real2 : test_real2.F90
	$(FC) $(FCFLAGS) $< -o $@ $(ELPA_INC) $(ELPA_LIB)

tests: numroc_fortran_test test_hybrid 

