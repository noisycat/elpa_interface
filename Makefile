CXX=mpic++
FC=mpif90
CXXFLAGS=-O0 -g -openmp -DDEBUG
FCFLAGS=$(CXXFLAGS) 
ELPA_LIB=-Wl,-rpath=/workspace/elpa/lib -L/workspace/elpa/lib -lelpa  -mkl=cluster
ELPA_INC=-I/workspace/elpa/include/elpa-2014.06.001 -I/workspace/elpa/include/elpa-2014.06.001/modules -mod /workspace/elpa/include/elpa-2014.06.001/modules
# main target
all : main tests

main : main.o  solve_provided.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(ELPA_LIB)

cleanobj:
	rm -f *.o

cleanexe:
	rm -f main

cleantests:
	rm -f numroc_fortran_test test_real2

cleanoutput:
	rm -f EVs_real2_out.txt EVs_c_out.txt


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

test_real2 : test_real2.F90
	$(FC) $(FCFLAGS) $< -o $@ $(ELPA_INC) $(ELPA_LIB)

tests: numroc_fortran_test test_real2

