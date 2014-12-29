CXX=mpic++
CXXFLAGS=-O2 -openmp
ELPA_LIB=-Wl,-rpath=/workspace/elpa/lib -L/workspace/elpa/lib -lelpa  -mkl=cluster
# main target
main : main.o elpa_interface.o solve_provided.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(ELPA_LIB)

clean:
	rm -f *.o
	rm -f main

# implicit rules
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

%.o : %.F90
	$(FC) $(FCFLAGS) -c $< -o $@ -I$(ELPA_INC)

# C++ has so many strange pitfalls. But we live and learn
main.o : main.cpp elpa_interface.hpp test_functions.hpp

