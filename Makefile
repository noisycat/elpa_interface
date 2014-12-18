ELPA_LIB=""
# main target
main : main.o elpa_interface.o
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(ELPA_LIB)

clean:
	rm *.o
	rm main

# implicit rules
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@ 

# C++ has so many strange pitfalls. But we live and learn
main.o : main.cpp elpa_interface.hpp

