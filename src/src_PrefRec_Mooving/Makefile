CXX = g++
CXXFLAGS = -mpopcnt -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 

Prefrec: Prefrec_init.o Prefrec_functions.o Prefrec_main.o
	$(CXX)  Prefrec_init.o Prefrec_functions.o Prefrec_main.o $(CXXFLAGS) -o Prefrec 

Prefrec_init.o: Prefrec_init.cpp Prefrec_init.h
	$(CXX) -c Prefrec_init.cpp $(CXXFLAGS) 
Prefrec_functions.o: Prefrec_functions.cpp Prefrec_init.h Prefrec_functions.h
	$(CXX) -c Prefrec_functions.cpp $(CXXFLAGS) 
Prefrec_main.o: Prefrec_main.cpp Prefrec_init.h Prefrec_functions.h
	$(CXX) -c Prefrec_main.cpp $(CXXFLAGS) 
