CXX = g++
CXXFLAGS = -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 

Prefrules: Prefrules_init.o Prefrules_coefficient.o Prefrules_functions.o Prefrules_main.o
	$(CXX)  Prefrules_init.o Prefrules_coefficient.o Prefrules_functions.o Prefrules_main.o $(CXXFLAGS) -o Prefrules 

Prefrules_init.o: Prefrules_init.cpp Prefrules_init.h
	$(CXX) -c Prefrules_init.cpp $(CXXFLAGS)

Prefrules_coefficient.o:	Prefrules_coefficient.cpp Prefrules_coefficient.h
	$(CXX) -c Prefrules_coefficient.cpp $(CXXFLAGS)

Prefrules_functions.o: Prefrules_functions.cpp Prefrules_init.h Prefrules_functions.h Prefrules_coefficient.h
	$(CXX) -c Prefrules_functions.cpp $(CXXFLAGS) 

Prefrules_main.o: Prefrules_main.cpp Prefrules_init.h Prefrules_functions.h
	$(CXX) -c Prefrules_main.cpp $(CXXFLAGS) 
