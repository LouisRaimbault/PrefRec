
CXX = g++
CXXFLAGS = -mpopcnt -msse4 -mabm -mmmx -mavx -mavx -mavx2 -mfma -O3 
CXXFLAGSB = -O3

all:	Prefrec Prefrules

Prefrec: Prefrec.cpp Prefrec.h 
	$(CXX)  $(CXXFLAGS) -o Prefrec -cpp  Prefrec.cpp  

Prefrules:	PrefRules.cpp PrefRules.h
	$(CXX)  $(CXXFLAGSB) -o Prefrules -cpp  PrefRules.cpp









