###################################################################
#  Makefile for homework 2
#
#  Abu Talha
#  SMU Mathematics
#  Math 6370
###################################################################

# compiler & flags
CXX = g++
CXXFLAGS = -O2 --std=c++11

# makefile targets
all: vec2d_test.exe vec2d_b_test.exe

vec2d_test.exe : vec2d_test.cpp vec2d.o GramSchmidt2d.o
	$(CXX) $(CXXFLAGS) $^  -o $@

vec2d.o : vec2d.cpp vec2d.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

vec2d_b_test.exe : vec2d_b_test.cpp vec2d_b.o GramSchmidt2d_b.o 
	$(CXX) $(CXXFLAGS) $^  -o $@

vec2d_b.o : vec2d_b.cpp vec2d_b.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

GramSchmidt2d.o : GramSchmidt2d.cpp vec2d.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

GramSchmidt2d_b.o : GramSchmidt2d_b.cpp vec2d_b.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	\rm -f *.o *.out a_data

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######
