CXX=g++
CXXFLAGS=-std=c++11 -pthread -O3

OBJ=su2_x.o

all: su2_x_test

tests: test_exp_su2

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

su2_x_test: su2_x_test.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $< -o $@

test_exp_su2: test_exp_su2.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(OBJ) $< -o $@
