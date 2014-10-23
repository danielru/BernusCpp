CXX=g++
FLAGS=-Winline -O3 -std=c++11

all: probe_bernus_functions.out probe_bernus.out test_bernus_steadystate.out

probe_bernus_functions.out: bernus_functions.hpp probe_bernus_functions.cpp
	$(CXX) $(FLAGS) probe_bernus_functions.cpp -o probe_bernus_functions.out

probe_bernus.out: probe_bernus.cpp bernus.hpp bernus.cpp bernus_functions.hpp
	$(CXX) $(FLAGS) probe_bernus.cpp -o probe_bernus.out

test_bernus_steadystate.out: test_bernus_steadystate.cpp bernus.cpp bernus.hpp bernus_functions.hpp
	$(CXX) $(FLAGS) test_bernus_steadystate.cpp -o test_bernus_steadystate.out

clean:
	rm -f *.out
	rm -f *.o