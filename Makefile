CXX=g++
FLAGS=-Winline -O3 -std=c++11

all: test_bernus.out probe_bernus_functions.out

test_bernus.out: test_bernus.cpp bernus.hpp bernus.cpp Iionmodel.hpp bernus_functions.hpp
	$(CXX) $(FLAGS) test_bernus.cpp -o test_bernus.out

probe_bernus_functions.out: bernus_functions.hpp probe_bernus_functions.cpp
	$(CXX) $(FLAGS) probe_bernus_functions.cpp -o probe_bernus_functions.out
clean:
	rm -f *.out
	rm -f *.o