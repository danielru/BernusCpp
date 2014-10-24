CXX=clang++
FLAGS=-Winline -O3 -std=c++11
INC=-Iinclude -Isrc

all: integrate_bernus.out

integrate_bernus.out: src/integrate_bernus.cpp  src/bernus.cpp include/bernus.hpp include/bernus_functions.hpp include/Iionmodel.hpp include/IionmodelFactory.hpp
	$(CXX) $(FLAGS) src/integrate_bernus.cpp -o integrate_bernus.out $(INC)


clean:
	rm -f *.out
	rm -f *.o