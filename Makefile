CXX=clang++
FLAGS=-Winline -O3 -std=c++11
INC=-Iinclude

all: integrate_bernus.out

build/bernus_functions.o: src/bernus_functions.cpp include/bernus_functions.hpp
	$(CXX) $(FLAGS) -c src/bernus_functions.cpp -o build/bernus_functions.o $(INC)

build/bernus.o: src/bernus.cpp include/bernus.hpp
	$(CXX) $(FLAGS) -c src/bernus.cpp -o build/bernus.o $(INC)

integrate_bernus.out: build/bernus_functions.o build/bernus.o include/Iionmodel.hpp include/IionmodelFactory.hpp
	$(CXX) $(FLAGS) build/bernus_functions.o build/bernus.o src/integrate_bernus.cpp -o integrate_bernus.out $(INC)


clean:
	rm -f *.out
	rm -f build/*.o