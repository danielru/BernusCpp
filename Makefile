CXX=clang++
FLAGS=-Winline -O3 -std=c++11 -Wfatal-errors -g
INC=-Iinclude

all: integrate_bernus.out

build/bernus_functions.o: build src/bernus_functions.C include/bernus_functions.h
	$(CXX) $(FLAGS) -c src/bernus_functions.C -o build/bernus_functions.o $(INC)

build/bernus.o: build src/bernus.C include/bernus.h
	$(CXX) $(FLAGS) -c src/bernus.C -o build/bernus.o $(INC)

integrate_bernus.out: build build/bernus_functions.o build/bernus.o include/Iionmodel.h include/IionmodelFactory.h src/integrate_bernus.C
	$(CXX) $(FLAGS) build/bernus_functions.o build/bernus.o src/integrate_bernus.C -o integrate_bernus.out $(INC)

build:
	mkdir build

clean:
	rm -f *.out
	rm -f build/*.o