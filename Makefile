CXX=g++
FLAGS=-Winline -O3 -std=c++11

all: 
	$(CXX) $(FLAGS) test_bernus.cpp -o test_bernus.out

clean:
	rm -f *.out
	rm -f *.o