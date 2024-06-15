CXX = g++
CXXFLAGS = -g3 -Wall -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -Ofast -llbfgs -lomp

all: geo

clean: 
	rm -f geo main.o

geo: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<