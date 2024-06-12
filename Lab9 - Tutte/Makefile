CXX = g++
CXXFLAGS = -g3 -Wall -pthread

all: main

clean: 
	rm -f main main.o

main: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<