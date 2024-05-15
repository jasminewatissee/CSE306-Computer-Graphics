CXX = g++
CXXFLAGS = -g3 -Wall

all: geo

clean: 
	rm -f geo main.o

geo: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<