CXX = g++
CXXFLAGS = -g3 -Wall -llbfgs

all: geo geo2

clean: 
	rm -f geo main.o
	rm -f geo2 maintd3.o

geo: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

geo2: maintd3.o
	$(CXX) $(CXXFLAGS) -o $@ $^

maintd3.o: main_td3.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<