CXX = g++
CXXFLAGS = -g3 -Wall -pthread

all: template

clean: 
	rm -f template template.o

template: template.o
	$(CXX) $(CXXFLAGS) -o $@ $^

template.o: template.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<