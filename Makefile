CXX = g++
CXXFLAGS = -g3 -Wall -pthread -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -lomp -Ofast

all: template

clean: 
	rm -f template template.o

template: template.o
	$(CXX) $(CXXFLAGS) -o $@ $^

template.o: template.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<