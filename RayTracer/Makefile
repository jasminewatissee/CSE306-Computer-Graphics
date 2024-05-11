CXX = g++
CXXFLAGS = -g3 -Wall -pthread -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -lomp -Ofast

all: raytracer soft_shadow

clean: 
	rm -f raytracer main.o
	rm -f soft_shadow soft_shadow.o

raytracer: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

soft_shadow: soft_shadow.o
	$(CXX) $(CXXFLAGS) -o $@ $^

soft_shadow.o: soft_shadow.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<