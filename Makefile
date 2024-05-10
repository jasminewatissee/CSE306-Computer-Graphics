CXX = g++
CXXFLAGS = -g3 -Wall -pthread -L/usr/local/opt/libomp/lib -I/usr/local/opt/libomp/include -Xpreprocessor -fopenmp -lomp -Ofast

all: template

clean: 
	rm -f template template.o

template: template.o
	$(CXX) $(CXXFLAGS) -o $@ $^

template.o: template.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

soft_shadow: template_soft_shadow.o
	$(CXX) $(CXXFLAGS) -o $@ $^

template_soft_shadow.o: template_soft_shadow.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<