CXX=icpc
LDFLAGS=-L${HDF5_ROOT}/lib -L${GSL_ROOT}/lib -lgsl -lgslcblas -lhdf5 -lz
CXXFLAGS=-O3 -g -pg -std=c++17 -I${HDF5_ROOT}/include -I${OPENMPI_ROOT}/include -I${GSL_ROOT}/include -Wall -Warray-bounds 

SRCS := ../src/m1_closure.cc ../src/m1_collisional.cc ../src/m1_fluxes.cc 
OBJECTS := m1_closure.o m1_collisional.o m1_fluxes.o 

all: objs test_fluxes test_implicit clean-obj

test_fluxes: $(OBJECTS) test_fluxes.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o test_fluxes

test_implicit: $(OBJECTS) test_implicit.o
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o test_implicit

%.o: %.cc
	$(CXX) $(CXXFLAGS) -DDEBUG_ROOTFINDING -DSTANDALONE -c $^ 

objs: $(SRCS)
	$(CXX) $(CXXFLAGS) -DDEBUG_ROOTFINDING -DSTANDALONE -c $(SRCS) 

clean-obj:
	rm *.o

clean:
	rm *.out
