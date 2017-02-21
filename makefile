#-----Macros---------------------------------

EXEC = volce3

INCLUDE_FLAGS = -Iusr/include -Isrc -Iz3-master/include
LIB_FLAGS = z3-master/lib/libz3.so -lglpk -larmadillo -llapack -lblas -lgfortran -lm -lz -lgmp -ldl

# set up compiler and options
CXX = g++
CXXFLAGS = -g $(INCLUDE_FLAGS) -O3 -std=c++11 -Wall

#-----File Dependencies----------------------

SRC = src/main.cpp src/parser.cpp src/error.cpp src/mk.cpp src/ineq.cpp src/solver.cpp src/print.cpp \
		src/vol.cpp src/polyvest.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

all: main

depend: 
	$(CXX) $(CXXFLAGS) -MM $(SRC) > .depend.tmp
	@rm -f .depend
	@sed "s/^/src\/&/g" .depend.tmp >> .depend
	@rm -f .depend.tmp
	
-include .depend 

main: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJ) $(LIB_FLAGS)

clean:
	rm -f $(OBJ) .depend $(EXEC) 
