# Possible compilers: g++, icpc
CXX    = g++

#
# DONT CHANGE ANYTHING BELOW
# UNLESS YOU KNOW WHAT YOU 
# ARE DOING !!!
#

# GNU compiler
ifeq ($(CXX),g++)
	CFLAGS = -O3 -Wall -fopenmp
	LDFLAGS = -fopenmp
endif

# Intel compiler
ifeq ($(CXX),icpc)
	CFLAGS = -fast -qopenmp
	LDFLAGS = -qopenmp
endif

HDR = $(wildcard *.h)
SRC = $(wildcard *.cc)
OBJ = $(patsubst %.cc,%.o,$(SRC))

TARGET = twodproblem

all: $(TARGET)

twodproblem: $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ)

$(OBJ): $(HDR)

.cc.o: $(HDR)
	$(CXX) $(CFLAGS) -c $*.cc
clean:
	rm -f *.o $(TARGET)

allclean:
	make clean
	rm -f *.vtk
