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

TARGET = run

all: $(TARGET)

run: $(OBJ)
	$(CXX) -o $(TARGET) $(OBJ) $(LDFLAGS)

%.o: %.cc $(HDR)
	$(CXX) $(CFLAGS) -c $<

clean:
	rm -f *.o $(TARGET)

allclean:
	make clean
	rm -f *.vtk

