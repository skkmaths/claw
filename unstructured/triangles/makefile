# Variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lgmsh -L/usr/local/lib
INCLUDES = -I/usr/local/include
MESH_TARGET = mesh.msh
# Target executable
TARGET = main

# Source files
SRCS = main.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Default target
all: $(TARGET)

# Rule to generate the .geo file (if needed) and create the .msh file
$(MESH_TARGET): mesh.geo
	gmsh -2 mesh.geo -o $(MESH_TARGET)
# Link the target executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Compile source files to object files
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean target to remove generated files
clean:
	rm -f $(TARGET) $(OBJS)

# Phony targets
.PHONY: all clean
