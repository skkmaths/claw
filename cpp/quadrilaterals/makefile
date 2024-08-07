# Compiler and Linker
CXX = g++
CXXFLAGS = -std=c++11 -Wall
INCLUDES = -I/opt/homebrew/include $(shell pkg-config --cflags vtk)
LDFLAGS = -L/opt/homebrew/lib $(shell pkg-config --libs vtk) -lgmsh
# Target executable
TARGET = run

# Source files
SRCS = main.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Mesh files
GEO_FILE = mesh.geo
MSH_FILE = mesh.msh

# Default target
all: generate $(TARGET)

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
.PHONY: all clean generate
