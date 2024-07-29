# Variables
CXX = g++
CXXFLAGS = -std=c++11 -Wall
LDFLAGS = -lgmsh  -L/opt/homebrew/lib
INCLUDES = -I/opt/homebrew/include

# Target executable
TARGET = main

# Source files
SRCS = main.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Default target
all: $(TARGET)

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
