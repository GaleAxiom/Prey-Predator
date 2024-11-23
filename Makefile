# Compiler
CXX = clang++
CXXFLAGS = -std=c++20 -Xpreprocessor -fopenmp -lomp -O3

# Source files
SRCS = solver.cpp
OBJS = $(SRCS:.cpp=.o)

# Target executable
TARGET = solver

# Build rules
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	rm -f $(OBJS) $(TARGET)