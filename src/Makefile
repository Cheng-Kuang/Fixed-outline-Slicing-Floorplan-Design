# Makefile

CXX      := g++
CXXFLAGS := -std=c++11 -O3

SRCS  := hw3.cpp
OBJS  := $(SRCS:.cpp=.o)
BIN_DIR := ../bin
TARGET  := $(BIN_DIR)/hw3

.PHONY: all clean

all: $(TARGET)

# Link into ../bin/hw3
$(TARGET): $(OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile .cpp ¡÷ .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
