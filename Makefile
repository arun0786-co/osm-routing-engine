CXX ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -march=native -Wall -Wextra -Wno-unused-parameter
LDLIBS = -lz -lpthread

LIB = src/pbf.o src/graph.o src/build.o src/ch.o

all: router

router: $(LIB) src/main.o src/server.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

test: tests/test_main
	./tests/test_main

tests/test_main: $(LIB) tests/test_main.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

src/%.o: src/%.cpp src/*.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

tests/%.o: tests/%.cpp src/*.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f src/*.o tests/*.o router tests/test_main

.PHONY: all clean test
