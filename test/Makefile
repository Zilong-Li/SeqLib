CXX ?= g++
CXXFLAGS = -std=c++11 -Wall -O2
LDFLAGS  = ../lib/libseqlib.a ../lib/libhts.a
LIBS     = -lz -lm -lbz2 -llzma
INC      = -I. -I.. -I../htslib

# don't output test-main.bin
OBJS = $(patsubst %.cpp, %.bin, $(filter-out test-main.cpp, $(wildcard *.cpp)))
$(info "compiling tests $(OBJS)" )

.PHONY: all test clean

all: $(OBJS)

test: $(OBJS)

test-main.o: test-main.cpp catch.hh
	$(CXX) $(CXXFLAGS) $< -c -o $@

%.bin: %.cpp test-main.o
	${CXX} ${CXXFLAGS} $< test-main.o -o $@ $(INC) $(LIBS) $(LDFLAGS)

clean:
	rm -f *.o *.bin
