# detect OS architecture and add flags
OS      := $(shell uname -s)
CXX     ?= g++
CXXFLAGS = -std=c++11 -Wall -O2
# https://github.com/samtools/htslib/issues/1176#issuecomment-722033848
LDFLAGS  = ../src/libseqlib.a
LIBS     = -lz -lm -lbz2 -llzma -lcurl -lhts -L../htslib
INC      = -I. -I.. -I../htslib

# we need to find libhts at runtime
ifeq ($(OS),Linux)
	LDFLAGS  += -Wl,-rpath,../htslib
	CXXFLAGS += -s
endif

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
