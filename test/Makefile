# detect OS architecture and add flags
CXX     ?= g++
CXXFLAGS = -std=c++11 -Wall -O2
# https://github.com/samtools/htslib/issues/1176#issuecomment-722033848
LDFLAGS += ../src/libseqlib.a -L../htslib -lhts
LIBS     = -lz -lm -lbz2 -llzma -lcurl -lpthread
INC      = -I. -I.. -I../htslib

# we need to find libhts at runtime
OS      := $(shell uname -s)
ifeq ($(OS),Linux)
	LDFLAGS  += -s -Wl,-rpath,../htslib -lhts
endif

OBJS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))

# don't output test-main.bin
BINS = $(patsubst %.cpp, %.bin, $(filter-out test-main.cpp, $(wildcard *.cpp)))
$(info "compiling tests $(BINS)" )

.PHONY: all test clean

all: $(BINS) $(OBJS)

test: $(BINS) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $< ${INC}

%.bin: %.o test-main.o
	${CXX} ${CXXFLAGS} -o $@ $< test-main.o $(LDFLAGS) $(LIBS)

clean:
	rm -f *.o *.bin
