SHELL = /bin/sh
CC = g++
CFLAGS = -std=c++14 -O3 -flto -Wall -Wextra
LD = g++
LFLAGS = -O3 -flto

prefix = ..
bindir = $(prefix)/bin
includedir = $(prefix)/include
libdir = $(prefix)/lib

SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXECS = $(SRCS:.cpp=)

all: $(EXECS)

$(EXECS): %: %.o
	$(LD) $(LFLAGS) $< -L$(libdir) -lreadcfg -o $@

test_config.o: test_config.cpp $(libdir)/libreadcfg.a
	$(CC) -c $(CFLAGS) -I$(includedir) $< -o $@

test_timestepping.o: test_timestepping.cpp $(includedir)/lasercool/timestepping.hpp
	$(CC) -c $(CFLAGS) -I$(includedir) $< -o $@

clean:
	rm -rf $(OBJS) $(EXECS)
