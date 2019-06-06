CC = g++
INCS = -I..
LPATHS = -L..
LIBS = -lreadcfg
CFLAGS = -c -std=c++14 -O3 -Wall -Wextra
LD = g++
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXECS = $(SRCS:.cpp=)

all: $(EXECS)

$(EXECS): %: %.o
	$(LD) $(LFLAGS) $< $(LPATHS) $(LIBS) -o $@

test_config.o: test_config.cpp
	$(CC) $(CFLAGS) $(INCS) $< -o $@

test_timestepping.o: test_timestepping.cpp ../timestepping.hpp
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS) $(EXECS)
