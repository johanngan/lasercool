CC = g++
INCS = -I..
LPATHS = -L..
LIBS = -lreadcfg
CFLAGS = -c -std=c++11 -O2 -Wall -Wextra
LD = g++
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = test_config

$(EXEC): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LPATHS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS)
