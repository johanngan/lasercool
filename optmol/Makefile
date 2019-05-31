CC = g++
INCS = -I../lib -I../vendor/pcg-cpp-0.98/include
LPATHS = -L../lib
LIBS = -lm -lreadcfg
CFLAGS = -c -std=c++11 -O3 -Wall -Wextra
LD = g++
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = optical_molasses

$(EXEC): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LPATHS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS) $(EXEC)
