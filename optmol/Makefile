CC = g++
CFLAGS = -c -std=c++11 -O3 -flto -Wall -Wextra
INCS = -I../lib -I../vendor/pcg-cpp-0.98/include
LPATHS = -L../lib
LIBS = -lm -lreadcfg
LD = g++
LFLAGS = -O3 -flto
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
EXEC = optical_molasses

$(EXEC): $(OBJS) mathutil.hpp
	$(LD) $(LFLAGS) $(OBJS) $(LPATHS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS) $(EXEC)
