CC = g++
CFLAGS = -c -std=c++14 -O3 -flto -Wall -Wextra
INCS = -I../lib
LPATHS = -L../lib
LIBS = -lm -lreadcfg
LD = g++
LFLAGS = -O3 -flto
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
TARGET = twolevel

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(OBJS) $(LPATHS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS) $(TARGET)
