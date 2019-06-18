CC = g++
CFLAGS = -c -std=c++14 -O3 -flto -Wall -Wextra
INCS = -I../lib
LPATHS = -L../lib
LIBS = -lm -lreadcfg
LD = g++
LFLAGS = -O3 -flto
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)
TARGETS = swapint swapmotion

all: $(TARGETS)

swapint: swapint.o HInt.o HSwap.o
swapmotion: swapmotion.o HMotion.o HSwap.o

$(TARGETS):
	$(LD) $(LFLAGS) $^ $(LPATHS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCS) $< -o $@

clean:
	rm -rf $(OBJS) $(TARGETS)
