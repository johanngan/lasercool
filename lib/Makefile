CC = g++
AR = ar
ARFLAGS = rvs
CFLAGS = -c -std=c++11 -O2 -Wall -Wextra
SRC = read_config.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = libreadcfg.a

$(TARGET): $(OBJ)
	$(AR) $(ARFLAGS) $@ $(OBJ)

$(OBJ): $(SRC)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJ) $(TARGET)