CC=g++
CFLAGS=-Wall -Wshadow -std=c++11
LDFLAGS=-lgsl -lgslcblas -lm

SRC = main.cpp
OBJ = main.o
DEPS = ../header_files/random_functions.hpp

all: $(OBJ) $(DEPS)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o myProg

main.o: main.cpp
	$(CC) $(CFLAGS) -g -c main.cpp $(DEPS)

#doc:
#	doxygen -g doxygenConfig.cfg

clean:
	rm -f *.o myProg

clean_all:
	rm -f *.o myProg
	rm -rf html/
	rm -rf latex/
