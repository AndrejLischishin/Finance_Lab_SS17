CC=g++
CFLAGS=-Wall -Wshadow -std=c++11
LDFLAGS=-lgsl -lgslcblas -lm

SRC = main.cpp
OBJ = main.o
DEPS = ../header_files/random_functions.hpp

all: $(OBJ) $(DEPS) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o worksheet1

main.o: main.cpp
	$(CC) $(CFLAGS) -g -c main.cpp $(DEPS)

doc:
	doxygen doxygenConfig.cfg

clean:
	rm -f *.o worksheet1

clean_all:
	rm -f *.o worksheet1
	rm -rf html/
	rm -rf latex/

