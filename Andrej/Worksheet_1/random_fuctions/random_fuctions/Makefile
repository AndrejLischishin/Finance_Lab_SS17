CC=g++
CFLAGS=-Wall -Wshadow -std=c++11
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = main.o random_functions.o

all: $(OBJ) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o myProg

main.o:main.cpp
	$(CC) $(CFLAGS) -g -c main.cpp random_functions.hpp
reject_sample.o:random_functions.cpp
	$(CC) $(CFLAGS) -g -c random_functions.cpp random_functions.hpp


doc:
	doxygen -g doxygenConfig.cfg

clean:
	rm -f *.o myProg

clean_all:
	rm -f *.o myProg
	rm -rf html/
	rm -rf latex/
