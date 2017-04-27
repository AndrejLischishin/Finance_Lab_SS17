CC=g++
CFLAGS=-Wall -Wshadow
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = main.o

all: $(OBJ) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o worksheet1

%.o: %.cpp
	$(CC) $(CFLAGS) -g -c $<

doc:
	doxygen doxygenConfig.cfg

clean:
	rm -f *.o main

clean_all:
	rm -f *.o main
	rm -rf html/
	rm -rf latex/
