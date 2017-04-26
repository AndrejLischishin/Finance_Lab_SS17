CC=g++
CFLAGS=-Wall -Wshadow
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = main.o

all: $(OBJ) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o myProg

%.o: %.cpp
	$(CC) $(CFLAGS) -g -c $<

doc:
	doxygen -g doxygenConfigFile.cfg

clean:
	rm -f *.o myProg

clean_all:
	rm -f *.o myProg
	rm -rf html/
	rm -rf latex/
