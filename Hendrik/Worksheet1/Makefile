CC=g++
CFLAGS=-Wall -Wshadow
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = worksheet1.o

all: $(OBJ) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o worksheet1

%.o: %.cpp
	$(CC) $(CFLAGS) -g -c $<

doc:
	doxygen doxygenConfig.cfg

clean:
	rm -f *.o worksheet1 

clean_all:
	rm -f *.o worksheet1
	rm -rf html/
	rm -rf latex/
