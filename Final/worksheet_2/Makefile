CC=g++
CFLAGS=-Wall -Wshadow -std=c++11
LDFLAGS=-lgsl -lgslcblas -lm

OBJ = main.o random_functions.o simulation_functions.o integration_functions.o

#-I /usr/local/lib/python2.7/site-packages/numpy/core/include -framework Python -I /#System/Library/Frameworks/Python.framework/Headers  -lpython2.7
#../header_files/matplotlibcpp.hpp
all: $(OBJ) doc
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS)  -o myProg

main.o:main.cpp
	$(CC) $(CFLAGS) -g -c main.cpp ../header_files/random_functions.hpp ../header_files/simulation_functions.hpp ../header_files/integration_functions.hpp

random_functions.o:../header_files/random_functions.cpp
	$(CC) $(CFLAGS) -g -c ../header_files/random_functions.cpp ../header_files/random_functions.hpp
simulation_functions.o:../header_files/simulation_functions.cpp
	$(CC) $(CFLAGS) -g -c ../header_files/simulation_functions.cpp ../header_files/simulation_functions.hpp
integration_functions.o:../header_files/integration_functions.cpp
	$(CC) $(CFLAGS) -g -c ../header_files/integration_functions.cpp ../header_files/integration_functions.hpp


doc:
	doxygen doxygenConfig.cfg


clean:
	rm -f *.o myProg

clean_all:
	rm -f *.o myProg
	rm -rf html/
	rm -rf latex/
