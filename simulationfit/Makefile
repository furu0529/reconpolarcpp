CC := g++
CFLAGS := -lfftw3 -lm -O3 -lintl -lnlopt -std=c++14 

main: main.cpp	
	$(CC) $(CFLAGS) -o main func.cpp main.cpp praxis.cpp

test:test.cpp
	$(CC) $(CFLAGS) -o test func.cpp test.cpp 

testgrid:testgrid.cpp
	$(CC) -g3 -o testgrid func.cpp testgrid.cpp -std=c++11 -lfftw3 -lm -W -Wall

testf:test.cpp
	$(CC) -O3 -o testf func.cpp test.cpp -std=c++11 -lfftw3 -lm

sub: func.cpp
	$(CC) -o func func.cpp -std=c++11 -std=c++1z

fit: fit.cpp func.hpp global.hpp	
	$(CC) $(CFLAGS) -o fit func.cpp fit.cpp praxis.cpp 
clean:
	$(RM) main func test testf testgrid fit *.so
