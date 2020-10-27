main: main.cpp	
	g++ -o main func.cpp main.cpp praxis.cpp -lfftw3 -lm -W -Wall -O3 -std=c++14
#	g++ -o main main.cpp -std=c++11 -I/Users/kokoorikunihiko/Downloads/matplotlib-cpp -I/usr//local/include/python3.8/ -I/usr//local/lib/python3.8/site-packages/numpy/core/include/ -lpython3.8 
#-std=c++11

pymain: func.cpp
	g++ -fPIC -shared -I/usr/local/include/python3.8 -I/usr/local/include/boost/python/ -L/usr/local/Cellar/python@3.8/3.8.6/Frameworks/Python.framework/Versions/3.8/lib/ -lboost_python38 func.cpp -o optim.so -std=c++14 -lfftw3 -lm -lpython3.8 -lintl -lboost_numpy38

python: intropython.cpp
	g++ -DPIC -shared -fPIC -I/usr/local/include/python3.8 -I/usr/local/include/boost/python/ -L/usr/local/Cellar/python@3.8/3.8.6/Frameworks/Python.framework/Versions/3.8/lib/ -lboost_python38 intropython.cpp -o MyClassModule.so -std=c++14 -lfftw3 -lm -lpython3.8 -lintl
test:test.cpp
	g++ -o test func.cpp test.cpp -std=c++11 -lfftw3 -lm -W -Wall -I/usr/local/include/python3.8 -I/usr/local/include/boost/python/ -L/usr/local/Cellar/python@3.8/3.8.6/Frameworks/Python.framework/Versions/3.8/lib/ -lboost_python38 -lpython3.8 -lintl -lboost_numpy38

testgrid:testgrid.cpp
	g++ -g3 -o testgrid func.cpp testgrid.cpp -std=c++11 -lfftw3 -lm -W -Wall

testf:test.cpp
	g++ -O3 -o testf func.cpp test.cpp -std=c++11 -lfftw3 -lm

sub: func.cpp
	g++ -o func func.cpp -std=c++11 -std=c++1z

example: example.cpp
	g++ example.cpp -std=c++17 -I/Users/kokoorikunihiko/Desktop/reserch/juniorreserch/sumiyakun/simulationnow/ceres/ceres-solver/include -I/Users/kokoorikunihiko/Desktop/reserch/juniorreserch/sumiyakun/simulationnow/ceres/ceres-solver/internal -I/Users/kokoorikunihiko/Desktop/reserch/juniorreserch/sumiyakun/simulationnow/ceres/ceres-bin -L/Users/kokoorikunihiko/Desktop/reserch/juniorreserch/sumiyakun/simulationnow/ceres/ceres-bin/lib -lceres -lgflags -lglog -lcxsparse -lcamd -lamd -lccolamd -lcolamd -lcholmod -lblas -llapack -lsuitesparseconfig

fit: fit.cpp	
	g++ -o fit func.cpp fit.cpp praxis.cpp -lfftw3 -lm -O3 -std=c++14 -I/usr/local/include/python3.8 -I/usr/local/include/boost/python/ -L/usr/local/Cellar/python@3.8/3.8.6/Frameworks/Python.framework/Versions/3.8/lib/ -lboost_python38 -lpython3.8 -lintl -lboost_numpy38 -lnlopt
clean:
	$(RM) main func test testf testgrid fit *.so
