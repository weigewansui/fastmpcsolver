mpc_test:
	g++ mpc_test.cpp -o mpctest -std=c++11 -larmadillo

mpc_debug:
	g++ -g mpc_test.cpp -o mpctest -std=c++11 -larmadillo

clean:
	rm -f mpctest mpc_debug mpc_test