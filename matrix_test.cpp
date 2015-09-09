#include <armadillo>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace arma;

int main() {

	mat A;
	A << 1 << 2 << 3 <<endr
	  << 4 << 5 << 6 <<endr;

	mat B;

	B << 7 <<endr
	  << 8 <<endr;

	A.insert_cols(3,B);

	A.print("A");
}