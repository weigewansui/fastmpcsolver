#include <stdio.h>
#include <armadillo>

using namespace std;
using namespace arma;

int main() {

	mat A = randu<mat>(5,5);


	timer.tic();

	mat B = inv(A);

	cout << "time taken = " << timer.toc() / double(N) << endl;


}