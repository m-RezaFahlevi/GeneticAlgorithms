#include <bits/stdc++.h>
#include <armadillo>
#define N 5
using namespace std;

const arma::mat A = arma::mat(N, N).randn();
const arma::vec b = arma::vec(N).randu();

arma::vec objective(arma::vec x_apprx) {
	arma::vec b_apprx = A * x_apprx;
	arma::vec x_vec = b - b_apprx;
	return x_vec.t() * x_vec;
}

int main(void) {
	arma::vec egsol = arma::vec(N).randu();
	cout << "Ax = b" << endl;
	A.print("A :");
	b.print("b :");
	egsol.print("egsol: ");
	objective(egsol).print("Objective b - b_apprx :");

	for (int i = 0; i < 1000; ++i) {
		egsol = arma::vec(N).randu();
		objective(egsol).t().print("Objective b - b_apprx :");
	}
}
