#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "misc.h"

#define M_PI 3.14159265358979323846 /* pi */

using namespace std;

Combinations::Combinations(int ni, int ne, int l) {
	// https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
	Ne = ne; L = l;
	m_arr = std::vector<int>(Ne);
	for (int i = 0; i < Ne; i++) { m_arr[i] = i; }  // array with component numbers for each m_i(v)
	m_data = std::vector<int>(ni*L);
	std::vector<int> temp(L, 0);
	unique_combinations(temp, 0, Ne-1, 0);
}

void Combinations::unique_combinations(std::vector<int> temp, int start, int end, int index) {
	/* std::vector<int> arr ---> Input Array
	std::vector<int> data ---> Temporary array to store current combination
	start & end ---> Starting and ending indexes in arr
	index ---> Current index in data
	r ---> Size of a combination to be printed */
	// This code is contributed by rathbhupendra
	if (index == L)	{ 
		for (int j = 0; j < L; j++) { m_data[index0 + j] = temp[j]; }
		index0 += L;
		return; 
	}

	// replace index with all possible elements. The condition "end-i+1 >= r-index" makes sure
	// that including one element at index will make a combination with remaining elements at remaining positions
	for (int i = start; i <= end && end - i + 1 >= L - index; i++) {
		temp[index] = m_arr[i];
		unique_combinations(temp, i+1, end, index+1);
	}
}

std::vector<std::complex<double>> cubic_roots(double a, double b, double c, double d) {
	// http://www.cplusplus.com/forum/beginner/234717/
	std::vector<std::complex<double>> Z(3);
	std::vector<double> real(3), imag(3);
	b /= a;
	c /= a;
	d /= a;

	double q = (3.*c - b*b) / 9.;
	double r = (b*(9.*c - 2*b*b) - 27.*d) / 54.;
	double D = q*q*q + r*r;
	double term1 = b / 3.;
	imag[0] = 0.;
	
	if (D > 0) {  // one real, two complex roots
		double s = r + sqrt(D);
		if (s < 0) { s = -pow(-s, 1./3.); }
		else { s = pow(s, 1./3.); }
		
		double t = r - sqrt(D);
		if (t < 0) { t = -pow(-t, 1./3.); }
		else { t = pow(t, 1./3.); }

		real[0] = -term1 + s + t;
		term1 += (s + t)/2.;
		real[1] = -term1; real[2] = -term1;
		term1 = sqrt(3.) * (s - t)/2.;
		imag[1] = term1; imag[2] = -term1;
	} else if (fabs(D)<1e-34) {  // all roots real, at least two equal
		imag[1] = 0.; imag[2] = 0.;
		double r13;
		if (r < 0) { r13 = -pow(-r, 1./3.); }
		else { r13 = pow(r, 1./3.); }
		real[0] = -term1 + 2.*r13;
		real[1] = -(r13 + term1); real[2] = -(r13 + term1);
	} else {  // all roots real and unequal
		imag[1] = 0.; imag[2] = 0.;
		q = -q;
		double dum1 = q*q*q;
		dum1 = acos(r/sqrt(dum1));
		double r13 = 2.*sqrt(q);
		real[0] = -term1 + r13*cos(dum1/3.);
		real[1] = -term1 + r13*cos((dum1 + 2.*M_PI) / 3.);
		real[2] = -term1 + r13*cos((dum1 + 4.*M_PI) / 3.);
	}

	for (int i = 0; i < 3; i++) { Z[i] = std::complex<double>(real[i], imag[i]); }

	return Z;
}

int factorial(int n) {
    if(n > 1) { return n * factorial(n - 1); }
    else { return 1; }
}

int searchString(std::vector<std::string> A, int size, std::string target)
{
	for (int j = 0; j < size; j++) {
		if (A[j] == target) {
			return j;
		}
	}
	cout << "Component or phase not found!" << endl;
	return 100;
}