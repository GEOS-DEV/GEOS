#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "eos.h"
#include "../global/misc.h"

namespace pr_par {
	double R{ 8.3145E-5 };
	std::vector<double> Tc{ 647.14,  304.10,  126.20,  373.53, 190.58,  305.32,  369.83,  425.12, 469.70, 507.60, 540.20 }; // critical temperature
	std::vector<double> Pc{ 220.50,  73.75,   34.00,   89.63,  46.04,   48.721,  42.481,  37.960, 33.701, 30.251, 27.40 };  // critical pressure
	std::vector<double> ac{ 0.328,   0.239,   0.0377,  0.0942, 0.012,   0.0995,  0.1523,  0.2002, 0.2515, 0.3013, 0.3495 }; // acentric factor
	std::vector<double> kij{ 0, 	 0.19014, 0.32547, 0.105,  0.47893, 0.5975,	 0.5612,  0.5569, 0.5260, 0.4969, 0.4880,   // H2O
                            0.19014, 0,		  -0.0462, 0.1093, 0.0936,  0.1320,  0.1300,  0.1336, 0.1454, 0.1167, 0.1209,   // CO2
                            0.32547, -0.0462, 0, 	   0.1475, 0.0291,  0.0082,  0.0862,  0.0596, 0.0917, 0.1552, 0.1206,   // N2
							0.105, 	 0.1093,  0.1475,  0,	   0.0912,  0.0846,  0.0874,  0.0564, 0.0655, 0.0465, 0.0191,   // H2S
							0.47893, 0.0936,  0.0291,  0.0912, 0, 		0.00518, 0.01008, 0.0152, 0.0193, 0.0258, 0.0148,   // C1
							0.5975,  0.1320,  0.0082,  0.0846, 0.00518, 0, 		 0,		  0, 	  0, 	  0,      0, 	    // C2
							0.5612,  0.1300,  0.0862,  0.0874, 0.01008, 0, 		 0, 	  0, 	  0,	  0, 	  0,	    // C3
							0.5569,  0.1336,  0.0596,  0.0564, 0.0152,  0,  	 0, 	  0, 	  0,	  0, 	  0,	    // nC4
							0.5260,  0.1454,  0.0917,  0.0655, 0.0193,  0, 		 0, 	  0, 	  0,	  0, 	  0,	    // nC5
							0.4969,  0.1167,  0.1552,  0.0465, 0.0258,  0, 		 0, 	  0, 	  0,	  0, 	  0,	    // nC6
							0.4880,  0.1209,  0.1206,  0.0191, 0.0148,  0, 		 0, 	  0, 	  0,	  0, 	  0,	    // nC7
							}; // binary interaction coefficients for PR
	int nc = Tc.size();
}

using namespace std;

PR::PR(std::vector<int> ci_, int phase_index) { //: EoS{ NC } {
	ci = ci_; NC = ci.size();
	phaseIndex = phase_index;
	ai = std::vector<double>(NC);
	aij = std::vector<double>(NC*NC);
	bi = std::vector<double>(NC);
}

void PR::init(double p, double T) {p = p;
	// calculate x-independent part of parameters
	for (int i = 0; i < NC; i++) {
		// attraction parameter
		double kappa = 0.37464 + 1.54226 * pr_par::ac[ci[i]] - 0.26992 * pow(pr_par::ac[ci[i]], 2);
		double alpha = pow(1 + kappa * (1 - sqrt(T / pr_par::Tc[ci[i]])), 2);
		ai[i] = 0.45724 * pow(pr_par::R, 2) * pow(pr_par::Tc[ci[i]], 2) / pr_par::Pc[ci[i]] * alpha;
	}
	for (int i = 0; i < NC; i++) {
		// interaction parameter
		for (int j = 0; j < NC; j++) {
			aij[j * NC + i] = sqrt(ai[i] * ai[j]) * (1 - pr_par::kij[pr_par::nc*ci[i] + ci[j]]);
		}

		// repulsion parameter
		bi[i] = 0.0778 * pr_par::R * pr_par::Tc[ci[i]] / pr_par::Pc[ci[i]];
	}
	return;
}

double PR::Z(double p, double T, std::vector<double> x) {
	// Calculate a, b, A and B parameters
	a = 0;
	b = 0;
	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {
			a += x[i] * x[j] * aij[j*NC + i];
		}
		b += x[i] * bi[i];
	}
	A = a * p / (pow(pr_par::R, 2) * pow(T, 2));
	B = b * p / (pr_par::R * T);

	// find roots of cubic eos
	std::vector<std::complex<double>> Z = cubic_roots(1, B - 1, A - 2 * B - 3 * pow(B, 2), -1 * (A * B - pow(B, 2) - pow(B, 3)));
	double z;
	if (fabs(Z[1].imag())>1e-34 && fabs(Z[2].imag())>1e-34) { // only one real root
		z = Z[0].real();
	}
	else {
		z = Z[0].real();
		if (phaseIndex == 0) { // 0 == v
			if (Z[1].real() > z) { z = Z[1].real(); }
			if (Z[2].real() > z) { z = Z[2].real(); }
		}
		else if (phaseIndex == 1) { // 1 == l
			if (Z[1].real() < z) { z = Z[1].real(); }
			if (Z[2].real() < z) { z = Z[2].real(); }
		}
	}
	return z;
}

std::vector<double> PR::fugacityCoefficient(double p, double T, std::vector<double> x) {
	double z = Z(p, T, x);

	std::vector<double> phi_ij(NC);
	for (int k = 0; k < NC; k++) {
		double sumxk{ 0 };
		for (int i = 0; i < NC; i++) {
			sumxk += x[i] * aij[i * NC + k];
		}
		phi_ij[k] = exp(bi[k] / b * (z - 1.0) - log(z - B) - A / (2.0 * sqrt(2.0) * B) * (2.0 * sumxk / a - bi[k] / b) * log((z + 2.414 * B) / (z - 0.414 * B)));
	}

	return phi_ij;
}
