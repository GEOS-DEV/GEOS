#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "eos.h"
#include "../global/misc.h"

namespace srk_par {
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
							}; // binary interaction coefficients for SRK
	std::vector<double> S2{ -0.201789, -0.004474, -0.011016, 0.010699, -0.012223, -0.012416, -0.003791, 0.003010, -0.000636, -0.007459, -0.003031 }; // SRK alpha parameter
    int nc = Tc.size();
}

using namespace std;

SRK::SRK(std::vector<int> ci_, int phase_index) {
	ci = ci_; NC = ci.size();
    phaseIndex = phase_index;
	ai = std::vector<double>(NC);
	aij = std::vector<double>(NC*NC);
	bi = std::vector<double>(NC);

}

void SRK::init(double p, double T) {p = p;
	// calculate x-independent part of parameters
    std::vector<double> alpha(NC);
	for (int i = 0; i < NC; i++) {
		// attraction parameter
        ai[i] = 0.42748*pow(srk_par::R*srk_par::Tc[ci[i]], 2)/srk_par::Pc[ci[i]];
		double S1;
		if (ci[i] == 0) { S1 = 1.244; }
		else { S1 = 0.48508 + 1.55171*srk_par::ac[ci[i]] - 0.15613*pow(srk_par::ac[ci[i]], 2); }
        alpha[i] = pow(1 + S1*(1-sqrt(T/srk_par::Tc[ci[i]])) + srk_par::S2[ci[i]]*(1-sqrt(T/srk_par::Tc[ci[i]]))/sqrt(T/srk_par::Tc[ci[i]]), 2);
	}
	for (int i = 0; i < NC; i++) {
		// interaction parameter
		for (int j = 0; j < NC; j++) {
            aij[j * NC + i] = sqrt(ai[i]*alpha[i]*ai[j]*alpha[j])*(1-(srk_par::kij[srk_par::nc*ci[i] + ci[j]]));
		}

		// repulsion parameter
        bi[i] = 0.086640*srk_par::R*srk_par::Tc[ci[i]]/srk_par::Pc[ci[i]];
	}
	return;
}

double SRK::Z(double p, double T, std::vector<double> x) {
	// Calculate a, b, A and B parameters
	a = 0;
	b = 0;
	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < NC; j++) {
			a += x[i] * x[j] * aij[j*NC + i];
		}
		b += x[i] * bi[i];
	}
	A = a * p / (pow(srk_par::R, 2) * pow(T, 2));
	B = b * p / (srk_par::R * T);

	// find roots of cubic eos
	std::vector<std::complex<double>> Z = cubic_roots(1, - 1, A-B-pow(B, 2), -A*B);
	double z;
	
	if ((fabs(Z[1].imag())>1e-34) && (fabs(Z[2].imag())>1e-34)) { // only one real root
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

std::vector<double> SRK::fugacityCoefficient(double p, double T, std::vector<double> x) {
	double z = Z(p, T, x);

	std::vector<double> phi_ij(NC);
	for (int k = 0; k < NC; k++) {
		double sumxk{ 0 };
		for (int i = 0; i < NC; i++) {
			sumxk += x[i] * aij[k * NC + i];
		}
        phi_ij[k] = exp((z-1)*bi[k]/b - log(z-B) - A/B*(2*sumxk/a - bi[k]/b)*log(1 + B/z));
	}
	return phi_ij;
}
