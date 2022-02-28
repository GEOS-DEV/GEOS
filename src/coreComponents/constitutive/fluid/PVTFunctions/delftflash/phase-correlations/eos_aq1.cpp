#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "eos.h"
#include "../global/misc.h"

using namespace std;

namespace aq1_par { // parameters for Ziabaksh correlation for Aq phase
	// std::vector<std::string> components{ "H2O", "CO2", "N2", "H2S", "CH4" };
	std::vector<double> labda{0, -0.0652869, -2.0939363, 1.03658689, -5.7066455E-1, -2.143686, 0.513068, 0.52862384, 0.52862384, 0.52862384, 0.52862384, // c1
							0, 1.6790636E-4, 3.1445269E-3, -1.1784797E-3, 7.2997588E-4, 2.598765E-3, -0.000958, -1.0298104E-3, -1.0298104E-3, -1.0298104E-3, -1.0298104E-3,  // c2
							0, 40.838951, 3.913916E2, -1.7754826E2, 1.5176903E2, 4.6942351E2, 0, 0, 0, 0, 0,        // c3
							0, 0, -2.9973977E-7, -4.5313285E-4, 3.1927112E-5, -4.6849541E-5, 0, 0, 0, 0, 0,           // c4
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,											    // c5	
							0, -3.9266518E-2, -1.5918098E-5, 0, -1.642651E-5, 0, 0, 0, 0, 0, 0, 		    // c6
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,											   // c7
							0, 2.1157167E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 								   // c8
							0, 6.5486487E-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 							   // c9
							0, 0, 0, 0.4775165E2, 0, 0, 0, 0, 0, 0, 0,								   // c10
							0, 0, 0, 0, 0, -8.4616602E-10, 0, 0, 0, 0, 0,								   // c11
							0, 0, 0, 0, 0, 1.095219E-6, 0, 0, 0, 0, 0 };								   // c12
	std::vector<double> ksi{0, -1.144624E-2, -6.3981858E-3, 0.010274152, -2.9990084E-3, -1.0165947E-2, -0.007485, 0.0206946, 0.0206946, 0.0206946, 0.0206946, // c1
							0, 2.8274958E-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 								   // c2
							0, 1.3980876E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 								   // c6
							0, -1.4349005E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0 };								   // c8
	std::vector<double> par{0, -0.114535, -0.008194, 0.77357854, -0.092248, -0.6091, -1.1471, -1.6849, -1.6849, -1.6849, -1.6849, // eta
							0, -5.279063, -5.175337, 0.27049433, -5.779280, -16.8037, -25.3879, -33.8492, -33.8492, -33.8492, -33.8492, // tau
							0, 6.187967, 6.906469, 0.27543436, 7.262730, 20.0628, 28.2616, 36.1457, 36.1457, 36.1457, 36.1457, // beta
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};  										   // Gamma
	double R{ 83.145 };
	double Tcw{ 647.14 };
	double Pcw{ 220.50 };
	int nc = labda.size()/12;
}

AQ1::AQ1(std::vector<int> ci_) { //: EoS{ NC } {
	ci = ci_; NC = ci.size();
	k_H = std::vector<double>(NC);
	labda_i = std::vector<double>(NC);
	ksi_i = std::vector<double>(NC);
}

void AQ1::init(double p, double T) {
	// Calculate composition-independent parameters of Aq-EoS
	double tau = 1 - T / aq1_par::Tcw;
	double P_s = aq1_par::Pcw * exp(aq1_par::Tcw / T * (tau * -7.85951783 + 1.84408259 * pow(tau, 1.5) - 11.7866497 * pow(tau, 3)
			+ 22.6807411 * pow(tau, 3.5) - 15.9618719 * pow(tau, 4.0) + 1.80122502 * pow(tau, 7.5)));
	double tc = T - 273.15;  // temperature in Celsius
	double V0 = (1 + 18.159725E-3 * tc) / (0.9998396 + 18.224944E-3 * tc - 7.922210E-6 * pow(tc, 2) - 55.44846E-9 * pow(tc, 3) + 149.7562E-12 * pow(tc, 4) - 393.2952E-15 * pow(tc, 5));
	double B = 19654.320 + 147.037 * tc - 2.21554 * pow(tc, 2) + 1.0478E-2 * pow(tc, 3) - 2.2789E-5 * pow(tc, 4);
	double A1 = 3.2891 - 2.3910E-3 * tc + 2.8446E-4 * pow(tc, 2) - 2.8200E-6 * pow(tc, 3) + 8.477E-9 * pow(tc, 4);
	double A2 = 6.245E-5 - 3.913E-6 * tc - 3.499E-8 * pow(tc, 2) + 7.942E-10 * pow(tc, 3) - 3.299E-12 * pow(tc, 4);
	double V = V0 - V0 * p / (B + A1 * p + A2 * pow(p, 2));  // volume of pure water at p[cm3 / g]
	double M = 18.0152;
	double f0_H2O = P_s * exp((p - P_s) * M * V / (aq1_par::R * T));  // fugacity of pure water[bar]
	double rho0_H2O = 1 / V;  // density of pure water at p[g / cm3]
	K0_H2O = pow(10.0, (-2.209 + 3.097E-2 * tc - 1.098E-4 * pow(tc, 2) + 2.048E-7 * pow(tc, 3)));  // equilibrium constant for H2O at 1 bar

	// Construct constant part of fugacity coefficients
	for (int i = 0; i < NC; i++) {
		if (ci[i] != 0) { // non-H2O components
			double dB = aq1_par::par[1*aq1_par::nc + ci[i]] + aq1_par::par[3*aq1_par::nc + ci[i]] * p + aq1_par::par[2*aq1_par::nc + ci[i]] * sqrt(1000 / T);
			k_H[i] = exp((1 - aq1_par::par[0*aq1_par::nc + ci[i]]) * log(f0_H2O) + aq1_par::par[0*aq1_par::nc + ci[i]] * log(aq1_par::R * T / M * rho0_H2O) + 2 * rho0_H2O * dB);
			labda_i[i] = aq1_par::labda[0*aq1_par::nc + ci[i]] + aq1_par::labda[1*aq1_par::nc + ci[i]] * T + aq1_par::labda[2*aq1_par::nc + ci[i]] / T + aq1_par::labda[3*aq1_par::nc + ci[i]] * p + aq1_par::labda[4*aq1_par::nc + ci[i]] / p + \
				aq1_par::labda[5*aq1_par::nc + ci[i]] * p / T + aq1_par::labda[6*aq1_par::nc + ci[i]] * T / pow(p, 2) + aq1_par::labda[7*aq1_par::nc + ci[i]] * p / (630 - T) + aq1_par::labda[8*aq1_par::nc + ci[i]] * T * log(p) \
				+ aq1_par::labda[9*aq1_par::nc + ci[i]] * p / pow(T, 2) + aq1_par::labda[10*aq1_par::nc + ci[i]] * pow(p, 2) * T + aq1_par::labda[11*aq1_par::nc + ci[i]] * p * T;
			ksi_i[i] = aq1_par::ksi[0*aq1_par::nc + ci[i]] + aq1_par::ksi[1*aq1_par::nc + ci[i]] * T + aq1_par::ksi[2*aq1_par::nc + ci[i]] * p / T + aq1_par::ksi[3*aq1_par::nc + ci[i]] * p / (630.0 - T);
		}
	}
	return;
}

std::vector<double> AQ1::fugacityCoefficient(double p, double T, std::vector<double> x) {
	// Construct fugacity coefficients
	x=x;
	std::vector<double> phi_ij(NC);
	double molality{ 0 };
	for (int i = 0; i < NC; i++) {
		if (ci[i] == 0) { //0 == H2O
			phi_ij[i] = K0_H2O * exp((p - 1.0) * V_H2O / (aq1_par::R * T)) / p;
		}
		else {
			double gamma = exp((2 * molality * labda_i[i]) + (molality * molality * ksi_i[i]));
			phi_ij[i] = k_H[i] * gamma / p;
		}
	}
	return phi_ij;
}