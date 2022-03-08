#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <unordered_map>

#include "eos.h"
#include "../global/global.h"

using namespace std;

// namespace aq1_par { // parameters for Ziabaksh correlation for Aq phase
// 	std::unordered_map<std::string, std::vector<double>> labda = {
// 		{"CO2", {-0.0652869, 1.6790636E-4, 40.838951, 0, 0, -3.9266518E-2, 0, 2.1157167E-2, 6.5486487E-6, 0, 0, 0}},
// 		{"N2", {-2.0939363, 3.1445269E-3, 3.913916E2, -2.9973977E-7, 0, -1.5918098E-5, 0, 0, 0, 0, 0, 0}},
// 		{"H2S", {1.03658689, -1.1784797E-3, -1.7754826E2, -4.5313285E-4, 0, 0, 0, 0, 0, 0.4775165E2, 0, 0}},
// 		{"C1", {-5.7066455E-1, 7.2997588E-4, 1.5176903E2, 3.1927112E-5, 0, -1.642651E-5, 0, 0, 0, 0, 0, 0}},
// 		{"C2", {-2.143686, 2.598765E-3, 4.6942351E2, -4.6849541E-5, 0, 0, 0, 0, 0, 0, -8.4616602E-10, 1.095219E-6}},
// 		{"C3", {0.513068, -0.000958, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		{"iC4", {0.52862384, -1.0298104E-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		{"nC4", {0.52862384, -1.0298104E-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"iC5", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC5", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC6", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC7", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC8", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC9", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 		// {"nC10", {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}},
// 	};
// 	std::unordered_map<std::string, std::vector<double>> ksi = {
// 		{"CO2", {-1.144624E-2, 2.8274958E-5, 1.3980876E-2, -1.4349005E-2}},
// 		{"N2", {-6.3981858E-3, 0, 0, 0}},
// 		{"H2S", {0.010274152, 0, 0, 0}},
// 		{"C1", {-2.9990084E-3, 0, 0, 0}},
// 		{"C2", {-1.0165947E-2, 0, 0, 0}},
// 		{"C3", {-0.007485, 0, 0, 0}},
// 		{"iC4", {0.0206946, 0, 0, 0}},
// 		{"nC4", {0.0206946, 0, 0, 0}},
// 	};
// 	std::unordered_map<std::string, double> eta = {{"CO2", -0.114535}, {"N2", -0.008194}, {"H2S", 0.77357854}, {"C1", -0.092248}, {"C2", -0.6091}, {"C3", -1.1471}, {"iC4", -1.6849}, {"nC4", -1.6849}};
// 	std::unordered_map<std::string, double> tau = {{"CO2", -5.279063}, {"N2", -5.175337}, {"H2S", 0.27049433}, {"C1", -5.779280}, {"C2", -16.8037}, {"C3", -25.3879}, {"iC4", -33.8492}, {"nC4", -33.8492}};
// 	std::unordered_map<std::string, double> beta = {{"CO2", 6.187967}, {"N2", 6.906469}, {"H2S", 0.27543436}, {"C1", 7.262730}, {"C2", 20.0628}, {"C3", 28.2616}, {"iC4", 36.1457}, {"nC4", 36.1457}};
// 	std::unordered_map<std::string, double> Gamma = {{"CO2", 0}, {"N2", 0}, {"H2S", 0}, {"C1", 0}, {"C2", 0}, {"C3", 0}, {"iC4", 0}, {"nC4", 0}};
	
// 	double R{ 83.145 };
// }

void AQ1::parameters(double p_, double T_) {
	// Calculate composition-independent parameters of Aq-EoS
	p = p_; T = T_;
	water_index = std::distance(components.begin(), std::find(components.begin(), components.end(), "H2O"));

	double tau = 1 - T / CompProp::prop["Tc"]["H2O"];
	double P_s = CompProp::prop["Pc"]["H2O"] * exp(CompProp::prop["Tc"]["H2O"] / T * (tau * -7.85951783 + 1.84408259 * pow(tau, 1.5) - 11.7866497 * pow(tau, 3)
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
	double logK0_H2O = -2.209 + 3.097E-2 * tc - 1.098E-4 * pow(tc, 2) + 2.048E-7 * pow(tc, 3);
	K0_H2O = pow(10., logK0_H2O);  // equilibrium constant for H2O at 1 bar

	// Construct constant part of fugacity coefficients
	k_H = std::vector<double>(NC);
	labda = std::vector<double>(NC);
	ksi = std::vector<double>(NC);
	for (int i = 0; i < NC; i++) 
	{
		std::string comp = components[i];
		if (comp != "H2O")  // non-H2O components
		{
			double dB = aq1_par::tau[comp] + aq1_par::Gamma[comp] * p + aq1_par::beta[comp] * sqrt(1000 / T);
			k_H[i] = exp((1 - aq1_par::eta[comp]) * log(f0_H2O) + aq1_par::eta[comp] * log(aq1_par::R * T / M * rho0_H2O) + 2 * rho0_H2O * dB);
			labda[i] = aq1_par::labda[comp][0] + aq1_par::labda[comp][1] * T + aq1_par::labda[comp][2] / T + aq1_par::labda[comp][3] * p + aq1_par::labda[comp][4] / p + \
				aq1_par::labda[comp][5] * p / T + aq1_par::labda[comp][6] * T / pow(p, 2) + aq1_par::labda[comp][7] * p / (630 - T) + aq1_par::labda[comp][8] * T * log(p) \
				+ aq1_par::labda[comp][9] * p / pow(T, 2) + aq1_par::labda[comp][10] * pow(p, 2) * T + aq1_par::labda[comp][11] * p * T;
			ksi[i] = aq1_par::ksi[comp][0] + aq1_par::ksi[comp][1] * T + aq1_par::ksi[comp][2] * p / T + aq1_par::ksi[comp][3] * p / (630.0 - T);
		}
	}
	return;
}

std::vector<double> AQ1::fugacityCoefficient(std::vector<double> x) {
	// Construct fugacity coefficients
	std::vector<double> phi(NC);

	// Find effective salt molalities
	if (NI > 0)
	{	
		m_a = 0.;
		m_c = 0.;
		m_ac = 0.;
		for (int i = 0; i < NI; i++)
		{	
			if (CompProp::charge[ions[i]] > 0)
			{
				m_c += m_i[i] * CompProp::charge[ions[i]];
				m_ac += m_i[i];
			}
			else
			{
				m_a += m_i[i];
			}
		}
		m_ac *= m_a;
		
		// Calculate total ions mole fraction in Aq phase to calculate actual x_H2O in Aq phase to approximate H2O activity
		double tot_moles = 0.;  // total number of moles
		for (int i = 0; i < NC; i++)
		{
			tot_moles += 55.509 * x[i] / x[water_index];
		}
		for (int i = 0; i < NI; i++)
		{
			tot_moles += m_i[i];
		}
		xs = 0.;  // total ions mole fraction
		for (int i = 0; i < NI; i++)
		{
			xs += m_i[i] / tot_moles;
		}
	}

	// Calculate fugacity coefficients
	for (int i = 0; i < NC; i++) 
	{
		if (components[i] == "H2O") 
		{
			phi[i] = K0_H2O * (1-xs) * exp((p - 1.0) * V_H2O / (aq1_par::R * T)) / p;
		}
		else 
		{
			double gamma = exp((2 * m_c * labda[i]) + (m_ac * ksi[i]));
			phi[i] = k_H[i] * gamma / p;
		}
	}

	return phi;
}