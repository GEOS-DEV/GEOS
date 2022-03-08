#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "eos.h"
#include "../global/global.h"

using namespace std;

namespace sI {
	// sI Hydrate parameters
	// Ballard (2002) - A.1.5
	std::unordered_map<std::string, std::vector<double>> a = {
		{"CO2", {15.8336435, 3.119, 0, 3760.6324, 1090.27777, 0, 0}},	
		{"H2S", {31.209396, -4.20751374, 0.761087, 8340.62535, -751.895, 182.905, 0}},
		{"C1", {27.474169, -0.8587468, 0, 6604.6088, 50.8806, 1.57577, -1.4011858}},
		{"C2", {14.81962, 6.813994, 0, 3463.9937, 2215.3, 0, 0}}
	};
	std::unordered_map<bool, std::vector<double>> a_N2 = {
		{false, {173.2164, -0.5996, 0, 24751.6667, 0, 0, 0, 1.441, -37.0696, -0.287334, -2.07405E-5, 0, 0}},  // without H2S
		{true, {71.67484, -1.75377, -0.32788, 25180.56, 0, 0, 0, 56.219655, -140.5394, 0, 8.0641E-4, 366006.5, 978852}}  // with H2S
	};
}

namespace sII {
	// sII Hydrate parameters
	// Ballard (2002) - A.1.5
	std::unordered_map<std::string, std::vector<double>> a = {
		{"CO2", {9.0242, 0, 0, -207.033, 0, 6.7588E-4, -6.992E-3, -6.0794E-4, -9.026E-2, 0, 0, 0, 0.0186833, 0, 0, 8.82E-5, 7.78015E-3, 0, 0}},
		{"N2", {1.78857, 0, -0.019667, -6.187, 0, 0, 0, 5.259E-5, 0, 0, 0, 0, 0, 0, 192.39, 0, 3.051E-5, 1.1E-7, 0}},
		{"H2S", {-6.42956, 0.06192, 0, 82.627, 0, -1.0718E-4, 0, 0, 3.493522, -0.64405, 0, 0, 0, -184.257, 0, -1.30E-6, 0, 0, 0}},
		{"C1", {-0.45872, 0, 0, 31.6621, -3.4028, -7.702E-5, 0, 0, 1.8641, -0.78338, 0, 0, 0, -77.6955, 0, -2.3E-7, -6.102E-5, 0, 0}},
		{"C2", {3.21799, 0, 0, -290.283, 181.2694, 0, 0, -1.893E-5, 1.882, -1.19703, -402.166, -4.897688, 0.0411205, -68.8018, 25.6306, 0, 0, 0, 0}},
		{"C3", {-7.51966, 0, 0, 47.056, 0, -1.697E-5, 7.145E-4, 0, 0, 0.12348, 79.34, 0, 0.0160778, 0, -14.684, 5.50E-6, 0, 0, 0}},
		{"nC4", {-37.211, 0.86564, 0, 732.2, 0, 0, 0, 1.9711E-3, -15.6144, 0, 0, -4.56576, 0, 0, 300.55350, 0, 0.0151942, -1.26E-6, 0}},
		{"iC4", {-9.55128, 0, 0, 0, 0, 0, 0.001251, 2.1036E-6, 2.40904, -2.75945, 0, 0, 0, 0, -0.28974, 0, -1.6476E-3, -1.0E-8, 0}},
	};
}

InitialK::InitialK(double p_, double T_, std::vector<double> zc, std::vector<std::string> comp) {
	p = p_; T = T_; z = zc;
	components = comp; NC = components.size();
	water_index = std::distance(components.begin(), std::find(components.begin(), components.end(), "H2O"));
}

std::vector<double> InitialK::K_initial(std::vector<std::string> phases, int j) {
	// Calculates initial K-values for phase j
	// Ballard (2002) only provides K-values based on reference vapour phase
	// This function derives K-values from any reference phase
	std::vector<double> K_j(NC);

	// Reference phase
	std::vector<double> k_0(NC);
	if (phases[0] == "V") 
	{
		for (int i = 0; i < NC; i++) 
		{ 
			k_0[i] = 1.; 
		}
	} 
	else if (phases[0] == "L") 
	{
		k_0 = vapour_liquid();
	} 
	else if (phases[0] == "Aq") 
	{
		k_0 = vapour_aqueous();
	} 
	else if (phases[0] == "sI") 
	{
		k_0 = vapour_sI();
	} 
	else if (phases[0] == "sII") 
	{
		k_0 = vapour_sII();
	}

	// J-th phase
	std::vector<double> k_j(NC);
	if (phases[j] == "V") 
	{
		for (int i = 0; i < NC; i++) 
		{ 
			k_j[i] = 1.; 
		}
	} 
	else if (phases[j] == "L") 
	{
		k_j = vapour_liquid();
	} 
	else if (phases[j] == "Aq") 
	{
		k_j = vapour_aqueous();
	} 
	else if (phases[j] == "sI") 
	{
		k_j = vapour_sI();
	} 
	else if (phases[j] == "sII") 
	{
		k_j = vapour_sII();
	}

	// Divide K_j = x_j/x_0 = k_0/k_j
	for (int i = 0; i < NC; i++) 
	{ 
		K_j[i] = k_0[i] / k_j[i]; 
	}
	
	return K_j;
}

std::vector<double> InitialK::vapour_liquid() {
	std::vector<double> K_(NC);
	// Wilson's equation for phase L
	// Ballard (2002) - Appendix A.1.1
	for (int ii = 0; ii < NC; ii++) 
	{
		std::string comp = components[ii];
		if (comp == "H2O") 
		{
			K_[ii] = (-133.67 + 0.63288 * T) / p + 3.19211E-3 * p;
		}
		else 
		{
			K_[ii] = CompProp::prop["Pc"][comp] / p * exp(5.373 * (1 + CompProp::prop["ac"][comp]) * (1 - CompProp::prop["Tc"][comp] / T));
		}
	}
	return K_;
}

std::vector<double> InitialK::vapour_aqueous() {
	// Initial K-value for phase Aq
	// Ballard (2002) - Appendix A.1.2
	std::vector<double> K_(NC);
	for (int ii = 0; ii < NC; ii++) 
	{
		std::string comp = components[ii];
		if (comp == "H2O") 
		{ 
			// Raoult's law for H2O
			double psat = exp(12.048399 - 4030.18245 / (T - 38.15));
			double j_inf = 1;
			K_[ii] = psat / p * j_inf;
		}
		else 
		{ 
			// Henry's law for solutes in dilute solution
			double x_iV = 1.;
			double H = CompProp::H0[comp] * exp(CompProp::dlnH[comp] * (1/T - 1/CompProp::T_0));
			double ca = H*p;
			double rho_Aq = 1100;
			double Vm = CompProp::prop["Mw"]["H2O"]*1E-3/rho_Aq;
			double x_iAq = ca*Vm;			
			K_[ii] = x_iV/x_iAq;
		}
	}

	return K_;
}

std::vector<double> InitialK::vapour_sI() {
	// Initial K-value for sI
	// Ballard (2002) - Appendix A.1.5
	std::vector<double> K_(NC);
	double x_wH = 0.88;
	
	for (int ii = 0; ii < NC; ii++) 
	{
		std::string comp = components[ii];
		if (comp == "H2O") 
		{
			// Kw_VAq
			double psat = exp(12.048399 - 4030.18245 / (T - 38.15));
			double j_inf = 1;
			double Kw_VAq = psat / p * j_inf;
			// Kw_IAq
			double p0 = 6.11657E-3;
			double T0 = 273.1576;
			double Ti = T0 - 7.404E-3*(p-p0) - 1.461E-6*pow(p-p0, 2.);
			double x_wAq = 1 + 8.33076E-3*(T-Ti) + 3.91416E-5*pow(T-Ti, 2.);
			double Kw_IAq = 1/x_wAq;
			K_[ii] = Kw_VAq/(x_wH*Kw_IAq);
		}
		else if (comp == "N2") 
		{
			std::vector<double> a;

			bool H2S = std::find(components.begin(), components.end(), "H2S") != components.end();
			a = sI::a_N2[H2S]; // depends on presence of H2S

			double Ki_wf = exp(a[0] + a[1]*log(p) + a[2]*pow(log(p), 2.) - (a[3] + a[4]*log(p) + a[5]*pow(log(p), 2.) + a[6]*pow(log(p), 3.))/T 
				+ a[7]/p + a[8]/(pow(p, 2)) + a[9]*T + a[10]*p + a[11]*log(p)/pow(T, 2) + a[12]/pow(T, 2));
			K_[ii] = Ki_wf/(1-x_wH);
		}
		else 
		{
			double Ki_wf = exp(sI::a[comp][0] + sI::a[comp][1]*log(p) + sI::a[comp][2]*pow(log(p), 2.) - (sI::a[comp][3] + sI::a[comp][4]*log(p) 
				+ sI::a[comp][5]*pow(log(p), 2.) + sI::a[comp][6]*pow(log(p), 3.))/T);
			K_[ii] = Ki_wf/(1-x_wH);
		}
	}
	return K_;
}

std::vector<double> InitialK::vapour_sII() {
	// Initial K-value for sII
	// Ballard (2002) - Appendix A.1.5
	std::vector<double> K_(NC);
	double x_wH = 0.90;
	
	for (int ii = 0; ii < NC; ii++) 
	{
		std::string comp = components[ii];
		if (comp == "H2O") 
		{
			// Kw_VAq
			double psat = exp(12.048399 - 4030.18245 / (T - 38.15));
			double j_inf = 1;
			double Kw_VAq = psat / p * j_inf;
			// Kw_IAq
			double p0 = 6.11657E-3;
			double T0 = 273.1576;
			double Ti = T0 - 7.404E-3*(p-p0) - 1.461E-6*pow(p-p0, 2.);
			double x_wAq = 1 + 8.33076E-3*(T-Ti) + 3.91416E-5*pow(T-Ti, 2.);
			double Kw_IAq = 1/x_wAq;
			K_[ii] = Kw_VAq/(x_wH*Kw_IAq);
		}
		else 
		{
			double Ki_wf = exp(sII::a[comp][0] + sII::a[comp][1]*T + sII::a[comp][2]*p + sII::a[comp][3]/T + sII::a[comp][4]/p + sII::a[comp][5]*p*T + sII::a[comp][6]*pow(T, 2) 
				+ sII::a[comp][7] * pow(p, 2) + sII::a[comp][8]*p/T + sII::a[comp][9]*log(p/T) + sII::a[comp][10]/pow(p, 2) + sII::a[comp][11]*T/p + sII::a[comp][12]*pow(T, 2)/p 
				+ sII::a[comp][13]*p/pow(T, 2) + sII::a[comp][14]*T/pow(p, 3) + sII::a[comp][15]*pow(T, 3) + sII::a[comp][16]*pow(p, 3)/pow(T, 2) + sII::a[comp][17]*pow(T, 4)) + sII::a[comp][18]*log(p);
			K_[ii] = Ki_wf/(1-x_wH);
		}
	}
	return K_;
}