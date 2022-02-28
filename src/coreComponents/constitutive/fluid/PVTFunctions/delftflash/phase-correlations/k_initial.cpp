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
	// components{ 		   "H2O", "CO2",      "N2",       "H2S",      "C1",      "C2",      "C3", "nC4", "nC5", "nC6", "nC7" };
	std::vector<double> a1{ 0,    15.8336435, 173.2164,   31.209396,   27.474169, 14.81962,  0,    0,     0,     0,     0 };
	std::vector<double> a2{ 0,	  3.119,	 -0.5996,    -4.20751374, -0.8587468, 6.813994,  0,    0,     0,     0,     0 };
	std::vector<double> a3{	0, 	  0,		  0,		  0.761087,    0,		  0, 	 	 0,    0,     0,     0,     0 };
	std::vector<double> a4{	0, 	  3760.6324,  24751.6667, 8340.62535,  6604.6088, 3463.9937, 0,    0,     0,     0,     0 };
	std::vector<double> a5{	0,    1090.27777, 0,	     -751.895,     50.8806,   2215.3,    0,    0,     0,     0,     0 };
	std::vector<double> a6{	0,    0, 		  0,	      182.905,     1.57577,   0,    	 0,    0,     0,     0,     0 };
	std::vector<double> a7{	0,    0,		  0,	      0,  		  -1.4011858, 0,    	 0,    0,     0,     0,     0 };
	// change N2 if H2S is present
	std::vector<double> N2_wo{ 173.2164, -0.5996, 0, 24751.6667, 0, 0, 0, 1.441, -37.0696, -0.287334, -2.07405E-5, 0, 0 };
	std::vector<double> N2_w{ 71.67484, -1.75377, -0.32788, 25180.56, 0, 0, 0, 56.219655, -140.5394, 0, 8.0641E-4, 366006.5, 978852 };
}

namespace sII {
	// sII Hydrate parameters
	// Ballard (2002) - A.1.5
	// components{ 		   "H2O", "CO2",      "N2",       "H2S",      "C1",      "C2",      "C3",     "nC4", 	 "nC5", 	"nC6", 		"nC7" };
	std::vector<double> a1{ 0,    9.0242, 	  1.78857,  -6.42956,   -0.45872, 	3.21799,   -7.51966,  -37.211,   -37.211,   -37.211,   -37.211 };
	std::vector<double> a2{ 0,    0, 		  0,   		 0.06192,    0, 		0,  		0,    	   0.86564,   0.86564,   0.86564, 	0.86564 };
	std::vector<double> a3{ 0,    0, 		 -0.019667,  0,   		 0, 		0,  		0,         0,     	  0,       	 0,       	0 };
	std::vector<double> a4{ 0,   -207.033,   -6.187,   	 82.627,     31.6621,  -290.283,  	47.056,    732.2,     732.2,     732.2,   	732.2 };
	std::vector<double> a5{ 0,    0, 		  0,   	 	 0,   	    -3.4028, 	181.2694,  	0,    	   0,         0,       	 0,       	0 };
	std::vector<double> a6{ 0,    6.7588E-4,  0,   		-1.0718E-4, -7.702E-5, 	0,  	   -1.697E-5,  0,     	  0,     	 0,     	0 };
	std::vector<double> a7{ 0,   -6.992E-3,   0,   		 0,   		 0, 		0,  		7.145E-4,  0,    	  0,     	 0,     	0 };
	std::vector<double> a8{ 0,   -6.0794E-4,  5.259E-5,  0,   		 0, 	   -1.893E-5,  	0,    	   1.9711E-3, 1.9711E-3, 1.9711E-3, 1.9711E-3 };
	std::vector<double> a9{ 0,   -9.026E-2,   0,   		 3.493522,   1.8641, 	1.882,  	0,    	  -15.6144,  -15.6144,  -15.6144,  -15.6144 };
	std::vector<double> a10{ 0,   0, 		  0,   		-0.64405,   -0.78338,  -1.19703,  	0.12348,   0,     	  0,     	 0,     	0 };
	std::vector<double> a11{ 0,   0, 		  0,   		 0,   		 0, 	   -402.166,  	79.34,     0,     	  0,     	 0,     	0 };
	std::vector<double> a12{ 0,   0, 		  0,   		 0,   		 0, 	   -4.897688,  	0,    	  -4.56576,  -4.56576,  -4.56576,  -4.56576 };
	std::vector<double> a13{ 0,   0.0186833,  0,   		 0,   	 	 0, 		0.0411205,  0.0160778, 0,     	  0,     	 0,     	0 };
	std::vector<double> a14{ 0,   0, 		  0,   		-184.257,   -77.6955,  -68.8018,  	0,    	   0,     	  0,     	 0,     	0 };
	std::vector<double> a15{ 0,   0, 		  192.39,    0,   		 0, 		25.6306,   -14.684,    300.55350, 300.55350, 300.55350, 300.55350 };
	std::vector<double> a16{ 0,   8.82E-5,    0,  		-1.30E-6,   -2.3E-7, 	0,  		5.50E-6,   0,     	  0,     	 0,     	0 };
	std::vector<double> a17{ 0,   7.78015E-3, 3.051E-5,  0,   		-6.102E-5, 	0,  		0,    	   0.0151942, 0.0151942, 0.0151942, 0.0151942 };
	std::vector<double> a18{ 0,   0, 		  1.1E-7,    0,   		 0, 		0,  		0,    	  -1.26E-6,  -1.26E-6,  -1.26E-6,  -1.26E-6 };
	// std::vector<double> a19{ 0,    0, 0,   0,   0, 0,  0,    0,     0,     0,     0 };
}

InitialK::InitialK(double p_, double T_, std::vector<double> z_, std::vector<int> ci_) {
	p = p_; T = T_; z = z_;
	ci = ci_; NC = ci.size();
	water_index = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 0));
}

std::vector<double> InitialK::K_initial(std::vector<int> pi, int j) {
	// Calculates initial K-values for phase j
	// Ballard (2002) only provides K-values based on reference vapour phase
	// This function derives K-values from any reference phase
	std::vector<double> K_j(NC);

	// Reference phase
	std::vector<double> k_0(NC);
	if (pi[0] == 0) {
		for (int i = 0; i < NC; i++) { k_0[i] = 1.; }
	} else if (pi[0] == 1) {
		k_0 = vapour_liquid();
	} else if (pi[0] == 2) {
		k_0 = vapour_aqueous();
	} else if (pi[0] == 3) {
		k_0 = vapour_sI();
	} else if (pi[0] == 4) {
		k_0 = vapour_sII();
	}

	// J-th phase
	std::vector<double> k_j(NC);
	if (pi[j] == 0) {
		for (int i = 0; i < NC; i++) { k_j[i] = 1.; }
	} else if (pi[j] == 1) {
		k_j = vapour_liquid();
	} else if (pi[j] == 2) {
		k_j = vapour_aqueous();
	} else if (pi[j] == 3) {
		k_j = vapour_sI();
	} else if (pi[j] == 4) {
		k_j = vapour_sII();
	}

	// Divide K_j = x_j/x_0 = k_0/k_j
	for (int i = 0; i < NC; i++) { K_j[i] = k_0[i] / k_j[i]; }
	
	return K_j;
}

std::vector<double> InitialK::vapour_liquid() {
	std::vector<double> K_(NC);
	// Wilson's equation for phase L
	// Ballard (2002) - Appendix A.1.1
	for (int ii = 0; ii < NC; ii++) {
		if (ci[ii] == 0) { // H2O
			K_[ii] = (-133.67 + 0.63288 * T) / p + 3.19211E-3 * p;
		}
		else {
			K_[ii] = CompProp::Pc[ci[ii]] / p * exp(5.373 * (1 + CompProp::ac[ci[ii]]) * (1 - CompProp::Tc[ci[ii]] / T));
		}
	}
	return K_;
}

std::vector<double> InitialK::vapour_aqueous() {
	// Initial K-value for phase Aq
	// Ballard (2002) - Appendix A.1.2
	std::vector<double> K_(NC);
	for (int ii = 0; ii < NC; ii++) {
		if (ci[ii] == 0) { // Raoult's law for H2O
			double psat = exp(12.048399 - 4030.18245 / (T - 38.15));
			double j_inf = 1;
			K_[ii] = psat / p * j_inf;
		}
		else { // Henry's law for solutes in dilute solution
			// modify for salt!!
			double x_iV = z[ii]/(1-z[water_index]);
			double H = CompProp::H0[ci[ii]] * exp(CompProp::dlnH[ci[ii]] * (1/T - 1/CompProp::T_0));
			double ca = H*p;
			double rho_Aq = 1100;
			double Vm = CompProp::Mw[0]*1E-3/rho_Aq;
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
	
	for (int ii = 0; ii < NC; ii++) {
		if (ci[ii] == 0) { // H2O
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
		else if (ci[ii] == 2) { // N2
			std::vector<double> a;
			int H2S = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 0));
			if (H2S) { a = sI::N2_w; } // depends on presence of H2S
			else { a = sI::N2_wo; }
			double Ki_wf = exp(a[0] + a[1]*log(p) + a[2]*pow(log(p), 2.) - (a[3] + a[4]*log(p) + a[5]*pow(log(p), 2.) + a[6]*pow(log(p), 3.))/T 
				+ a[7]/p + a[8]/(pow(p, 2)) + a[9]*T + a[10]*p + a[11]*log(p)/pow(T, 2) + a[12]/pow(T, 2));
			K_[ii] = Ki_wf/(1-x_wH);
		}
		else {
			double Ki_wf = exp(sI::a1[ci[ii]] + sI::a2[ci[ii]]*log(p) + sI::a3[ci[ii]]*pow(log(p), 2.) - (sI::a4[ci[ii]] + sI::a5[ci[ii]]*log(p) 
				+ sI::a6[ci[ii]]*pow(log(p), 2.) +  + sI::a7[ci[ii]]*pow(log(p), 3.))/T);
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
	
	for (int ii = 0; ii < NC; ii++) {
		if (ci[ii] == 0) {
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
		else {
			double Ki_wf = exp(sII::a1[ci[ii]] + sII::a2[ci[ii]]*T + sII::a3[ci[ii]]*p + sII::a4[ci[ii]]/T + sII::a5[ci[ii]]/p + sII::a6[ci[ii]]*p*T + sII::a7[ci[ii]]*pow(T, 2) 
				+ sII::a8[ci[ii]] * pow(p, 2) + sII::a9[ci[ii]]*p/T + sII::a10[ci[ii]]*log(p/T) + sII::a11[ci[ii]]/pow(p, 2) + sII::a12[ci[ii]]*T/p + sII::a13[ci[ii]]*pow(T, 2)/p 
				+ sII::a14[ci[ii]]*p/pow(T, 2) + sII::a15[ci[ii]]*T/pow(p, 3) + sII::a16[ci[ii]]*pow(T, 3) + sII::a17[ci[ii]]*pow(p, 3)/pow(T, 2) + sII::a18[ci[ii]]*pow(T, 4));  //+ sII::a19[ci[i]]*log(p);
			K_[ii] = Ki_wf/(1-x_wH);
		}
	}
	return K_;
}