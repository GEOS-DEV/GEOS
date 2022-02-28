#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "../flash/ssi.h"
#include "phase.h"
#include "eos.h"
#include "../global/global.h"
#include "../global/misc.h"

using namespace std;

EoS* Phase::setEoS(std::vector<int> ci, int eos1) {
	if (phase_index == 0 || phase_index == 1) {
		if (eos1 == 0) {
			return new PR(ci, phase_index);
		}
		else {
			return new SRK(ci, phase_index);
		}
	} else if (phase_index == 2) {
		if (eos1 == 0) {
			return new AQ1(ci);
		} 
		else {
			return new AQ2(ci);
		}
	}
	else { //if (phase_index == 3 || phase_index == 4 || phase_index == 5) {
		return new Hydrate(ci, phase_index);
	}
}

PhaseProperties::PhaseProperties(std::vector<std::string> comp, std::string ph, int phase_eos) {
	NC = comp.size();
	ci = std::vector<int>(NC); Mw = std::vector<double>(NC);
	for (int i = 0; i < NC; i++) {
		ci[i] = searchString(Composition::components, Composition::compsize, comp[i]); // finds index of each component for correlations
		Mw[i] = CompProp::Mw[ci[i]];
	}
	phaseIndex = searchString(Composition::phases, Composition::phasesize, ph);
	
	eos = setEoS(ci, phaseIndex, phase_eos);
}

EoS* PhaseProperties::setEoS(std::vector<int> cii, int phase_index, int phase_eos) {
	if (phase_index == 0 || phase_index == 1) {
		if (phase_eos == 0) {
			return new PR(cii, phase_index);
		}
		else {
			return new SRK(cii, phase_index);
		}
	} else if (phase_index == 2) {
		if (phase_eos == 0) {
			return new AQ1(cii);
		} 
		else {
			return new AQ2(cii);
		}
	}
	else { //if (phase_index == 3 || phase_index == 4 || phase_index == 5) {
		return new Hydrate(cii, phase_index);
	}
}

double PhaseProperties::molar_weight(std::vector<double> x) {
	double M = 0;
	for (int i = 0; i < NC; i++) { M += x[i] * Mw[i]; }
	return M;
}

double PhaseProperties::density(double p, double T, std::vector<double> x) {
	double rho; // [kg/m3]

	// for V and L hydrocarbon phases, use PR/SRK EoS
	if (phaseIndex == 0 || phaseIndex == 1) {
		eos->init(p, T);
		double M = molar_weight(x); // [g/mol]
		double z = eos->Z(p, T, x);
		double R = 8.3145E-5;
		double Vm = z * R * T / p; // [m3/mol]
		rho = M / 1000 / Vm; // [kg/mol] / [m3/mol] = [kg/m3]
	}
	else { // if (phaseIndex == 2) {
		// find index of water and CO2 first
		bool CO2 = false; int CO2_index;
		for (int i = 0; i < NC; i++) { 
			if (ci[i] == 0) { water_index = i; }
			else if (ci[i] == 1) { CO2 = true; CO2_index = i; }
		}

		// density of pure water
		// a = (a1 * (Tc / 100) ** 2 + a2 * (Tc / 100) + a3) / (a4 * (Tc / 100) ** 2 + a5 * (Tc / 100) + 1)
		double tc = T - 273.15;  // Temp in [Celcius]
		double p0 = 700;  // reference pressure of 70 MPa

		double rho_w0 = (-0.127213 * pow((tc / 100), 2.) + 0.645486 * (tc / 100) + 1.03265) / (-0.070291 * pow((tc / 100), 2.) + 0.639589 * (tc / 100) + 1);  // density of pure water at 70Mpa
    	double Ew_w = (4.221 * pow((tc / 100), 2.) - 3.478 * (tc / 100) + 6.221) / (0.5182 * pow((tc / 100), 2.) - 0.4405 * (tc / 100) + 1);
    	double Fw_w = (-11.403 * pow((tc / 100), 2.) + 29.932 * (tc / 100) + 27.952) / (0.20684 * pow((tc / 100), 2.) + 0.3768 * (tc / 100) + 1);
    	double Iw = (1 / Ew_w) * log(abs(Ew_w * (p / p0) + Fw_w));
    	double Iw0 = (1 / Ew_w) * log(abs(Ew_w * (p0 / p0) + Fw_w));

    	rho = 1000 * rho_w0 * exp(Iw - Iw0); // pure water, kg/m3

		// if (std::any_of(ci.begin(), ci.end(), [](int i){return i==1;})) {
		if (CO2) {  // CO2
			// Ref. Garcia et al. (2001)
			double rho_b = rho;
			double M = molar_weight(x)/1000;  // molar weight kg/mol
        	// Apparent molar volume of dissolved CO2
        	double V_app = (37.51 - 9.585e-2 * tc + 8.740e-4 * pow(tc, 2.) - 5.044e-7 * pow(tc, 3.)) * 1e-6;  // in [m3/mol]
        	// double M_b = 1000 / (1000 / 18.01 + Cm) + 58.44 * Cm / (1000 / 18.01 + Cm);
			double M_b = Mw[water_index]/1000;  // molar weight kg/mol

        	// density of aq phase
        	rho = 1 / ((x[CO2_index] * V_app) / M + (M_b * x[water_index]) / (rho_b * M));  // in [kg/m3]
		}
	}
	return rho;
}

double PhaseProperties::viscosity(double p, double T, std::vector<double> x) {
	double mu; // cP

	// for V phase
	if (phaseIndex == 0) {
		double M = molar_weight(x);
		double rho = density(p, T, x);
		double T_ran = T*1.8;  // rankine scale
		double a = (9.379 + 0.0160 * M) * pow(T_ran, 1.5) / (209.2 + 19.26 * M + T_ran);
        double b = 3.448 + 0.01009 * M + (986.4 / T_ran);
        double c = 2.447 - 0.2224 * b;
		rho *= 0.0624279606;  // convert rho from [kg/m3] -> [lb/ft^3]
		mu = 1E-4 * a * exp(b * pow((rho / 62.43), c));
	}
	else if (phaseIndex == 1) {
		mu = 0.5;
	}
	else if (phaseIndex == 2) {
		double tc = T-273.15;
    	// S = Cm * 58.40/1e3;
		double S = 0;
    	mu = 0.1 + 0.333*S + (1.65 + 91.9*pow(S, 3.)) * exp(-1*(0.42*pow((pow(S, 0.8) - 0.17), 2.) + 0.045) * pow(tc, 0.8));

		bool CO2 = false; int CO2_index;
		for (int i = 0; i < NC; i++) { 
			if (ci[i] == 1) { CO2 = true; CO2_index = i; }
		}
		if (CO2) {
			mu *= 1+4.65*pow(x[CO2_index], 1.0134);
		}
	}
	else {
		mu = 1.0;
	}

	return mu;
}

// double PhaseProperties::enthalpy(double p, double T, std::vector<double> x) {
// 	double H;
	
// }

// double PhaseProperties::heat_capacity(double p, double T, std::vector<double> x) {
// 
// }

// double PhaseProperties::thermal_conductivity(double p, double T, std::vector<double> x) {
// 
// }

// std::vector<double> PhaseProperties::fugacity(double p, double T, std::vector<double> x) {
// 	std::vector<double> f(NC);

// 	eos->init(p, T);
// 	std::vector<double> phi = eos->fugacityCoefficient(p, T, x);

// 	for (int i = 0; i < NC; i++) { f[i] = phi[i] * x[i] * p; }

// 	return f;
// }