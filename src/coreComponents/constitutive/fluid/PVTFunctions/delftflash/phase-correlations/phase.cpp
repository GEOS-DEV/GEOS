#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "phase.h"
#include "eos.h"
#include "../global/global.h"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
using namespace std;

namespace H_par {
	std::unordered_map<std::string, std::vector<double>> aik = {
		{"H2O", {31.040960, -39.142208, 37.969528, -21.837491, 7.422515, -1.381789, 0.108807, -12.077118, 3.391051, -0.584521, 0.058993, -0.003130, 6.5746E-5}},
		{"CO2", {-1.8188731, 12.903022, -9.6634864, 4.2251879, -1.042164, 0.12683515, -0.49939675E-2, 2.4950242, -0.8272375, 0.15372481, -0.015861243, 0.8601715E-3, -0.19222165E-4}},
		{"N2", {-9.2984251, 20.007476, -16.763488, 8.6903787, -2.7510686, 0.48793873, -0.037167758, 4.0387289, -0.30781129, -0.19090602, 0.06465393, -0.0082736889, 0.00039772373}},
		{"H2S", {49.7193149691389, 49.0872303016932, 22.3724117815171, -3.24712589494E2, 0.525395845693E3, -0.36598306734E3, 0.990252722795E2, -0.93016407531E2, 0.624862496660E2, -0.23034052163E2, 4.98632802419611, -0.59734764063354, 0.030712261233025}},
		{"C1", {-133.833552230702, -44.3388396062382, 214.421008792969, 62.3762952145953, -821.181789165805, 1123.86680403028, -507.399195660801, 177.100692102039, -114.06659331988, 44.0472526925844, -10.2011315911762, 1.3056224566755, -0.07103768338394}},
	};
	std::unordered_map<std::string, double> H_int = {{"H2O", 9908.00}, {"CO2", 0.210820700252E4}, {"N2", 0.642444280732E4}, {"H2S", -2.9785008714E2}, {"C1", 28351.109188499}};
}

EoS* Phase::setEoS(std::vector<std::string> components, std::string eos_) {
	if (eos_ == "PR") 
	{
		return new PR(components);
	}
	else if (eos_ == "SRK") 
	{
		return new SRK(components);
	}
	else if (eos_ == "AQ1") 
	{
		return new AQ1(components, ions);
	}
	else if (eos_ == "AQ2") 
	{
		return new AQ2(components, ions);
	}
	else if (eos_ == "AQ3")
	{
		return new AQ3(components, ions);
	}
	else 
	{ // if (eos_ == "VdWP") {
		return new VdWP(components, phase);
	}
}

Properties::Properties(std::vector<std::string> comp, std::string ph, std::string eos_) {
	components = comp; phase = ph; NC = comp.size();
	Mw = std::vector<double>(NC);
	for (int i = 0; i < NC; i++) 
	{
		Mw[i] = CompProp::prop["Mw"][comp[i]];
		if (comp[i] == "H2O") { H2O = true; water_index = i; }
	}
	
	eos = setEoS(components, eos_);
};

// std::vector<double> ComponentProp::get(std::string prop) {
// 	std::vector<double> properties(NC);
// 	for (int i = 0; i < NC; i++) { properties[i] = CompProp::prop[prop][components[i]]; }
	
// 	return properties;
// }

double ComponentProp::mix_MW(std::vector<double> x) {
	double M = 0;
	for (int i = 0; i < NC; i++) 
	{ 
		M += x[i] * CompProp::prop["Mw"][components[i]]; 
	}
	return M;
}

double Properties::density(double p, double T, std::vector<double> x, double m) {
	double rho; // [kg/m3]

	// for V and L hydrocarbon phases, use PR/SRK EoS
	if (phase == "V" || phase == "L") 
	{
		eos->parameters(p, T);
		ComponentProp prop(components);
		double M = prop.mix_MW(x); // [g/mol]
		double z = eos->calc_z(x);
		double R = 8.3145E-5;
		double Vm = z * R * T / p; // [m3/mol]
		rho = M / 1000 / Vm; // [kg/mol] / [m3/mol] = [kg/m3]
	}
	else 
	{ // if (phase == "Aq") {
		// find index of CO2 first
		bool CO2 = false; int CO2_index;
		for (int i = 0; i < NC; i++) 
		{ 
			if (components[i] == "CO2") { CO2 = true; CO2_index = i; }
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

		if (CO2) 
		{
			// Ref. Garcia et al. (2001)
			double rho_b = rho;
			ComponentProp prop(components);
			double M = prop.mix_MW(x)/1000;  // molar weight kg/mol

        	// Apparent molar volume of dissolved CO2
        	double V_app = (37.51 - 9.585e-2 * tc + 8.740e-4 * pow(tc, 2.) - 5.044e-7 * pow(tc, 3.)) * 1e-6;  // in [m3/mol]
        	double M_b = 1000 / (1000 / 18.01 + m) + 58.44 * m / (1000 / 18.01 + m);  // molar weight kg/mol
			M_b /= 1000;

        	// density of aq phase
        	rho = 1 / ((x[CO2_index] * V_app) / M + (M_b * x[water_index]) / (rho_b * M));  // in [kg/m3]
		}
	}
	return rho;
}

double Properties::viscosity(double p, double T, std::vector<double> x, double rho, double m) {
	double mu; // cP

	// for V phase
	if (phase == "V") 
	{
		ComponentProp prop(components);
		double M = prop.mix_MW(x);
		double T_ran = T*1.8;  // rankine scale
		double a = (9.379 + 0.0160 * M) * pow(T_ran, 1.5) / (209.2 + 19.26 * M + T_ran);
        double b = 3.448 + 0.01009 * M + (986.4 / T_ran);
        double c = 2.447 - 0.2224 * b;
		rho *= 0.0624279606;  // convert rho from [kg/m3] -> [lb/ft^3]
		mu = 1E-4 * a * exp(b * pow((rho / 62.43), c));
	}
	else if (phase == "L") 
	{
		mu = 0.5;
	}
	else if (phase == "Aq") 
	{
		double tc = T-273.15;
    	double S = m * 58.40/1e3;
    	mu = 0.1 + 0.333*S + (1.65 + 91.9*pow(S, 3.)) * exp(-1*(0.42*pow((pow(S, 0.8) - 0.17), 2.) + 0.045) * pow(tc, 0.8));

		bool CO2 = false; int CO2_index;
		for (int i = 0; i < NC; i++) 
		{ 
			if (components[i] == "CO2") { CO2 = true; CO2_index = i; }
		}
		if (CO2) 
		{
			mu *= 1+4.65*pow(x[CO2_index], 1.0134);
		}
	}
	else 
	{
		mu = 1.0;
	}

	return mu;
}

double Properties::enthalpy(double p, double T, std::vector<double> x, double m) {
	double H;

	// for V and L hydrocarbon phases, use PR/SRK EoS
	if (phase == "V" || phase == "L") 
	{
		// Guo (2019): H_gas = H_ideal + H_deviation
		// H_ideal
		double T_  = 1000;
		double tau = T/T_;
		double R = 8.3145;

		std::vector<double> ai(13);
		double H_int;
		for (int k = 0; k < NC; k++) 
		{
			for (int i = 0; i < 13; i++) 
			{ 
				ai[i] += x[k]*H_par::aik[components[k]][i]; 
			}
			H_int += x[k] * H_par::H_int[components[k]];
		}

		double H_ideal = 0;
		for (int i = 0; i < 7; i++) 
		{ 
			H_ideal += ai[i]/(i+1) * pow(tau, i+1); 
		}
		H_ideal += ai[7] * log(tau);
		for (int i = 8; i < 13; i++) 
		{ 
			H_ideal += ai[i]/(7-i) * pow(1/tau, i-7); 
		}
		H_ideal = H_ideal*R*T_ + R*H_int;

		// H_dev (from EoS)
		double H_dev = eos->calc_H(p, T, x);

		// // H_vap_H2O
		// double H_vap_H2O = 0;
		// if (H2O) { H_vap_H2O = 40.8 * x[water_index]; }

		H = (H_ideal + H_dev)/1000;  // - H_vap_H2O;
	}
	else 
	{ // if (phase == "Aq") {
		// Guo (2019): H_aq = H_brine + H_dev_gas 
		int CO2_index;
		for (int i = 0; i < NC; i++) 
		{ 
			if (components[i] == "CO2") 
			{ 
				CO2_index = i; 
			}
		}

		double tc = T-273.15;

		// Pure water enthalpy and salt (Keenan, Keyes, Hill and Moore)
		double Hw = 0.12453e-4 * pow(tc, 3) - 0.4513e-2 * pow(tc, 2) + 4.81155 * tc - 29.578;  // kJ/kg
		double Hs = (-0.83624e-3 * pow(tc, 3) + 0.16792 * pow(tc, 2) - 25.9293 * tc) * (4.184 / 58.44);  // kJ/kg

		// Michaelidis/Lorenz: H_sat_brine = x1*Hw + x2*Hs + H_dev_salt
		double x1, x2, H_dev_salt{ 0 };
		if ((tc > 100.) && (tc < 300.))
		{
			x1 = 1000 / (1000 + 58.44 * m);  // mass fraction of water [-]
			x2 = 58.44 * m / (1000 + 58.44 * m);  // mass fraction of salt [-]
			if (m > 0.)
			{
				std::vector<double> aij = {-9633.6, -4080.0, 286.49, 166.58, 68.577, -4.6856, -0.90963, -0.36524, 0.0249667, 0.17965e-2, 0.71924e-3, -0.4900e-4};
				for (int ii = 0; ii < 4; ii++) 
				{
					for (int jj = 0; jj < 3; jj++) 
					{
						H_dev_salt += aij[ii*3 + jj] * pow(tc, ii) * pow(m, jj);
					}
				}
				H_dev_salt *= 4.184 / (1000 + 58.44 * m);
			}
		}
		else if ((tc > 0.) && (tc < 100.))
		{
			x1 = 1;
			x2 = 0;
			if (m > 0.)
			{
				std::vector<double> bij = {0.2985, -7.819E-2, 3.479E-4, -1.203E-6, -7.257E-2, 2.169E-3, -1.809E-5, 5.910E-8, 1.071E-3, -3.343E-5, 3.45E-7, -1.131E-9};
				for (int ii = 0; ii < 3; ii++) 
				{
					for (int jj = 0; jj < 4; jj++)
					{
						H_dev_salt += bij[ii*4 + jj] * pow(m, ii) * pow(tc, jj);
					}
				}
			}
		}
		else { cout << "T out of bounds for Aq enthalpy\n"; }

		double H_sat_brine = x1*Hw + x2*Hs + H_dev_salt;

		// temperature, pressure and salinity effects on brine enthalpy
		// Spivey, McCain and North (2004): H_brine ~= H_sat_brine + (V - T(dV/dT))*(p-psat)
		// dV/dT ~= (V2-V)/(T2-T)
		double T2 = T + 0.001;
		double rho_b = density(p, T, x, m);
		double rho_b2 = density(p, T2, x, m);
		double M_b = Mw[water_index];  // g/mol
		if (m > 0.)
		{	
			double x_s = m * x[water_index] / 55.509;
			double M_s = CompProp::prop["Mw"]["NaCl"];
			M_b = (x[water_index]*Mw[water_index] + x_s*M_s)/(x[water_index] + x_s);
		}
		double V = 1 / rho_b;  // m3/kg
		double V2 = 1 / rho_b2;
		double dVdT = (V2 - V) / (T2-T);

		// Antoine equation for psat
		double A, B, C;
		if (tc < 100.)
		{
			A = 8.07131; B = 1730.63; C = 233.426;
		}
		else if (tc < 370)
		{
			A = 8.140191; B = 1810.94; C = 244.485;
		}
		else { cout << "T out of bounds for Antoine\n"; }
		double psat = 0.00133322 * pow(10., A - B / (tc + C));  // bar

		double H_brine = H_sat_brine + (V - tc*dVdT) * (p - psat)*1e5/1e3;  // kJ/kg
		H_brine *= M_b/1000;  // kJ/kg -> kJ/mol

		// Enthalpy change with gas dissolution in Aq phase
		// Guo (2019): H_dev_gas = d(lnH)/d(1/T) * R or Battistelli (1997): H_dev_gas = (log(Hcoeff2) - log(Hcoeff1)) / (T2 - T) * R * pow(T, 2)
		eos->parameters(p, T);
		eos->setMolality(m_i);
		std::vector<double> phi = eos->fugacityCoefficient(x);

		eos->parameters(p, T2);
		std::vector<double> phi2 = eos->fugacityCoefficient(x);

		double Hcoeff1 = phi[CO2_index]/p;
		double Hcoeff2 = phi2[CO2_index]/p;

		double R = 8.3145;
		// double H_diss_g = (log(Hcoeff2) - log(Hcoeff1)) / (1/T2 - 1/T)*R;  // J/mol
		double H_diss_g = (log(Hcoeff2) - log(Hcoeff1)) / (T2 - T) * R * pow(T, 2);  // J/mol
		H_diss_g /= 1000;  // kJ/mol

		H = x[water_index] * H_brine + x[CO2_index] * H_diss_g;
	}
	return H;
}

double Properties::thermal_conductivity(double p, double T, std::vector<double> x, double rho, double m) {
	double kappa; // kJ/m/day/K

	if (phase == "V") 
	{
		std::vector<double> A = {105.161, 0.9007, 0.0007, 3.5e-15, 3.76e-10, 0.75, 0.0017};
		kappa = (A[0] + A[1] * rho + A[2] * pow(rho, 2) + A[3] * pow(rho, 3) * pow(T, 3) + A[4] * pow(rho, 4) + A[5] * T + A[6] * pow(T, 2)) / sqrt(T);
		kappa *= 1e-3 / 1000 * 3600 * 24; // convert from W/m/K to kJ/m/day/K
	}
	else if (phase == "L") 
	{
		kappa = 1;
	}
	else  //if (phase == "Aq") {}
	{
		double S = m * 58.44;
		double T_d = T / 300;
        double cond_aq = 0.797015 * pow(T_d, -0.194) - 0.251242 * pow(T_d, -4.717) + 0.096437 * pow(T_d, -6.385) - 0.032696 * pow(T_d, -2.134);
        kappa = (cond_aq / (0.00022 * S + 1)) +0.00005*(p-50);
	}
	return kappa;
}

void KineticHydrate::separateIons(std::vector<std::string> species) {
	// check if ions present in species vector
	int ns = species.size();
	for (int i = 0; i < ns; i++)
	{
		if (CompProp::species_type[species[i]] == 'i')
		{	
			ions.push_back(species[i]);
		}
		else
		{
			components.push_back(species[i]);
		}
	}
	NC = components.size();
	NI = ions.size();
	return;
}

double KineticHydrate::fwH(double p, double T, std::vector<double> f0) {
	h_eos->parameters(p, T);

	// If ions are present, remove them from f0 vector
	if (NI > 0)
	{
		std::vector<double> f0_(NC);
		for (int i = 0; i < NC; i++)
		{
			f0_[i] = f0[i];
		}
		f0 = f0_;
	}

	double f_wH = h_eos->fwH(f0);

	return f_wH;
}
