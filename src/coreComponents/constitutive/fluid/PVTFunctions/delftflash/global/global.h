#include <unordered_map>
#pragma once

#define M_PI 3.14159265358979323846 // pi
#define M_NA 6.022140857e23 // Avogadro's number [mol-1]
#define M_kB 1.38064852e-23 // Boltzmann constant [J/K]

namespace CompProp {
	extern std::unordered_map<std::string, char> species_type;
	extern std::unordered_map<std::string, int> charge;
	extern std::unordered_map<std::string, std::unordered_map<std::string, double>> prop;
	extern std::unordered_map<std::string, std::unordered_map<std::string, double>> bic;
	
	// Henry's constant	
	// https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015.pdf
	extern double T_0;
	extern std::unordered_map<std::string, double> H0; // from Sander (2015)
	extern std::unordered_map<std::string, double> dlnH; // from Sander (2015)
}

namespace aq1_par {
	extern std::unordered_map<std::string, std::vector<double>> labda, ksi;
	extern std::unordered_map<std::string, double> eta, tau, beta, Gamma;
	extern double R;
}

namespace aq2_par {
	extern double R, T_0, P_0;
    extern std::unordered_map<std::string, double> gi_0, hi_0, gi_00, hi_00;
    extern std::unordered_map<std::string, std::vector<double>> hi_a;
    
    extern std::unordered_map<std::string, double> omega, cp1, cp2;
    extern std::unordered_map<std::string, std::vector<double>> vp;
    extern std::vector<double> eps;

    extern std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Bca, Cca, Dca, B;

    extern double e, eps0;
}