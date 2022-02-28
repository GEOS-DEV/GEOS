#pragma once

namespace Composition {
	extern std::vector<std::string> components;
	extern int compsize;
	extern std::vector<std::string> phases;
	extern int phasesize;
}

namespace CompProp {
	// components{ 		    "H2O",   "CO2",   "N2",    "H2S",  "C1",   "C2",   "C3",   "nC4",  "nC5",  "nC6",  "nC7" };
	extern std::vector<double> Tc;  // critical temperature
	extern std::vector<double> Pc;   // critical pressure
	extern std::vector<double> ac;  // acentric factor
	extern std::vector<double> Mw; // molecular weight
	
	// Henry's constant	
	// https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015.pdf
	extern double T_0;
	extern std::vector<double> H0; // from Sander (2015)
	extern std::vector<double> dlnH; // from Sander (2015)
}