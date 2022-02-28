#include <iostream>
#include <vector>
#include <complex>

#include "global.h"

namespace Composition {
	std::vector<std::string> components{ "H2O", "CO2", "N2", "H2S", "C1", "C2", "C3", "nC4", "nC5", "nC6", "nC7" };
	int compsize = components.size();
	std::vector<std::string> phases{ "V", "L", "Aq", "sI", "sII", "sH", "I", "S" };
	int phasesize = phases.size();
}

namespace CompProp {
	// components{ 		    "H2O",   "CO2",   "N2",    "H2S",  "C1",   "C2",   "C3",   "nC4",  "nC5",  "nC6",  "nC7" };
	std::vector<double> Tc{ 647.14,  304.10,  126.20,  373.53, 190.58, 305.32, 369.83, 425.12, 469.70, 507.60, 540.20 };  // critical temperature
	std::vector<double> Pc{ 220.50,  73.75,   34.00,   89.63,  46.04,  48.721, 42.481, 37.960, 33.701, 30.251, 27.40 };   // critical pressure
	std::vector<double> ac{ 0.328,   0.239,   0.0377,  0.0942, 0.012,  0.0995, 0.1523, 0.2002, 0.2515, 0.3013, 0.3495 };  // acentric factor
	std::vector<double> Mw{ 18.015,  44.01,   28.013,  34.10,  16.043, 30.07,  44.097, 58.124, 72.151, 86.178, 100.205 }; // molecular weight
	
	// Henry's constant	
	// https://acp.copernicus.org/articles/15/4399/2015/acp-15-4399-2015.pdf
	double T_0 = 298.15;
	std::vector<double> H0{ 0, 		 33.0, 	  0.64,    100.0,  1.4,	   1.9,    1.5,    1.2,    0.8,    0.61,   0.44 }; // from Sander (2015)
	std::vector<double> dlnH{ 0, 	 2400, 	  1600,    2100,   1900,   2400,   2700,   3100,   3400,   3800,   4100 }; // from Sander (2015)
}