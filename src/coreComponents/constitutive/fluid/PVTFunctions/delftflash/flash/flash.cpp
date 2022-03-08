#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>

#include "flash.h"
#include "../phase-correlations/phase.h"
#include "../global/global.h"

using namespace std;

Flash::Flash(std::vector<std::string> species, std::unordered_map<std::string, std::string> eos_) {
	// Split components from ions
	separateIons(species);

	eos = eos_;
	f = std::vector<double>(NC);
	for (int i = 0; i < NC; i++) 
	{ 
		if (components[i] == "H2O") 
		{ 
			H2O_idx = i; 
		}
	}
}

void Flash::separateIons(std::vector<std::string> species) {
	// check if ions present in components vector
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

void Flash::normalizeIons(std::vector<double> zc) {
	// If ions are present, normalize mole fractions and store ion mole fractions for molality
	if (NI > 0)
	{
		z = std::vector<double>(NC);
		z_i = std::vector<double>(NI);
		zi = 0;
		for (int i = 0; i < NI; i++)
		{
			z_i[i] = zc[NC + i];
			zi += z_i[i];
		}
		for (int i = 0; i < NC; i++)
		{
			z[i] = zc[i] / (1-zi);
		}
	}
	else
	{
		z = zc;
	}

	return;
}

void Flash::renormalizeIons() {
	// If ions are present, re-normalize the phase fractions and phase compositions to include ions in Aq phase
	if (NI > 0)
	{	
		int NS = NC + NI;
		std::vector<double> x_(NP*NS, 0.);
		std::vector<double> f_(NS, 0.);
		for (int j = 0; j < NP; j++)
		{
			// 1. Regard salt as separate phase, normalize other phase fractions for this
			V[j] *= (1-zi);
			
			// 2. Copy phase composition from flash results
			for (int i = 0; i < NC; i++)
			{
				x_[j*NS + i] = x[j*NC + i];
			}

			// 3. Modify Aq phase V and x
			if (phases[j] == "Aq")
			{
				// 3a. Add salt fraction to Aq phase fraction
				V[j] += zi;
				
				// 3b. Calculate phase mole fraction of salt components
				double xi = 0;
				for (int i = 0; i < NI; i++)
				{
					x_[j*NS + NC + i] = z_i[i] / V[j];
					xi += z_i[i] / V[j];
				}
				// 3c. Normalize other components
				for (int i = 0; i < NC; i++)
				{
					x_[j*NS + i] *= (1-xi);
					f_[i] = f[i];
				}
			}
		}
		x = x_;
		f = f_;
	}
	return;
}

void Flash::calculateMolality() {
	// Calculate molality of ions
	if (NI > 0)
	{
		m_i = std::vector<double>(NI, 0.);
		for (int i = 0; i < NI; i++)
		{
			double z_H2O = V[Aq_idx] * phaseVector[Aq_idx].x_j[H2O_idx] * (1-zi);  // mole fraction of aqueous H2O
			m_i[i] = 55.509 * z_i[i] / z_H2O;
		}
	}

	phaseVector[Aq_idx].eos->setMolality(m_i);

	return;
}

bool Flash::checkConvergence() {
    double tolerance{ 1E-4 };
	double df;
	for (int i = 0; i < NC; i++) 
	{
		for (int j = 0; j < (NP - 1); j++) 
		{
			for (int jj = j + 1; jj < NP; jj++) 
			{
				df = abs((phaseVector[j].f[i] - phaseVector[jj].f[i]) / phaseVector[j].f[i]); // relative fugacity difference
				if (df > tolerance) 
				{
					return false; // some fugacity has not yet converged
				}
			}
		}
	}
	return true; // all fugacities converged within tolerance
}

std::vector<std::string> Flash::checkPositivePhases() {
	int phaseSum = 0;
    for (int j = 0; j < NP; j++) 
	{ 
		if (V[j] > 0) 
		{ 
			phaseSum++; 
		}
	}

    if (phaseSum == NP) 
	{ 
		return phases; 
	}
    else if (phaseSum > 1) 
	{
        std::vector<std::string> positivePhases(phaseSum);
        int k = 0;
        for (int j = 0; j < NP; j++) 
		{
            if (V[j] > 0) 
			{
                positivePhases[k] = phases[j];
                k++;
            }
        }
        return positivePhases;
    }
	else 
	{ 
		return singlePhase(z); 
	}
}

std::vector<std::string> Flash::singlePhase(std::vector<double> z) {
	std::vector<std::string> positivePhase(1);
	for (int j = 0; j < NP; j++) 
	{
		if (V[j] > 0) 
		{
			positivePhase[0] = phases[j];

			calculateMolality();
			std::vector<double> phi = phaseVector[j].eos->fugacityCoefficient(z);
			for (int i = 0; i < NC; i++) 
			{
				phaseVector[j].x_j[i] = z[i];
				phaseVector[j].f[i] = phi[i] * phaseVector[j].x_j[i] * p;
				f[i] = phaseVector[j].f[i];
			}
			x = z;

			if (NI > 0)
			{
				// If ions are present, and single phase Aq (?), include 
				std::vector<double> x_(NC + NI);
				for (int i = 0; i < NC; i++)
				{
					x_[i] = x[i]*(1-zi);
				}
				for (int i = 0; i < NI; i++)
				{
					x_[NC + i] = z_i[i];
				}
				x = x_;
			}
		}
	}
	V = {1.};
	
	return positivePhase;
}
