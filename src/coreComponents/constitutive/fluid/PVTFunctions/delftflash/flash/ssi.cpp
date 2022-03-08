#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "flash.h"
#include "rr.h"
#include "../phase-correlations/eos.h"
#include "../phase-correlations/phase.h"

using namespace std;

NPhaseSSI::NPhaseSSI(std::vector<std::string> species, std::unordered_map<std::string, std::string> eos_, std::vector<std::string> ph) : Flash(species, eos_) { 
	phases = ph;
	NP = phases.size();
	// Each phase is stored as Phase class in phaseVector
	phaseVector = new Phase[NP];
	for (int j = 0; j < NP; j++) 
    {
		phaseVector[j].phase = ph[j];
		phaseVector[j].x_j = std::vector<double>(NC);
		phaseVector[j].f = std::vector<double>(NC);
		if (ph[j] == "Aq") 
		{ 
			Aq_idx = j;
			phaseVector[j].ions = ions;
		}
		phaseVector[j].eos = phaseVector[j].setEoS(components, eos_[ph[j]]);
	}
}

std::vector<std::string> NPhaseSSI::runFlash(double p_, double T_, std::vector<double> zc) {
	// Run NP flash
	p = p_; T = T_;
    normalizeIons(zc);
	K = std::vector<double>((NP-1)*NC);

    // Initial K-values
	InitialK initial_k(p, T, zc, components);
    std::vector<double> k(NC);
    for (int j = 1; j < NP; j++) 
    {
        k = initial_k.K_initial(phases, j);
        for (int i = 0; i < NC; i++) 
        {
            K[(j-1)*NC + i] = k[i];
        }
    }
    // Initialize EoS parameters
    for (int j = 0; j < NP; j++) 
    {
		phaseVector[j].eos->parameters(p, T);
    }

    // Converge SSI
	convergenceLoop();

    // Include ions in flash results
    renormalizeIons();
	
	return checkPositivePhases();
}

void NPhaseSSI::convergenceLoop() {
	RRbn rr(z, K);
    bool converged = false;
    std::vector<double> phi0, phi;

    // SSI convergence loop: update fugacity coefficients after each iteration
    int iter = 0;
    while (!converged)
    {
        V = rr.solveRR(K);
        x = rr.getx(V);
	    for (int j = 0; j < NP; j++) 
        {
		    for (int i = 0; i < NC; i++) 
            {
			    phaseVector[j].x_j[i] = x[j*NC + i];
		    }
        }
        calculateMolality();

        // Update fugacity coefficients
        phi0 = phaseVector[0].eos->fugacityCoefficient(phaseVector[0].x_j);
	    for (int i = 0; i < NC; i++) 
        { 
            f[i] = phi0[i] * phaseVector[0].x_j[i] * p; 
            phaseVector[0].f[i] = f[i];
        }

        for (int j = 1; j < NP; j++) 
        {
            phi = phaseVector[j].eos->fugacityCoefficient(phaseVector[j].x_j);
            for (int i = 0; i < NC; i++) 
            {
                phaseVector[j].f[i] = phi[i] * phaseVector[j].x_j[i] * p;
                K[(j-1)*NC + i] = phi0[i] / phi[i];
            }
        }
        converged = checkConvergence();
        
        iter++;
        if (iter > 100)
        {
            converged = true;
            cout << "Flash unable to converge\n";
        }
    }
}
