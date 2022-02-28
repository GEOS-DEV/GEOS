#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>

#include "ssi.h"
#include "rr.h"
#include "../phase-correlations/eos.h"
#include "../phase-correlations/phase.h"
#include "../global/misc.h"
#include "../global/global.h"

using namespace std;

SSI::SSI(double p_, double T_, std::vector<double> z_, std::vector<int> ci_, std::vector<int> pi_, std::vector<int> eos_) {
	// Initialize member variables
	p = p_; T = T_; z = z_;
	ci = ci_; NC = ci.size();
	pi = pi_; NP = pi.size();
	
	eos = eos_;
	f = std::vector<double>(NC);

	// Each phase is stored as Phase class in phaseVector
	phaseVector = new Phase[NP];
	for (int j = 0; j < NP; j++) {
		phaseVector[j].phase_index = pi[j];
		phaseVector[j].x_j = std::vector<double>(NC);
		phaseVector[j].f = std::vector<double>(NC);
	}
}

std::vector<double> SSI::npFlash() {
	// First try NP flash
	V = runSSI();
	int countPhases = 0;
	for (int j = 0; j < NP; j++) { if (V[j] > 0) { countPhases++; } }

    if (countPhases == NP) { // NP-phase behaviour confirmed
		return V;
	} else if (NP == 2 && countPhases == 1) { // single-phase behaviour
		singlePhase();
        return V;
	} else {
		cout << "SSI Flash did not converge.\n";
		return V;
	}
}

std::vector<double> SSI::runSSI() {
    std::vector<double> K((NP-1)*NC);

    // Initial K-values
	InitialK initial_k(p, T, z, ci);
    std::vector<double> k(NC);
    for (int j = 1; j < NP; j++) {
        k = initial_k.K_initial(pi, j);
        for (int i = 0; i < NC; i++) {
            K[(j-1)*NC + i] = k[i];
        }
    }
    // Initialize EoS
    for (int j = 0; j < NP; j++) {
        phaseVector[j].eos = phaseVector[j].setEoS(ci, eos[j]);
		phaseVector[j].eos->init(p, T);
    }

    // Converge SSI
	convergenceLoop(K);
    
    return V;
}

void SSI::convergenceLoop(std::vector<double> K) {
	RRbn rr(z, K);
    V = rr.solveRR(K);
    x = rr.getx(V);
	for (int j = 0; j < NP; j++) {
		for (int i = 0; i < NC; i++) {
			phaseVector[j].x_j[i] = x[j*NC + i];
		}
	}

    // Update fugacity coefficients
	std::vector<double> phi0, phi;
    phi0 = phaseVector[0].eos->fugacityCoefficient(p, T, phaseVector[0].x_j);
	for (int i = 0; i < NC; i++) { phaseVector[0].f[i] = phi0[i] * phaseVector[0].x_j[i] * p; }

    for (int j = 1; j < NP; j++) {
        phi = phaseVector[j].eos->fugacityCoefficient(p, T, phaseVector[j].x_j);
        for (int i = 0; i < NC; i++) {
            phaseVector[j].f[i] = phi[i] * phaseVector[j].x_j[i] * p;
            K[(j-1)*NC + i] = phi0[i] / phi[i];
        }
    }
    bool converged = checkConvergence();

    // SSI convergence loop
    while (!converged) {
        // Solve RR for updated K-values
        V = rr.solveRR(K);
        x = rr.getx(V);
		for (int j = 0; j < NP; j++) {
			for (int i = 0; i < NC; i++) {
				phaseVector[j].x_j[i] = x[j*NC + i];
			}
		}

        // Update fugacity coefficients
		phi0 = phaseVector[0].eos->fugacityCoefficient(p, T, phaseVector[0].x_j);
		for (int i = 0; i < NC; i++) { phaseVector[0].f[i] = phi0[i] * phaseVector[0].x_j[i] * p; }

		for (int j = 1; j < NP; j++) {
			phi = phaseVector[j].eos->fugacityCoefficient(p, T, phaseVector[j].x_j);
			for (int i = 0; i < NC; i++) {
				phaseVector[j].f[i] = phi[i] * phaseVector[j].x_j[i] * p;
				K[(j-1)*NC + i] = phi0[i] / phi[i];
				f[i] = phaseVector[j].f[i];
			}
		}
        converged = checkConvergence();
    }
}

bool SSI::checkConvergence() {
    double tolerance{ 1E-4 };
	double df;
	for (int i = 0; i < NC; i++) {
		for (int j = 0; j < (NP - 1); j++) {
			for (int jj = j + 1; jj < NP; jj++) {
				df = abs((phaseVector[j].f[i] - phaseVector[jj].f[i]) / phaseVector[j].f[i]); // relative fugacity difference
				if (df > tolerance) {
					return false; // some fugacity has not yet converged
				}
			}
		}
	}
	return true; // all fugacities converged within tolerance
}

std::vector<int> SSI::checkPositivePhases() {
	int phaseSum = 0;
    for (int j = 0; j < NP; j++) { if (V[j] > 0) { phaseSum++; } }

    if (phaseSum == NP) { return pi; }
    else {
        std::vector<int> positivePhases(phaseSum);
		// std::vector<int> eos_(phaseSum);
        int k = 0;
        for (int j = 0; j < NP; j++) {
            if (V[j] > 0) {
                positivePhases[k] = pi[j];
				// eos_[k] = eos[j];
                k++;
            }
        }
		// eos = eos_;
        return positivePhases;
    }
}

void SSI::singlePhase() {
	for (int j = 0; j < NP; j++) {
		if (V[j] > 0) {
			std::vector<double> phi = phaseVector[j].eos->fugacityCoefficient(p, T, z);
			for (int i = 0; i < NC; i++) {
				phaseVector[j].x_j[i] = z[i];
				phaseVector[j].f[i] = phi[i] * phaseVector[j].x_j[i] * p;
				f[i] = phaseVector[j].f[i];
			}
			x = phaseVector[j].x_j;
		}
	}
	pi = checkPositivePhases();
	V = {1.};
	
	return;

}
