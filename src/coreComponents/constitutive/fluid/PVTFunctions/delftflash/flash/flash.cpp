#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <chrono>

#include "ssi.h"
#include "flash.h"
#include "../phase-correlations/eos.h"
#include "../phase-correlations/phase.h"
#include "../global/misc.h"
#include "../global/global.h"

using namespace std;

void Flash::getIndices(std::vector<std::string> components, std::vector<std::string> phases) {
	ci = std::vector<int>(NC);
	for (int i = 0; i < NC; i++) {
		ci[i] = searchString(Composition::components, Composition::compsize, components[i]); // finds index of each component for correlations
	}
	pi = std::vector<int>(NP);
	for (int j = 0; j < NP; j++) {
		pi[j] = searchString(Composition::phases, Composition::phasesize, phases[j]); // finds index of each phase for correlations
	}
	return;
}

std::vector<std::string> Flash::getPhases() {
	std::vector<std::string> phases(NP);
	for (int j = 0; j < NP; j++) {
		phases[j] = Composition::phases[pi[j]];
	}
	
	return phases;
}

bool SSIFlash::runFlash(double pres, double Temp, std::vector<double> zc, std::vector<std::string> comp, std::vector<std::string> ph, std::vector<int> phase_eos) {
	p = pres; T = Temp; z = zc; eos = phase_eos;
	NC = comp.size(); NP = ph.size();
	getIndices(comp, ph); // finds ci and pi indices for components and phases used in correlations
	
	// Run NP flash
	SSI ssi(p, T, z, ci, pi, eos); // Runs recursive multi-stage SSI
	V = ssi.npFlash();
	NP = V.size();
	x = ssi.getx();
	pi = ssi.getpi();
	std::vector<std::string> phases = getPhases();
	
	return true;
}

bool SSIHydrateFlash::runFlash(double pres, double Temp, std::vector<double> zc, std::vector<std::string> comp, std::vector<std::string> ph, std::vector<int> phase_eos) { // flash with kinetic hydrate
	p = pres; T = Temp; z = zc; eos = phase_eos;
	NC = comp.size(); NP = ph.size();
	getIndices(comp, ph); // finds indices for components and phases used in correlations
	
	// First run regular non-hydrate flash
	SSI ssi(p, T, z, ci, pi, eos); // Runs recursive multi-stage SSI
	V = ssi.npFlash();
	NP = V.size();
	x = ssi.getx();
	pi = ssi.getpi();
	f = ssi.getf();
	water_index = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 0));
	
	std::vector<std::string> phases = getPhases();

	// Then calculate fugacity of water in the hydrate phase 
	H_phases = std::vector<bool>(1, false);
	f_wH = std::vector<double>(1);

	Hydrate sI(ci, 3);
	sI.init(p, T);
	sI.hydrateFugacity(T, f);
	f_wH[0] = sI.getfH();
	x_H = sI.getx();
	
	return true;
}
