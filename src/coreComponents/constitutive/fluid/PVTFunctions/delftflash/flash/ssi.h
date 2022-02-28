#pragma once

class Phase;

class SSI {
protected:
	double p, T;
	int NC, NP;
	std::vector<double> z, V, x, f;
	std::vector<int> ci, pi, eos;
	Phase* phaseVector;

public:
	SSI(double p_, double T_, std::vector<double> z_, std::vector<int> ci_, std::vector<int> pi_, std::vector<int> eos_);
	std::vector<double> npFlash();
	// Phase* &getPhaseVector() { return phaseVector; }
	std::vector<double> &getx() { return x; }
	std::vector<int> &getpi() { return pi; }
	std::vector<double> &getf() { return f; }

protected:
	std::vector<double> runSSI();
	virtual void convergenceLoop(std::vector<double> K);
	bool checkConvergence();
	std::vector<int> checkPositivePhases();
	void singlePhase();
		
};
