#include <unordered_map>
#include <memory>
#include "../phase-correlations/phase.h"
#pragma once

class Flash {
	// Abstract Flash class
protected:
	double p, T;
    int NC, NP;
	std::vector<double> z, V, x, f, K;
	std::vector<std::string> components, ions, phases;
	std::unordered_map<std::string, std::string> eos;
	// std::vector<std::unique_ptr<Phase>> phaseVector;
	Phase* phaseVector;

	int NI{ 0 };  // number of ion species
	double zi = 0;  // total ion mole fraction
	std::vector<double> z_i, m_i; // mole fraction and molality of ions
	int Aq_idx, H2O_idx;  // for calculating salt molality

public:
	Flash(std::vector<std::string> species, std::unordered_map<std::string, std::string> eos_);
	
	virtual std::vector<std::string> runFlash(double p_, double T_, std::vector<double> zc) = 0;
	
	std::vector<double> &getV() { return V; }
	std::vector<double> &getx() { return x; }
	std::vector<double> &getf() { return f; }
	std::vector<double> &getMolality() { return m_i; }
	
protected:
	// void normalize0();
	void separateIons(std::vector<std::string> species);
	void normalizeIons(std::vector<double> zc);
	void calculateMolality();
	void renormalizeIons();
	
	virtual void convergenceLoop() = 0;
	bool checkConvergence();
	std::vector<std::string> checkPositivePhases();
	std::vector<std::string> singlePhase(std::vector<double> z);

};

class NPhaseSSI : public Flash 
{
public:
	NPhaseSSI(std::vector<std::string> species, std::unordered_map<std::string, std::string> eos_, std::vector<std::string> ph); 

	std::vector<std::string> runFlash(double p_, double T_, std::vector<double> zc) override;

private:
	virtual void convergenceLoop() override;
};
