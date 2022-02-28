#pragma once

class EoS;

class Phase {
public:
	int phase_index; // index for fugacity correlations
	std::vector<double> x_j; // phase composition
	std::vector<double> f;  // component fugacity in phase j
	EoS* eos; // pointer to EoS base class, which points to derived class (specific EoS) later

public: // protected: rather not have fully public class
	EoS* setEoS(std::vector<int> ci, int eos);

};

class PhaseProperties {
private:
	int NC, phaseIndex;
	int water_index;
	std::vector<int> ci;
	std::vector<double> Mw;
	EoS* eos;

public:
	PhaseProperties(std::vector<std::string> comp, std::string phase, int eos);
	double molar_weight(std::vector<double> x);
	double density(double p, double T, std::vector<double> x);
	double viscosity(double p, double T, std::vector<double> x);
	// double enthalpy(double p, double T, std::vector<double> x);
	// double heat_capacity();
	// double thermal_conductivity();
	// std::vector<double> fugacity(double p, double T, std::vector<double> x);
private:
	EoS* setEoS(std::vector<int> ci, int phase_index, int eos);

};