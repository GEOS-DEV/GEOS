#include <vector>
#include "eos.h"
#include "../global/global.h"
#pragma once
//#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-function"
//#endif
class Phase {
public:
	std::string phase;
	std::vector<double> x_j; // phase composition
	std::vector<double> f;  // component fugacity in phase j
	std::vector<std::string> ions;
	EoS* eos; // pointer to EoS base class, which points to derived class (specific EoS) later

public: // protected: rather not have fully public class
	EoS* setEoS(std::vector<std::string> components, std::string eos_);

};

class Properties : public Phase
{
protected:
	int NC, water_index;
	bool H2O{ false };
	std::vector<double> m_i;  // molality [mol/kg]
	std::vector<std::string> components;
	std::vector<double> Mw;

public:
	Properties(std::vector<std::string> comp, std::string ph, std::string eos_);

protected:
	double density(double p, double T, std::vector<double> x, double m=0.);
	double viscosity(double p, double T, std::vector<double> x, double rho, double m=0.);
	double enthalpy(double p, double T, std::vector<double> x, double m=0.);
	double thermal_conductivity(double p, double T, std::vector<double> x, double rho, double m=0.);
	std::vector<double> fugacity(double p, double T, std::vector<double> x, double m=0.);
};

class Density : public Properties
{
public:
	Density(std::vector<std::string> comp, std::string ph, std::string eos_) : Properties(comp, ph, eos_) {

	}
	double evaluate(double p, double T, std::vector<double> x, double m=0.) { return density(p, T, x, m); };
};

class Viscosity : public Properties
{
public:
	Viscosity(std::vector<std::string> comp, std::string ph, std::string eos_) : Properties(comp, ph, eos_) {

	}
	double evaluate(double p, double T, std::vector<double> x, double rho, double m=0.) { return viscosity(p, T, x, rho, m); };
};

class Enthalpy : public Properties
{
public:
	Enthalpy(std::vector<std::string> comp, std::string ph, std::string eos_) : Properties(comp, ph, eos_) {

	}
	double evaluate(double p, double T, std::vector<double> x, double m=0.) { return enthalpy(p, T, x, m); };
};

class ThermalConductivity : public Properties
{
public:
	ThermalConductivity(std::vector<std::string> comp, std::string ph, std::string eos_) : Properties(comp, ph, eos_) {	}
	double evaluate(double p, double T, std::vector<double> x, double rho, double m=0.) { return thermal_conductivity(p, T, x, rho, m); };
};

class KineticHydrate
{
private:
	int NC, NI;
	std::vector<std::string> components, ions;
	VdWP* h_eos;

public:
	KineticHydrate(std::vector<std::string> species, std::string ph) {
		separateIons(species);
		h_eos = new VdWP(components, ph);
	}
	double fwH(double p, double T, std::vector<double> f0);
	std::vector<double> xH() { return h_eos->xH(); };

private:
	void separateIons(std::vector<std::string> species);
};

class ComponentProp
{
private:
	int NC;
	std::vector<std::string> components;
public:
	ComponentProp(std::vector<std::string> comp) { components = comp; NC = comp.size(); }
	double mix_MW(std::vector<double> x);
	// std::vector<double> get(std::string prop);
	// void set(std::string prop, std::vector<double>);
};
