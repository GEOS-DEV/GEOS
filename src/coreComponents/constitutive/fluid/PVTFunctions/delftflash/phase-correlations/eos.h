#include <complex>
#include <vector>
#include "../global/global.h"
#pragma once
//#if defined(__GNUC__)
#pragma GCC diagnostic ignored "-Wunused-parameter"
//#endif
class EoS
{
protected:
	double p, T, g;
	int NC;
	std::string phase;
	std::vector<std::string> components;
	
	EoS(std::vector<std::string> comp) 
	{ 
		components = comp; NC = comp.size(); 
	}
public:
	virtual ~EoS() {}
	virtual void parameters(double p_, double T_) = 0;
	virtual double calc_z(std::vector<double> x) { x=x;return 0; }
	virtual double calc_H(double p_, double T_, std::vector<double> x) { return 0; }
	virtual void setMolality(std::vector<double> mi) { return; }
	
	virtual std::vector<double> fugacityCoefficient(std::vector<double> x) = 0;
	// double gibbs(std::vector<double> x, std::vector<double> phi) {
	// 	g = 0;
	// 	for (int i = 0; i < NC; i++) { g += x[i] * log(x[i] * phi[i])}
	// }
};

class TwoParCubic : public EoS
{
protected:
	double R = 8.3145e-5;
	double A, B, omegaA, omegaB, d1, d2;
	std::vector<double> Ai, Bi, Aij, kappa;
    
	TwoParCubic(std::vector<std::string> comp) : EoS(comp) { }
	
	virtual void init() = 0;
	
	std::vector<std::complex<double>> cubic_roots(double a, double b, double c, double d);
public:
	virtual void parameters(double p_, double T_);

	double calc_z(std::vector<double> x) override;
	double calc_H(double p_, double T_, std::vector<double> x) override;

    std::vector<double> fugacityCoefficient(std::vector<double> x) override;
};

class PR : public TwoParCubic
{
private:
	void init() override;
public:
	PR(std::vector<std::string> comp) : TwoParCubic(comp) { 
		init(); 
	}
};

class SRK: public TwoParCubic
{
private:
	void init() override;
public:
	SRK(std::vector<std::string> comp) : TwoParCubic(comp) { 
		init(); 
	}
};

class AQ : public EoS
{
protected:
	int NI;
	int water_index; // composition index of H2O and salt components;
	std::vector<std::string> ions;
	std::vector<double> m_s, m_i;

	AQ(std::vector<std::string> comp, std::vector<std::string> ions_) : EoS(comp)
	{
		ions = ions_; NI = ions_.size();
	}

public:
	void setMolality(std::vector<double> mi) override { m_i = mi; }
};

class AQ1 : public AQ
{
private:
	double V_H2O = 18.1;
	double K0_H2O;
	double m_a{ 0. }, m_c{ 0. }, m_ac{ 0. }, xs{ 0. };
	std::vector<double> k_H, labda, ksi;

public:
	AQ1(std::vector<std::string> comp, std::vector<std::string> ions) : AQ(comp, ions) { }
	
	void parameters(double p_, double T_) override;
	std::vector<double> fugacityCoefficient(std::vector<double> x) override;
};

class AQ2 :public AQ
{
private:
	double pp, TT; // needed for integrals
	std::string comp; // needed for integrals
	std::vector<double> gi, hi, vi, gi0; // gibbs energy, enthalpy, volume, gibbs energy of ideal gas
	double eps, A_DH; // ionic contribution coefficients
	std::vector<double> gamma_P1, gamma_P2, B_ca, C_ca, D_ca;
	
	int NS;
	std::vector<std::string> species;

public:
	AQ2(std::vector<std::string> comp, std::vector<std::string> ions) : AQ(comp, ions) 
	{
		species = comp;
		if (ions.size() > 0)
		{
			species.insert(species.end(), ions.begin(), ions.end());
		}
		NS = species.size();
	}
	
	void parameters(double p_, double T_) override;
	std::vector<double> fugacityCoefficient(std::vector<double> x) override;

private:
	double Aq2_integrals(double a, double b, int steps, int function);
	double Aq2_functions(double x, int function);
	double xt(double x);
	std::vector<double> ln_a(std::vector<double> m);

};

class AQ3 : public AQ
{
private:
	// AQ2 parameters for H2O fugacity
	double pp, TT; // needed for integrals
	std::string comp; // needed for integrals
	double giw, hiw, viw, giw0; // H2O gibbs energy, enthalpy, volume, gibbs energy of ideal gas
	double eps, A_DH; // ionic contribution coefficients
	std::vector<double> gamma_P1, gamma_P2, B_ca, C_ca, D_ca;

	int NS;
	std::vector<std::string> species;

	// AQ1 parameters for solute fugacity
	double m_a{ 0. }, m_c{ 0. }, m_ac{ 0. }, xs{ 0. };
	std::vector<double> k_H, labda, ksi;

public:
	AQ3(std::vector<std::string> comp, std::vector<std::string> ions) : AQ(comp, ions) 
	{
		species = comp;
		if (ions.size() > 0)
		{
			species.insert(species.end(), ions.begin(), ions.end());
		}
		NS = species.size();
	}
	
	void parameters(double p_, double T_) override;
	std::vector<double> fugacityCoefficient(std::vector<double> x) override;

private:
	double Aq2_integrals(double a, double b, int steps, int function);
	double Aq2_functions(double x, int function);
	double xt(double x);
	double ln_a(std::vector<double> m);

};

class VdWP :public EoS
{
private:
	int water_index;
	double pp, TT;
	std::vector<double> theta_im;
	double g_w0, g_beta;

	// member variables needed for integrals
	std::string comp;
	int cage_index, shell_index, R1_index; // cage, shell indices and index of innermost shell in cage, needed for integrals
	std::vector<int> zm, zn, n_shells; // #waters in cage, #waters in shell, #shells in cage
	std::vector<double> Rn, Nm; // radius of shells, #guests per unit cell
	int n_cages; // #cages in hydrate structure
	double nH2O; // #H2O molecules in hydrate structure

public:
	VdWP(std::vector<std::string> comp_, std::string hydrate_type) : EoS(comp_) 
	{ 
		phase = hydrate_type; 
	}
	void parameters(double p_, double T_) override;

	std::vector<double> xH();
	double fwH(std::vector<double> f0);

	std::vector<double> fugacityCoefficient(std::vector<double> f) override;
	
private:
	double dmuH(std::vector<double> f);
	double hydrate_integrals(double a, double b, int steps, int function);
	double hydrate_functions(double x, int function);
	double pot(double x);
	double dvH(double x);
	
};

class InitialK {
private:
	double p, T;
	int NC, water_index;
	std::vector<double> z;
	std::vector<std::string> components;
public:
	InitialK(double p_, double T_, std::vector<double> zc, std::vector<std::string> comp);
	std::vector<double> K_initial(std::vector<std::string> phases, int j);

private:
	std::vector<double> vapour_aqueous();
	std::vector<double> vapour_liquid();
	std::vector<double> vapour_sI();
	std::vector<double> vapour_sII();

};