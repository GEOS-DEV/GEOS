#pragma once

class EoS
{
public:
	EoS() {};
	virtual void init(double p, double T) = 0;
	virtual std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) = 0;
	virtual double Z(double p, double T, std::vector<double> x) {p = p; T= T; x = x; return 0.;};
};

class PR :public EoS 
{
private:
	int NC, phaseIndex;
	std::vector<int> ci;
	double a, b, A, B;
	std::vector<double> ai, aij, bi;

public:
	PR(std::vector<int> ci_, int phase_index);
	
	void init(double p, double T) override;
	double Z(double p, double T, std::vector<double> x) override;
	std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) override;
};

class SRK :public EoS 
{
private:
	int NC, phaseIndex;
	std::vector<int> ci;
	double a, b, A, B;
	std::vector<double> ai, aij, bi;

public:
	SRK(std::vector<int> ci_, int phase_index);

	void init(double p, double T) override;
	double Z(double p, double T, std::vector<double> x) override;
	std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) override;
};

class AQ1 :public EoS 
{
private:
	int NC;
	std::vector<int> ci;
	double K0_H2O;
	double V_H2O{ 18.1 };
	std::vector<double> k_H, labda_i, ksi_i;
	bool salt{ false };
	double x_salt{ 0 };
public:
	AQ1(std::vector<int> ci_);

	void init(double p, double T) override;

	std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) override;

};

class AQ2 :public EoS 
{
private:
	double pp, TT; // needed for integrals
	int NC, index; // needed for integrals
	std::vector<int> ci;
	std::vector<double> gi; // gibbs energy
	std::vector<double> hi; // enthalpy
	std::vector<double> vi; // volume
	std::vector<double> gi0; // gibbs energy of ideal gas
	bool salt{ false };
	double x_salt{ 0 };
	int water_index, salt_index; // composition index of H2O and salt components
	double B_ca, C_ca, D_ca, eps, A_DH; // ionic contribution coefficients

public:
	AQ2(std::vector<int> ci_);

	void init(double p, double T) override;

private:
	double Aq2_integrals(double a, double b, int steps, int function);
	double h_io(double x);
	double h_i(double x);
	double v_i(double x);
	double xt(double x);
	std::vector<double> ln_a(std::vector<double> m, std::vector<int> ci);

public:
	std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) override;

};

class Hydrate :public EoS
{
private:
	int NC, water_index, pi; // stores type: 0) sI; 1) sII; 2) sH
	double pp, TT;
	std::vector<int> ci; // component indices
	std::vector<double> theta_im;
	std::vector<double> x_H, f;
	double g_w0;
	double g_beta;
	double f_wH;

	// type specific parameters
	int comp_index, cage_index, shell_index; // component, cage and shell indices, needed for integrals
	int R1_index; // index of innermost shell in cage
	std::vector<int> zm; // #waters in cage
	std::vector<int> zn; // #waters in shell
	std::vector<int> n_shells; // #shells in cage
	std::vector<double> Rn; // radius of shells
	std::vector<double> dr_im; // repulsive constants of components in cages
	std::vector<double> Nm; // #guests per unit cell
	int n_cages; // #cages in hydrate structure
	double nH2O; // #H2O molecules in hydrate structure
public:
	Hydrate(std::vector<int> ci_, int phase_index);

	void init(double p, double T) override;
	void hydrateFugacity(double T, std::vector<double> f);
	void xH();
	std::vector<double> &getx() { return x_H; }
	double &getfH() { return f_wH; }

	std::vector<double> fugacityCoefficient(double p, double T, std::vector<double> x) override;
	
private:
	double dmuH(std::vector<double> f);
	double Hyd_integrals(double a, double b, int steps, int function);
	double h_io(double x);
	double h_beta(double x);
	double v_beta(double x);
	double pot(double x);
	double dvH(double x);
};

class InitialK {
private:
	double p, T;
	int NC, water_index;
	std::vector<double> z;
	std::vector<int> ci; // component indices
public:
	InitialK(double p_, double T_, std::vector<double> z_, std::vector<int> ci_);
	std::vector<double> K_initial(std::vector<int> pi, int j);

private:
	std::vector<double> vapour_aqueous();
	std::vector<double> vapour_liquid();
	std::vector<double> vapour_sI();
	std::vector<double> vapour_sII();

};