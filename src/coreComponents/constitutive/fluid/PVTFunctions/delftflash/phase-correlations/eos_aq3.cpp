#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <unordered_map>

#include "eos.h"
#include "../global/global.h"

using namespace std;

void AQ3::parameters(double p_, double T_) {
    p = p_; pp = p_; 
    T = T_; TT = T_; // initialize for reference p & T in integrals

    water_index = std::distance(components.begin(), std::find(components.begin(), components.end(), "H2O"));

    comp = "H2O";
    // ideal gas Gibbs energy
    double hio = Aq2_integrals(aq2_par::T_0, T, 10, 0);  // eq. 3.3
    double gio = aq2_par::gi_00["H2O"]/(aq2_par::R * aq2_par::T_0);
    giw0 = gio - hio;  // eq. 3.2

    // Gibbs energy of aqueous phase
    giw = aq2_par::gi_0["H2O"]/(aq2_par::R*aq2_par::T_0);
    hiw = Aq2_integrals(aq2_par::T_0, T, 10, 1);
    viw = Aq2_integrals(aq2_par::P_0, p, 10, 2);

    // Initialize vectors
    gamma_P1 = std::vector<double>(NS*NS, 0.);
    gamma_P2 = std::vector<double>(NS*NS, 0.);

    // salt related numbers
    B_ca = std::vector<double>(NI*NI, 0.);
    C_ca = std::vector<double>(NI*NI, 0.);
    D_ca = std::vector<double>(NI*NI, 0.);
    for (int j = 0; j < NI; j++)
    {
        int cj = CompProp::charge[ions[j]];
        if (cj > 0)  // only cations (+ charge)
        {
            for (int k = 0; k < NI; k++)
            {
                int ck = CompProp::charge[ions[k]];
                if (ck < 0)  // only anions (- charge)
                {
                    B_ca[j*NI + k] = aq2_par::Bca[ions[j]][ions[k]][0] + aq2_par::Bca[ions[j]][ions[k]][1] * TT + aq2_par::Bca[ions[j]][ions[k]][2] * pow(TT, 2);
                    C_ca[j*NI + k] = aq2_par::Cca[ions[j]][ions[k]][0] + aq2_par::Cca[ions[j]][ions[k]][1] * TT + aq2_par::Cca[ions[j]][ions[k]][2] * pow(TT, 2);
                    D_ca[j*NI + k] = aq2_par::Dca[ions[j]][ions[k]][0] + aq2_par::Dca[ions[j]][ions[k]][1] * TT + aq2_par::Dca[ions[j]][ions[k]][2] * pow(TT, 2);
                }
            }
        }
    }
    
    eps = 0;
    for (int i = 0; i < 3; i++)
    {   
        double c = aq2_par::eps[i*3 + 0] + aq2_par::eps[i*3 + 1] * pp + aq2_par::eps[i*3 + 2] * pow(pp, 2);
        eps += c * pow(TT, i);
    }
    double rho_s = 1000;
    A_DH = pow(pow(aq2_par::e, 2)/(aq2_par::eps0*eps*aq2_par::R*TT), 1.5) * pow(M_NA, 2)/(8*M_PI) * sqrt(2*rho_s); // kg^0.5/mol^0.5

    // Calculate composition-independent parameters of Aq-EoS
	double tau = 1 - T / CompProp::prop["Tc"]["H2O"];
	double P_s = CompProp::prop["Pc"]["H2O"] * exp(CompProp::prop["Tc"]["H2O"] / T * (tau * -7.85951783 + 1.84408259 * pow(tau, 1.5) - 11.7866497 * pow(tau, 3)
			+ 22.6807411 * pow(tau, 3.5) - 15.9618719 * pow(tau, 4.0) + 1.80122502 * pow(tau, 7.5)));
	double tc = T - 273.15;  // temperature in Celsius
	double V0 = (1 + 18.159725E-3 * tc) / (0.9998396 + 18.224944E-3 * tc - 7.922210E-6 * pow(tc, 2) - 55.44846E-9 * pow(tc, 3) + 149.7562E-12 * pow(tc, 4) - 393.2952E-15 * pow(tc, 5));
	double B = 19654.320 + 147.037 * tc - 2.21554 * pow(tc, 2) + 1.0478E-2 * pow(tc, 3) - 2.2789E-5 * pow(tc, 4);
	double A1 = 3.2891 - 2.3910E-3 * tc + 2.8446E-4 * pow(tc, 2) - 2.8200E-6 * pow(tc, 3) + 8.477E-9 * pow(tc, 4);
	double A2 = 6.245E-5 - 3.913E-6 * tc - 3.499E-8 * pow(tc, 2) + 7.942E-10 * pow(tc, 3) - 3.299E-12 * pow(tc, 4);
	double V = V0 - V0 * p / (B + A1 * p + A2 * pow(p, 2));  // volume of pure water at p[cm3 / g]
	double M = 18.0152;
	double f0_H2O = P_s * exp((p - P_s) * M * V / (aq1_par::R * T));  // fugacity of pure water[bar]
	double rho0_H2O = 1 / V;  // density of pure water at p[g / cm3]
	// double logK0_H2O = -2.209 + 3.097E-2 * tc - 1.098E-4 * pow(tc, 2) + 2.048E-7 * pow(tc, 3);
	// K0_H2O = pow(10., logK0_H2O);  // equilibrium constant for H2O at 1 bar

	k_H = std::vector<double>(NC);
	labda = std::vector<double>(NC);
	ksi = std::vector<double>(NC);
	for (int i = 0; i < NC; i++) 
	{
		comp = components[i];
		if (comp != "H2O")  // non-H2O components
		{
			double dB = aq1_par::tau[comp] + aq1_par::Gamma[comp] * p + aq1_par::beta[comp] * sqrt(1000 / T);
			k_H[i] = exp((1 - aq1_par::eta[comp]) * log(f0_H2O) + aq1_par::eta[comp] * log(aq1_par::R * T / M * rho0_H2O) + 2 * rho0_H2O * dB);
			labda[i] = aq1_par::labda[comp][0] + aq1_par::labda[comp][1] * T + aq1_par::labda[comp][2] / T + aq1_par::labda[comp][3] * p + aq1_par::labda[comp][4] / p + \
				aq1_par::labda[comp][5] * p / T + aq1_par::labda[comp][6] * T / pow(p, 2) + aq1_par::labda[comp][7] * p / (630 - T) + aq1_par::labda[comp][8] * T * log(p) \
				+ aq1_par::labda[comp][9] * p / pow(T, 2) + aq1_par::labda[comp][10] * pow(p, 2) * T + aq1_par::labda[comp][11] * p * T;
			ksi[i] = aq1_par::ksi[comp][0] + aq1_par::ksi[comp][1] * T + aq1_par::ksi[comp][2] * p / T + aq1_par::ksi[comp][3] * p / (630.0 - T);
		}
	}
	return;
}

double AQ3::xt(double x) {
    double eps_ = 0;
    std::vector<double> c(3);
    for (int i = 0; i < 3; i++)
    {   
        c[i] = aq2_par::eps[i*3 + 0] + aq2_par::eps[i*3 + 1] * pp + aq2_par::eps[i*3 + 2] * pow(pp, 2);
        eps_ += c[i] * pow(x, i);
    }
    double dln = (c[1] + 2 * c[2] * x) / (c[0] + c[1] * x + c[2] * pow(x, 2));
    double d2ln = ((c[0] + c[1] * x + c[2] * pow(x, 2))*c[2] - pow(c[1] + 2*c[2]*x, 2)) / (pow(c[0] + c[1] * x + c[2] * pow(x, 2), 2));
    double dln2 = pow(dln, 2);
    return 1/eps_ * (d2ln - dln2) * x;
}

double AQ3::Aq2_functions(double x, int function) {
    if (function == 0)  // hio
    {
        double hi0 = aq2_par::hi_00[comp] + aq2_par::hi_a[comp][0]*(x-aq2_par::T_0) + 1/2.0 * aq2_par::hi_a[comp][1]*(pow(x, 2)-pow(aq2_par::T_0, 2)) + 1/3.0 * aq2_par::hi_a[comp][2]*(pow(x, 3)-pow(aq2_par::T_0, 3)) + 1/4.0 * aq2_par::hi_a[comp][3]*(pow(x, 4)-pow(aq2_par::T_0, 4));
        return hi0 / (aq2_par::R * pow(x, 2));
    } 
    else if (function == 1)  // hi
    {
        if (comp == "H2O")  // pure water enthalpy
        {
            double hw = aq2_par::hi_0[comp] + 8.712 * aq2_par::R * (x - aq2_par::T_0) + 1 / 2.0 * 0.125E-2 * aq2_par::R * (pow(x, 2) - pow(aq2_par::T_0, 2)) - 1 / 3.0 * 0.018E-5 * aq2_par::R * (pow(x, 3) - pow(aq2_par::T_0, 3));
            return hw / (aq2_par::R * pow(x, 2));
        }
        else  // solute molar enthalpy
        {
            double xt = Aq2_integrals(aq2_par::T_0, x, 20, 3);
            double int_cp = aq2_par::cp1[comp]*(x-aq2_par::T_0) - aq2_par::cp1[comp]*(1/(x-228.)-1/(aq2_par::T_0-228.)) + aq2_par::omega[comp]*xt;
            return (aq2_par::hi_0[comp] + int_cp)/(aq2_par::R * pow(x, 2));
        }
    } 
    else  // vi
    {
        if (comp == "H2O") 
        {
            double vw = ((31.1251 - 2.46176E-2 * x + 8.69425E-6 * pow(x, 2) - 6.03348E-10 * pow(x, 3)) +
                        (-1.14154E-1 + 2.15663E-4 * x - 7.96939E-8 * pow(x, 2) + 5.57791E-12 * pow(x, 3)) * TT +
                        (3.10034E-4 - 6.48160E-7 * x + 2.45391E-10 * pow(x, 2) - 1.72577E-14 * pow(x, 3)) * pow(TT, 2) +
                        (-2.48318E-7 + 6.47521E-10 * x - 2.51773E-13 * pow(x, 2) + 1.77978E-17 * pow(x, 3)) * pow(TT, 3)) * 1E-6;  // m3/mol
            return vw / (aq2_par::R * 1E-5 * TT);  // m3/mol / (m3.bar/K.mol). K
        }
        else 
        {
            double eps_ = 0;
            std::vector<double> b(3);
            for (int i = 0; i < 3; i++)
            {   
                b[i] = aq2_par::eps[0*3 + i] + aq2_par::eps[1*3 + i] * TT + aq2_par::eps[2*3 + i] * pow(TT, 2);
                eps_ += b[i] * pow(x, i);
            }
            double dedp = b[1] + 2*x*b[2];
            double vi_ = aq2_par::vp[comp][0] + aq2_par::vp[comp][1]/(2600+x) + (aq2_par::vp[comp][2] + aq2_par::vp[comp][3]/(2600+x))*1/(TT-228) - aq2_par::omega[comp]/(pow(eps_, 2)) * dedp;  // J/mol.bar
            return vi_ / (aq2_par::R * TT);  // J/mol.bar / (J/mol.K).K
        }
    }
}

double AQ3::Aq2_integrals(double a, double b, int steps, int function) {
    // integrals solved numerically with simpson's rule
    double s = 0;
    double h = (b-a)/steps;
    double x;

    switch (function)
    {
    case 0: // ideal gas enthalpy
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = T
            s += h*((Aq2_functions(x, 0) + 4*Aq2_functions(x+h/2, 0) + Aq2_functions(x+h, 0))/6); // hio
        };
        break;
    case 1: // hi
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = T
            s += h*((Aq2_functions(x, 1) + 4*Aq2_functions(x+h/2, 1) + Aq2_functions(x+h, 1))/6); // hi
        };
        break;
    case 2: // vi
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = p
            s += h*((Aq2_functions(x, 2) + 4*Aq2_functions(x+h/2, 2) + Aq2_functions(x+h, 2))/6); // vi
        };
        break;
    case 3: // XT integral (for h_i)
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = T
            s += h*((xt(x) + 4*xt(x+h/2) + xt(x+h))/6);
        }
        break;
    }
    return s;
}

double AQ3::ln_a(std::vector<double> m_s) {
    // Activity of water and solutes
    double lna = 0.;
    
    double I = 0.; // ionic strength
    // Calculate ionic strength        
    for (int i = 0; i < NI; i++)
    {
        I += 0.5 * pow(CompProp::charge[ions[i]], 2)*m_s[NC + i];
    }
    
    // Calculate solute interactions (m-m or m-i or i-i)
    // gamma_P1 = std::unordered_map<std::string, std::unordered_map<std::string, double>>;
    for (int j = 0; j < NS; j++) 
    {
        // Pitzer
        if (aq2_par::B.find(species[j]) != aq2_par::B.end())
        {
            for (int k = 0; k < NS; k++)
            {
                if (aq2_par::B[species[j]].find(species[k]) != aq2_par::B[species[j]].end())
                {
                    std::vector<double> b = aq2_par::B[species[j]][species[k]];
                    gamma_P1[j*NS + k] = b[0] + b[1] * T + b[2] * p;  // B0
                    if (NI > 0)
                    {
                        // Pitzer contribution of ions
                        double B1 = b[3];
                        gamma_P1[j*NS + k] += B1/(2*I) * (1 - (1 + 2*sqrt(I)) * exp(-2*sqrt(I)));
                        gamma_P2[j*NS + k] = B1 * (1 - (1 + 2*sqrt(I) + 2*I) * exp(-2*sqrt(I)));
                    }
                }
            }
        }
    }
    
    // Calculate activity of components
    // int i = water_index;
    if (NI > 0)
    {
        // Ionic contribution to water activity [eq. 3.22]
        double numerator = 0;
        double c_denominator = 0;
        double a_denominator = 0;
        for (int j = 0; j < NI; j++)
        {   
            // numerator
            int cj = CompProp::charge[ions[j]];
            if (cj > 0)  // only cations (+ charge)
            {
                for (int k = 0; k < NI; k++)
                {
                    int ck = CompProp::charge[ions[k]];
                    if (ck < 0)  // only anions (- charge)
                    {
                        // eq. 3.25
                        double cc_ca = fabs(cj*ck);  // |charge_j * charge_k|
                        double gamma_Bca = (0.13816 + 0.6*B_ca[j*NI + k])*I*cc_ca/1.5;
                        gamma_Bca *= (1+3*I/cc_ca)/pow(1+3*I/(2*cc_ca), 2) - log(1+3*I/(2*cc_ca))/(3*I/(2*cc_ca));
                        gamma_Bca += 2/cc_ca * (0.5*B_ca[j*NI + k]*pow(I, 2) + 2./3.*C_ca[j*NI + k]*pow(I, 3) + 3./4.*D_ca[j*NI + k]*pow(I, 4));
                        // eq. 3.24
                        double gamma_DH = 2*A_DH * ((1-pow(1+sqrt(I), 2))/(1+sqrt(I)) + 2*log(1+sqrt(I)));
                        // eq. 3.23
                        double gamma_Ica = -CompProp::prop["Mw"]["H2O"]/1000*(2*I/cc_ca + gamma_DH + gamma_Bca);
                        numerator += m_s[NC + j] * m_s[NC + k] * pow(cj, 2) * pow(ck, 2) * gamma_Ica;
                    }
                }
            }
            // denominator
            if (cj > 0)
            {
                c_denominator += m_s[NC + j] * pow(cj, 2);
            }
            else
            {
                a_denominator += m_s[NC + j] * pow(cj, 2);
            }
        }
        lna += numerator/(c_denominator*a_denominator);
    }

    // Contribution of molecular/ionic species to water activity [eq. 3.22]
    double molecular{ 0. };
    for (int j = 0; j < NS; j++) 
    { 
        // all molecular/ionic species
        for (int k = 0; k < NS; k++)
        {
            double term = gamma_P1[j*NS + k];
            if (NI > 0)  // P2 contribution only non-zero if ions are present
            {
                term -= gamma_P2[j*NS + k]/(2*I);
            }
            molecular += m_s[j] * m_s[k] * term;
        }
    }
            
    for (int j = 0; j < NC; j++)
    {
        // all molecular species
        molecular += m_s[j];
    }
    lna -= CompProp::prop["Mw"]["H2O"]/1000*molecular;
    
    return lna;
}

std::vector<double> AQ3::fugacityCoefficient(std::vector<double> x) {
    // Construct fugacity coefficients
    std::vector<double> phi(NC);

    // Calculate molality of molecular species
    std::vector<double> m_s(NS);
    for (int i = 0; i < NC; i++) 
    {
        if (species[i] == "H2O")
        {
            m_s[i] = 0.;
        }
        else
        {
            m_s[i] = 55.509 * x[i] / x[water_index];
        }
    }
    // molality of ionic species
    for (int i = 0; i < NI; i++)
    {
        m_s[i+NC] = m_i[i];
    }

    // Find effective salt molalities
	if (NI > 0)
	{	
		m_a = 0.;
		m_c = 0.;
		m_ac = 0.;
		for (int i = 0; i < NI; i++)
		{	
			if (CompProp::charge[ions[i]] > 0)
			{
				m_c += m_i[i] * CompProp::charge[ions[i]];
				m_ac += m_i[i];
			}
			else
			{
				m_a += m_i[i];
			}
		}
		m_ac *= m_a;
    }

    double lna = ln_a(m_s);
    
    for (int i = 0; i < NC; i++) 
    {
        if (components[i] == "H2O")
        {
            double mu = giw - hiw + viw + lna;  // eq. 3.5/3.10
            phi[i] = 1 * exp(mu - giw0)/(x[i]*p);
        }
        else
        {
            double gamma = exp((2 * m_c * labda[i]) + (m_ac * ksi[i]));
			phi[i] = k_H[i] * gamma / p;
        }
    }
    
    return phi;
}