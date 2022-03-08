#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <numeric>
#include <unordered_map>

#include "eos.h"
#include "../global/global.h"

namespace hyd_par {
    // Following hydrate fugacity model described in Ballard (2002)
	double R{ 8.3145 }; // gas constant
    double T_0{ 298.15 }; // reference temperature
    double P_0{ 1 }; // reference pressure [bar]

    // water and empty hydrate constants
    double gw_00{ -228700 }; // formation Gibbs energy of H2O in ideal gas
    double hw_00{ -242000 }; // enthalpy of formation of H2O in ideal gas
    std::vector<double> hw_a = { 3.8747*R, 0.0231E-2*R, 0.1269E-5*R, -0.4321E-9*R }; // ideal gas heat capacity parameters of H2O in ideal gas
    
    std::unordered_map<std::string, double> g_beta0 = {{"sI", -235537.85}, {"sII", -235627.53}, {"sH", -235491.02}}; // formation Gibbs energy of pure hydrate phases (sI, sII, sH)
    std::unordered_map<std::string, double> h_beta0 = {{"sI", -291758.77}, {"sII", -292044.10}, {"sH", -291979.26}}; // enthalpy of formation of pure hydrate phases (sI, sII, sH)
    std::vector<double> h_beta_a = { 0.735409713 * R, 1.4180551E-2 * R, -1.72746E-5 * R, 63.5104E-9 * R }; // heat capacity parameters of pure hydrate phases
    
    std::unordered_map<std::string, double> v_0{{"sI", 22.712}, {"sII", 22.9456}, {"sH", 24.2126}}; // molar volume of pure hydrate phases at reference p, T (sI, sII, sH)
    std::unordered_map<std::string, std::vector<double>> v_a = {
        {"sI", {3.384960E-4, 5.400990E-7, -4.769460E-11, 3E-5}},   // molar volume parameters of pure sI hydrate
        {"sII", {2.029776E-4, 1.851168E-7, -1.879455E-10, 3E-6}},
        {"sH", {3.575490E-4, 6.294390E-7, 0, 3E-7}}};

    // cage parameters
    double l0 = 1E-13;
    std::unordered_map<std::string, double> a = {{"sI", 25.74}, {"sII", 260.}, {"sH", 0.}}; // parameter a in eq. 4.39 (sI, sII, sH)
    std::unordered_map<std::string, double> b = {{"sI", -481.32}, {"sII", -68.64}, {"sH", 0.}}; // parameter b in eq. 4.39 (sI, sII, sH)
    std::unordered_map<std::string, double> nH2O = {{"sI", 46.}, {"sII", 136.}, {"sH", 34.}}; // total number of H2O molecules in hydrate structure (sI, sII, sH)
    std::unordered_map<std::string, double> a_0 = {{"sI", 11.99245}, {"sII", 17.1}, {"sH", 11.09826}}; // standard lattice parameters for sI, sII, sH [eq. 4.42, p. ]
    std::unordered_map<std::string, int> n_cages = {{"sI", 2}, {"sII", 2}, {"sH", 3}};

    std::unordered_map<std::string, std::vector<int>> zm = {{"sI", {20, 24}}, {"sII", {20, 28}}, {"sH", {20, 20, 36}}}; // #waters in cage sI
    std::unordered_map<std::string, std::vector<int>> zn = {{"sI", {8, 12, 8, 4, 8, 4}}, {"sII", {2, 6, 12, 12, 12, 4}}, {"sH", {20, 20, 36}}}; // #water in layers
    std::unordered_map<std::string, std::vector<int>> shells = {{"sI", {2, 4}}, {"sII", {3, 3}}, {"sH", {1, 1, 1}}};
    std::unordered_map<std::string, std::vector<double>> Nm = {{"sI", {2, 6}}, {"sII", {16, 8}}, {"sH", {3, 2, 1}}}; // number per unit cell
    std::unordered_map<std::string, std::vector<double>> Rn {
        {"sI", {3.83E-10, 3.96E-10, 4.06E-10, 4.25E-10, 4.47E-10, 4.645E-10}},  // radius of layers S0, S1, L0, L1, L2, L3 (sI)
        {"sII", {3.748E-10, 3.845E-10, 3.956E-10, 4.635E-10, 4.715E-10, 4.729E-10}},  // radius of layers S0, S1, S2, L0, L1, L2 (sII)
        {"sH", {3.91E-10, 4.06E-10, 5.71E-10}}};  // radius of layers S, M, L (sH)

    // guest parameters
    std::unordered_map<std::string, double> ai = {{"CO2", 0.6805E-10}, {"N2", 0.3526E-10},  {"H2S", 0.3600E-10}, {"C1", 0.3834E-10}, {"C2", 0.5651E-10}, {"C3", 0.6502E-10}, {"nC4", 0.9379E-10}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0}, }; // hard core radius
    std::unordered_map<std::string, double> sigma = {{"CO2", 2.97638E-10}, {"N2", 3.13512E-10}, {"H2S", 3.10000E-10}, {"C1", 3.14393E-10}, {"C2", 3.24693E-10}, {"C3", 3.41670E-10}, {"nC4", 3.51726E-10}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0}, }; // soft core radius
    std::unordered_map<std::string, double> eik = {{"CO2", 175.405}, {"N2", 127.426}, {"H2S", 212.047}, {"C1", 155.593}, {"C2", 188.181}, {"C3", 192.855}, {"nC4", 197.254}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0}, }; // potential well depth/k
    std::unordered_map<std::string, double> Di = {{"CO2", 4.603}, {"N2", 4.177}, {"H2S", 4.308}, {"C1", 4.247}, {"C2", 5.076}, {"C3", 5.745}, {"nC4", 6.336}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0}, }; // molecular diameter of component i
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> dr_im = {
        {"sI", {{"CO2", {0, 5.8282E-3}}, {"N2", {1.7377E-2, 0}}, {"H2S", {1.7921E-2, 0}}, {"C1", {1.7668E-2, 0}}, {"C2", {0, 1.5773E-2}}, {"C3", {0, 2.9839E-2}}, {"nC4", {0, 0}}, {"nC5", {0, 0}}, {"nC6", {0, 0}}, {"nC7", {0, 0}},}}, // repulsive const in cage S, L, sI
        {"sII", {{"CO2", {2.2758E-3, 1.2242E-2}}, {"N2", {2.0652E-3, 1.1295E-2}}, {"H2S", {2.1299E-3, 1.1350E-2}}, {"C1", {2.0998E-3, 1.1383E-2}}, {"C2", {2.5097E-3, 1.4973E-2}}, {"C3", {0,  2.5576E-2}}, {"nC4", {0, 3.6593E-2}}, {"nC5", {0, 0}}, {"nC6", {0, 0}}, {"nC7", {0, 0}},}}
    }; // repulsive const in cage S, L, sII
    
    // std::vector<double> dr_im3 = { 0 }; // volume of sH is assumed constant [p. 66]
    std::unordered_map<std::string, std::unordered_map<std::string, double>> kappa_iH = {  // compressibility parameters p. 118
        {"sI", {{"CO2", 1E-6}, {"N2", 1.1E-5}, {"H2S", 5E-6}, {"C1", 1E-5}, {"C2", 1E-8}, {"C3", 1E-7}, {"nC4", 0}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0},}},
        {"sII", {{"CO2", 1E-5}, {"N2", 1.1E-5}, {"H2S", 1E-5}, {"C1", 5E-5}, {"C2", 1E-7}, {"C3", 1E-6}, {"nC4", 1E-8}, {"nC5", 0}, {"nC6", 0}, {"nC7", 0},}}};
}

using namespace std;

void VdWP::parameters(double p_, double T_) {
    // Fugacity of water in hydrate phase following modified VdW-P (Ballard, 2002)
    // Initializes all composition-independent parameters:
    // - Ideal gas Gibbs energy of water [eq. 3.1-3.4]
    // - Gibbs energy of water in empty hydrate lattice [eq. 3.47]
    // Each time calculating fugacity of water in hydrate, evaluate:
    // - Contribution of cage occupancy to total energy of hydrate [eq. 3.44]
    // - Activity of water in the hydrate phase [eq. 4.38]
    // - Chemical potential of water in hydrate [eq. 4.35]
    p = p_; pp = p_; T = T_; TT = T_;
    water_index = std::distance(components.begin(), std::find(components.begin(), components.end(), "H2O"));

    zm = hyd_par::zm[phase];
    zn = hyd_par::zn[phase];
    n_shells = hyd_par::shells[phase];
    Rn = hyd_par::Rn[phase];
    Nm = hyd_par::Nm[phase];
    n_cages = hyd_par::n_cages[phase];
    nH2O = hyd_par::nH2O[phase];

    // Ideal gas Gibbs energy
    double g_wo = hyd_par::gw_00/(hyd_par::R * hyd_par::T_0); // Gibbs energy of formation of water in ideal gas at reference conditions [constant]
    double h_wo = hydrate_integrals(hyd_par::T_0, TT, 20, 0); // molar enthalpy of formation of water in ideal gas [eq. 3.3]
    g_w0 = g_wo - h_wo; // Gibbs energy of formation of water in ideal gas [eq. 3.2]

    // Gibbs energy of water in empty hydrate lattice
    double g_beta0 = hyd_par::g_beta0[phase] / (hyd_par::R * hyd_par::T_0); // Gibbs energy of formation of empty lattice at reference conditions [constant]
    double h_beta = hydrate_integrals(hyd_par::T_0, TT, 20, 1); // molar enthalpy of formation of empty lattice [eq. 4.40]
    double v_beta = hydrate_integrals(hyd_par::P_0, pp, 20, 2); // molar volume of standard hydrate [eq. 4.41]
    g_beta = g_beta0 - h_beta + v_beta; // Gibbs energy of formation of water in the standard state [eq. 3.47]
    
    return;
}

double VdWP::fwH(std::vector<double> f0) {
    // Calculate fugacity of water in the hydrate
    // - Contribution of cage occupancy to total energy of hydrate [eq. 3.44]
    // - Activity of water in the hydrate [eq. 4.38]
    // - Chemical potential of water in the hydrate [eq. 4.35]

    // Contribution of cage occupancy to total energy of hydrate
    double dmu_H = dmuH(f0);

    // Activity of water in the hydrate phase
    TT = hyd_par::T_0;
    double dv = dvH(hyd_par::P_0); // molar volume difference between standard hydrate and the real hydrate at reference conditions
    double dg_w0b = hyd_par::a[phase]*1E6*dv; // eq. 4.39, used in eq. 4.38
    double dh_w0b = hyd_par::b[phase]*1E6*dv; // eq. 4.39, used in eq. 4.38
    TT = T;
    double int_dvH = hydrate_integrals(hyd_par::P_0, pp, 20, 4); // integral of volume change over pressure in eq. 4.38
    double lnj = dg_w0b/(hyd_par::R*hyd_par::T_0) + dh_w0b/hyd_par::R * (1/TT - 1/hyd_par::T_0) + int_dvH; // activity coefficient of water in hydrate [eq. 4.38]

    // Fugacity of water in hydrate
    double mu_wH = g_beta + dmu_H + lnj; // Chemical potential of water in hydrate [eq. 4.35]
    double f_w0 = 1; // fugacity of ideal gas at reference pressure
    double f_wH = f_w0*exp(mu_wH - g_w0); // eq. 4.47

    return f_wH;
}

double VdWP::dmuH(std::vector<double> f0) {
    // Calculate Langmuir constant of each guest i in each cage m
    std::vector<double> C_im(n_cages*NC); // contains Langmuir constants for each component i in each cage m
    for (int i = 0; i < NC; i++) 
    {
        comp = components[i];
        if (comp != "H2O") 
        {
            R1_index = 0;
            for (int m = 0; m < n_cages; m++) 
            {
                cage_index = m;
                double Cim_int = hydrate_integrals(hyd_par::l0, Rn[R1_index]-hyd_par::ai[comp], 20, 3); // integral in [m^3]
                C_im[NC*m + i] = 4 * M_PI / (M_kB * TT) * Cim_int * 1E5;
                R1_index += n_shells[m];
            }
        }
    }
    
    // Fractional occupancy of cage m by component i
    theta_im = std::vector<double>(n_cages*NC);
    for (int m = 0; m < n_cages; m++) 
    {
        double sum_cf{ 0 };
        for (int j = 0; j < NC; j++) 
        {
            sum_cf += C_im[NC*m + j] * f0[j];
        }
        for (int i = 0; i < NC; i++)
        {
            theta_im[NC*m + i] = C_im[NC*m + i] * f0[i] / (1 + sum_cf);
        }
    }

    // Contribution of cage occupancy to chemical potential
    double dmu_H{ 0 };
    for (int m = 0; m < n_cages; m++) 
    {
        double sumtheta{ 0 };
        for (int i = 0; i < NC; i++) 
        {
            sumtheta += theta_im[NC*m + i];
        }
        dmu_H += Nm[m] / nH2O * log(1-sumtheta);
    }
    return dmu_H;
}

double VdWP::dvH(double x) {
    // cage occupancy volume change

    // Function f(theta_im) [eq. 6.18/6.19]
    // fractional occupancy average molecular diameter
    double sum_theta_D{ 0 }, sum_theta{ 0 };
    for (int i = 0; i < NC; i++) 
    {
        for (int m = 0; m < n_cages; m++) 
        {
            sum_theta_D += theta_im[NC*m + i] * hyd_par::Di[components[i]];
            sum_theta += theta_im[NC*m + i];
        }
    }
    double D_ = sum_theta_D / sum_theta;

    std::vector<double> f_im(NC*n_cages);
    for (int i = 0; i < NC; i++) 
    {
        if (components[i] != "H2O") 
        {
            // For first cage (type 5^12), use
            f_im[0 + i] = (1 + zm[0] / nH2O) * theta_im[0 + i] / (1 + zm[0] / nH2O * theta_im[0 + i]) * exp(hyd_par::Di[components[i]] - D_);
            // For other cages, use
            for (int m = 1; m < n_cages; m++) 
            {
                f_im[NC*m + i] = (1 + zm[m] / nH2O) * theta_im[NC*m + i] / (1 + zm[m] / nH2O * theta_im[NC*m + i]);
            }
        }
    }

    // Kappa & v0
    double v0;
    double kappa{ 0. };

    if (phase == "sH") 
    { 
        kappa = 1E-7; v0 = hyd_par::v_0["sH"]; 
    }
    else 
    {
        // V_0(x) [eq. 4.42]
        double sum_m{ 0 };
        for (int m = 0; m < n_cages; m++) 
        {
            double sum_i{ 0 };
            for (int i = 0; i < NC; i++) 
            {
                if (components[i] != "H2O") 
                {
                    sum_i += f_im[NC*m + i] * hyd_par::dr_im[phase][components[i]][m];
                }
            }
            sum_m += Nm[m] * sum_i;
        }

        for (int i = 0; i < NC; i++) 
        {
            if (components[i] != "H2O") 
            {
                kappa += hyd_par::kappa_iH[phase][components[i]] * theta_im[NC + i];
            }
        }
        v0 = pow((hyd_par::a_0[phase] + sum_m), 3) * M_NA / nH2O * 1E-24; // compositional dependence of molar volume [4.42]
    }

    // Molar volume of hydrate [eq. 4.41]
    double vH = v0 * exp(hyd_par::v_a[phase][0]*(TT-hyd_par::T_0) + hyd_par::v_a[phase][1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[phase][2]*pow((TT-hyd_par::T_0), 3) - kappa*(x-hyd_par::P_0)) * 1E-6; // eq. 4.41

    // Volume of empty hydrate lattice [eq. 4.41]
    double vB = hyd_par::v_0[phase] * exp(hyd_par::v_a[phase][0]*(TT-hyd_par::T_0) + hyd_par::v_a[phase][1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[phase][2]*pow((TT-hyd_par::T_0), 3) - hyd_par::v_a[phase][3]*(x-hyd_par::P_0)) * 1E-6;
    return vH - vB;
}

double VdWP::hydrate_functions(double x, int function) {
    if (function == 0)  // hio
    {
        // ideal gas enthalpy [eq. 3.3]
        double hi0 = hyd_par::hw_00 + hyd_par::hw_a[0]*(x-hyd_par::T_0) + 1/2*hyd_par::hw_a[1]*(pow(x, 2)-pow(hyd_par::T_0, 2)) + 1/3*hyd_par::hw_a[2]*(pow(x, 3)-pow(hyd_par::T_0, 3)) + 1/4*hyd_par::hw_a[3]*(pow(x, 4)-pow(hyd_par::T_0, 4));
        return hi0 / (hyd_par::R * pow(x, 2));
    }
    else if (function == 1)  // h_beta
    {
        // molar enthalpy of empty hydrate lattice [eq. 4.40]
        double h_beta_ = hyd_par::h_beta0[phase] + hyd_par::h_beta_a[0]*(x-hyd_par::T_0) + 1/2.0 * hyd_par::h_beta_a[1]*(pow(x, 2)-pow(hyd_par::T_0, 2)) + 1/3.0 * hyd_par::h_beta_a[2]*(pow(x, 3)-pow(hyd_par::T_0, 3)) + 1/4.0 * hyd_par::h_beta_a[3]*(pow(x, 4)-pow(hyd_par::T_0, 4));
        return h_beta_ / (hyd_par::R * pow(x, 2));
    }
    else  // v_beta
    {
        // molar volume of empty hydrate lattice [eq. 4.41]
        double v_beta_ = hyd_par::v_0[phase] * exp(hyd_par::v_a[phase][0]*(TT-hyd_par::T_0) + hyd_par::v_a[phase][1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[phase][2]*pow((TT-hyd_par::T_0), 3) - hyd_par::v_a[phase][3]*(x-hyd_par::P_0)) * 1E-6;
        return v_beta_ / (hyd_par::R * 1E-5 * TT);
    }
}

double VdWP::pot(double x) {
    // hydrate cage cell potential omega(r) [eq. 3.43a]
    double om{ 0 };
    // hard core, soft core radius and potential well depth of guest molecule
    double ai = hyd_par::ai[comp]; double sigmai = hyd_par::sigma[comp]; double eik = hyd_par::eik[comp];

    // Loop over shell layers of cage m [eq. 4.45]
    for (int n = 0; n < n_shells[cage_index]; n++) 
    {
        double Rn_ = Rn[R1_index + n];
        int zn_ = zn[R1_index + n];
        double delta10 = 1. / 10.0 * (pow((1. - x / Rn_ - ai / Rn_), -10.) - pow((1. + x / Rn_ - ai / Rn_), -10.));
        double delta11 = 1. / 11.0 * (pow((1. - x / Rn_ - ai / Rn_), -11.) - pow((1. + x / Rn_ - ai / Rn_), -11.));
        double delta4 = 1. / 4.0 * (pow((1. - x / Rn_ - ai / Rn_), -4.) - pow((1. + x / Rn_ - ai / Rn_), -4.));
        double delta5 = 1. / 5.0 * (pow((1. - x / Rn_ - ai / Rn_), -5.) - pow((1. + x / Rn_ - ai / Rn_), -5.));
        om += 2.0 * zn_ * (pow(sigmai, 12.) / (pow(Rn_, 11.) * x) * (delta10 + ai / Rn_ * delta11) -
                         pow(sigmai, 6.) / (pow(Rn_, 5.) * x) * (delta4 + ai / Rn_ * delta5));
        
    }
    // term in integral for Langmuir constant C_im [eq. 3.42]
    double pot = exp(-eik / TT * om) * pow(x, 2);

    return pot;
}

double VdWP::hydrate_integrals(double a, double b, int steps, int function) {
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
            s += h*((hydrate_functions(x, 0) + 4*hydrate_functions(x+h/2, 0) + hydrate_functions(x+h, 0))/6);
        };
        break;
    case 1: // molar enthalpy of empty hydrate lattice
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = T
            s += h*((hydrate_functions(x, 1) + 4*hydrate_functions(x+h/2, 1) + hydrate_functions(x+h, 1))/6);
        };
        break;
    case 2: // molar volume of empty hydrate lattice
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = p
            s += h*((hydrate_functions(x, 2) + 4*hydrate_functions(x+h/2, 2) + hydrate_functions(x+h, 2))/6);
        };
        break;
    case 3: // hydrate cage cell potential
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = r
            s += h*((pot(x) + 4*pot(x+h/2) + pot(x+h))/6);
            if (pot(x+h) < 1E-200) {break;} // otherwise integral blows up
        };
        break;
    case 4: // cage occupancy volume change
        for (int i = 0; i < steps; i++) 
        {
            x = a + h*i; // x = p
            s += h*((dvH(x) + 4*dvH(x+h/2) + dvH(x+h))/6) / (hyd_par::R*1E-5*TT);
        }
        break;
    }
    return s;
}

std::vector<double> VdWP::xH() {
    std::vector<double> x_H(NC);

    double denominator{ 1. }; // denominator of eq. 3.50
    for (int m = 0; m < n_cages; m++) 
    {
        for (int i = 0; i < NC; i++) 
        {
            denominator += Nm[m]/nH2O * theta_im[NC*m + i];
        }
    }
    
    x_H[water_index] = 1.;
    for (int i = 0; i < NC; i++) 
    {
        if (components[i] != "H2O") 
        {
            double numerator{ 0 }; // numerator of eq. 3.50
            for (int m = 0; m < n_cages; m++) 
            { 
                numerator += Nm[m]/nH2O * theta_im[NC*m + i]; 
            }
            x_H[i] = numerator / denominator;
            x_H[water_index] -= numerator / denominator;
        }
    }
    return x_H;
}

std::vector<double> VdWP::fugacityCoefficient(std::vector<double> f) {
    f[water_index] = fwH(f);
    std::vector<double> x_H = xH(); // calculate phase composition

    std::vector<double> phi(NC);
    for (int i = 0; i < NC; i++) 
    {
        phi[i] = f[i]/(x_H[i] * p);
    }
    return phi;
}