#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <numeric>

#include "eos.h"
#include "../global/misc.h"

#define M_PI 3.14159265358979323846 /* pi */

namespace hyd_par {
    // Following hydrate fugacity model described in Ballard (2002)
	double R{ 8.3145 }; // gas constant
    double T_0{ 298.15 }; // reference temperature
    double P_0{ 1 }; // reference pressure [bar]
    double N_A{ 6.022E23 }; // Avogadro number
    double k_B{ 1.380649E-23 }; // Boltzmann constant

    // water and empty hydrate constants
    double gw_00{ -228700 }; // formation Gibbs energy of H2O in ideal gas
    double hw_00{ -242000 }; // enthalpy of formation of H2O in ideal gas
    std::vector<double> hw_a = { 3.8747*R, 0.0231E-2*R, 0.1269E-5*R, -0.4321E-9*R }; // ideal gas heat capacity parameters of H2O in ideal gas
    
    std::vector<double> g_beta0 = { -235537.85, -235627.53, -235491.02 }; // formation Gibbs energy of pure hydrate phases (sI, sII, sH)
    std::vector<double> h_beta0{ -291758.77, -292044.10, -291979.26 }; // enthalpy of formation of pure hydrate phases (sI, sII, sH)
    std::vector<double> h_beta_a = { 0.735409713 * R, 1.4180551E-2 * R, -1.72746E-5 * R, 63.5104E-9 * R }; // heat capacity parameters of pure hydrate phases
    
    std::vector<double> v_0{ 22.712, 22.9456, 24.2126 }; // molar volume of pure hydrate phases at reference p, T (sI, sII, sH)
    std::vector<double> v_a = { 3.384960E-4, 5.400990E-7, -4.769460E-11, 3E-5,   // molar volume parameters of pure sI hydrate
                                2.029776E-4, 1.851168E-7, -1.879455E-10, 3E-6,   // sII
                                3.575490E-4, 6.294390E-7, 0,             3E-7 }; // sH

    // cage parameters
    double l0 = 1E-13;
    std::vector<double> l = { 3.83E-10, 3.748E-10, 3.91E-10 }; // 
    std::vector<double> a = { 25.74, 260.0, 0. }; // parameter a in eq. 4.39 (sI, sII, sH)
    std::vector<double> b = { -481.32, -68.64, 0. }; // parameter b in eq. 4.39 (sI, sII, sH)
    std::vector<double> nH2O = { 46.0, 136.0, 34.0 }; // total number of H2O molecules in hydrate structure (sI, sII, sH)
    std::vector<double> a_0 = { 11.99245, 17.1, 11.09826 }; // standard lattice parameters for sI, sII, sH [eq. 4.42, p. ]

    std::vector<double> Nm1 = { 2, 6 }; // number per unit cell (sI)
    std::vector<double> Nm2 = { 16, 8 }; // number per unit cell (sII)
    std::vector<double> Nm3 = { 3, 2, 1 }; // number per unit cell (sH)
    std::vector<int> zm1 = { 20, 24 }; // #waters in cage sI
    std::vector<int> zm2 = { 20, 28 }; // #waters in cage sII
    std::vector<int> zm3 = { 20, 20, 36 }; // #waters in cage sH
    
    std::vector<int> zn3 = { 20, 20, 36 };  // #water in layers S, M, L (sH)
    
    std::vector<double> Rn3 = { 3.91E-10, 4.06E-10, 5.71E-10 };  // radius of layers S, M, L (sH)
    
    std::vector<int> shells3 = { 1, 1, 1 }; // number of layers in each cage sH

    std::vector<int> zn1 = { 8, 12, 8, 4, 8, 4 };  // #waters in layers S0, S1, L0, L1, L2, L3 (sI)
    std::vector<int> zn2 = { 2, 6, 12, 12, 12, 4 };  // #water in layers S0, S1, S2, L0, L1, L2 (sII)
    std::vector<double> Rn1 = { 3.83E-10, 3.96E-10, 4.06E-10, 4.25E-10, 4.47E-10, 4.645E-10 };  // radius of layers S0, S1, L0, L1, L2, L3 (sI)
    std::vector<double> Rn2 = { 3.748E-10, 3.845E-10, 3.956E-10, 4.635E-10, 4.715E-10, 4.729E-10 };  // radius of layers S0, S1, S2, L0, L1, L2 (sII)
    std::vector<int> shells1 = { 2, 4 }; // number of layers in each cage sI
    std::vector<int> shells2 = { 3, 3 }; // number of layers in each cage sII

    // std::vector<double> Rn1 = { 3.908E-10, 4.326E-10 };  // radius of layers S, L (sI)
    // std::vector<double> Rn2 = { 3.902E-10, 4.683E-10 };  // radius of layers S, L (sII)
    // std::vector<int> zn1 = { 20, 24 };  // #waters in layers S, L (sI)
    // std::vector<int> zn2 = { 20, 28 };  // #water in layers S, L (sII)
    // std::vector<int> shells1 = { 1, 1 }; // number of layers in each cage sH
    // std::vector<int> shells2 = { 1, 1 }; // number of layers in each cage sH

    // guest parameters
    std::vector<std::string> comp{ "H2O", "CO2",   "N2",        "H2S",       "C1",        "C2",        "C3",        "nC4",       "nC5", "nC6", "nC7" };
    std::vector<double> ai = {     0, 0.6805E-10,  0.3526E-10,  0.3600E-10,  0.3834E-10,  0.5651E-10,  0.6502E-10,  0.9379E-10,  0, 0, 0 }; // hard core radius
    std::vector<double> sigma = {  0, 2.97638E-10, 3.13512E-10, 3.10000E-10, 3.14393E-10, 3.24693E-10, 3.41670E-10, 3.51726E-10, 0, 0, 0 }; // soft core radius
    std::vector<double> eik = {    0, 175.405,     127.426,     212.047,     155.593,     188.181,     192.855,     197.254,     0, 0, 0 }; // potential well depth/k
    std::vector<double> Di = {     0, 4.603,       4.177,       4.308,       4.247,       5.076,       5.745,       6.336,       0, 0, 0 }; // molecular diameter of component i
    std::vector<double> dr_im1 = { 0, 0,           1.7377E-2,   1.7921E-2,   1.7668E-2,   0,           0,           0,           0, 0, 0,   // repulsive const in cage S, sI
                                   0, 5.8282E-3,   0,           0,           0,           1.5773E-2,   2.9839E-2,   0,           0, 0, 0 }; // repulsive const in cage L, sI
    std::vector<double> dr_im2 = { 0, 2.2758E-3,   2.0652E-3,   2.1299E-3,   2.0998E-3,   2.5097E-3,   0,           0,           0, 0, 0,   // repulsive const in cage S, sII
                                   0, 1.2242E-2,   1.1295E-2,   1.1350E-2,   1.1383E-2,   1.4973E-2,   2.5576E-2,   3.6593E-2,   0, 0, 0 }; // repulsive const in cage L, sII
    std::vector<double> dr_im3 = { 0 }; // volume of sH is assumed constant [p. 66]
    std::vector<double> kappa_iH = { 0, 1E-6,      1.1E-5,      5E-6,        1E-5,        1E-8,        1E-7,        0,           0, 0, 0,
                                     0, 1E-5,      1.1E-5,      1E-5,        5E-5,        1E-7,        1E-6,        1E-8,        0, 0, 0 }; // compressibility parameters p. 118
    int nc = comp.size();
}

// namespace langmuir {
//     // Out of range langmuir constants according to C_im = A*exp(B/T)
//     std::vector<std::string> comp{ "H2O", "CO2",   "N2",        "H2S",       "C1",        "C2",        "C3",        "nC4",       "nC5", "nC6", "nC7" };
//     std::vector<double> A_im1 = { 0, 7.7765E-7,    3.9496E-5,   2.3444E-5,   8.3453E-5,   0,           0,           0, 0, 0, 0,
//                                   0, 520.5579E-7,  25.6897E-5,  7.2080E-5,   116.6313E-5, 3.5164E-6,   5.5707E-8,   0, 0, 0, 0 };
//     std::vector<double> A_im2 = { 0, 7.9970E-7,    4.8836E-5,   3.6415E-5,   5.4792E-5,   0,           0,           0, 0, 0, 0,
//                                   0, 6907.0012E-7, 201.3238E-5, 758.3575E-5, 829.8039E-5, 727.2717E-6, 597.9850E-8, 0, 0, 0, 0 };
//     std::vector<double> B_im1 = { 0, 2976.629,     2869.400,    4463.910,    2901.747,    0,           0,           0, 0, 0, 0,
//                                   0, 4674.690,     2680.372,    4073.045,    2959.901,    4226.997,    3537.025,    0, 0, 0, 0 };
//     std::vector<double> B_im2 = { 0, 2277.757,     2679.423,    3073.324,    2546.660,    0,           0,           0, 0, 0, 0,
//                                   0, 3370.363,     2226.480,    2495.937,    2629.194,    4440.484,    7118.782,    0, 0, 0, 0 };
//     int nc = comp.size();
// }

using namespace std;

Hydrate::Hydrate(std::vector<int> ci_, int phase_index) {
    // Fugacity of water in hydrate phase following modified VdW-P (Ballard, 2002)
    // Initializes all composition-independent parameters:
    // - Ideal gas Gibbs energy of water [eq. 3.1-3.4]
    // - Gibbs energy of water in empty hydrate lattice [eq. 3.47]
    // Each time calculating fugacity of water in hydrate, evaluate:
    // - Contribution of cage occupancy to total energy of hydrate [eq. 3.44]
    // - Activity of water in the hydrate phase [eq. 4.38]
    // - Chemical potential of water in hydrate [eq. 4.35]

    ci = ci_; NC = ci.size();
    pi = phase_index-3;
    f = std::vector<double>(NC);
    x_H = std::vector<double>(NC);
    water_index = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 0));

    // Define type-specific hydrate parameters
    if (pi == 0) { // sI
        zm = hyd_par::zm1; zn = hyd_par::zn1; n_shells = hyd_par::shells1;
        Rn = hyd_par::Rn1; dr_im = hyd_par::dr_im1; Nm = hyd_par::Nm1;
    } else if (pi == 1) { // sII
        zm = hyd_par::zm2; zn = hyd_par::zn2; n_shells = hyd_par::shells2;
        Rn = hyd_par::Rn2; dr_im = hyd_par::dr_im2; Nm = hyd_par::Nm2;
    } else { // sH
        zm = hyd_par::zm3; zn = hyd_par::zn3; n_shells = hyd_par::shells3;
        Rn = hyd_par::Rn3; dr_im = hyd_par::dr_im3; Nm = hyd_par::Nm3;
    }
    n_cages = zm.size();
    nH2O = hyd_par::nH2O[pi];

}

void Hydrate::init(double p, double T) {
    pp = p; TT = T;

    // Ideal gas Gibbs energy
    double g_wo = hyd_par::gw_00/(hyd_par::R * hyd_par::T_0); // Gibbs energy of formation of water in ideal gas at reference conditions [constant]
    double h_wo = Hyd_integrals(hyd_par::T_0, TT, 20, 0); // molar enthalpy of formation of water in ideal gas [eq. 3.3]
    g_w0 = g_wo - h_wo; // Gibbs energy of formation of water in ideal gas [eq. 3.2]

    // Gibbs energy of water in empty hydrate lattice
    double g_beta0 = hyd_par::g_beta0[pi] / (hyd_par::R * hyd_par::T_0); // Gibbs energy of formation of empty lattice at reference conditions [constant]
    double h_beta = Hyd_integrals(hyd_par::T_0, TT, 20, 1); // molar enthalpy of formation of empty lattice [eq. 4.40]
    double v_beta = Hyd_integrals(hyd_par::P_0, pp, 20, 2); // molar volume of standard hydrate [eq. 4.41]
    g_beta = g_beta0 - h_beta + v_beta; // Gibbs energy of formation of water in the standard state [eq. 3.47]
    
    return;
}

void Hydrate::hydrateFugacity(double T, const std::vector<double> f1) {
    // Calculate fugacity of water in the hydrate
    // - Contribution of cage occupancy to total energy of hydrate [eq. 3.44]
    // - Activity of water in the hydrate [eq. 4.38]
    // - Chemical potential of water in the hydrate [eq. 4.35]

    // Contribution of cage occupancy to total energy of hydrate
    double dmu_H = dmuH(f1);

    // Activity of water in the hydrate phase
    TT = hyd_par::T_0;
    double dv = dvH(hyd_par::P_0); // molar volume difference between standard hydrate and the real hydrate at reference conditions
    double dg_w0b = hyd_par::a[pi]*1E6*dv; // eq. 4.39, used in eq. 4.38
    double dh_w0b = hyd_par::b[pi]*1E6*dv; // eq. 4.39, used in eq. 4.38
    TT = T;
    double int_dvH = Hyd_integrals(hyd_par::P_0, pp, 20, 4); // integral of volume change over pressure in eq. 4.38
    double lnj = dg_w0b/(hyd_par::R*hyd_par::T_0) + dh_w0b/hyd_par::R*(1/TT - 1/hyd_par::T_0) + int_dvH; // activity coefficient of water in hydrate [eq. 4.38]

    // Fugacity of water in hydrate
    double mu_wH = g_beta + dmu_H + lnj; // Chemical potential of water in hydrate [eq. 4.35]
    double f_w0 = 1; // fugacity of ideal gas at reference pressure
    f_wH = f_w0*exp(mu_wH - g_w0); // eq. 4.47
    
    xH(); // calculate phase composition
}

double Hydrate::dmuH(std::vector<double> f1) {
    // Calculate Langmuir constant of each guest i in each cage m
    std::vector<double> C_im(n_cages*NC); // contains Langmuir constants for each component i in each cage m
    for (int i = 0; i < NC; i++) {
        comp_index = ci[i];
        if (comp_index != 0) {
            R1_index = 0;
            for (int m = 0; m < n_cages; m++) {
                cage_index = m;
                double Cim_int = Hyd_integrals(hyd_par::l0, Rn[R1_index]-hyd_par::ai[comp_index], 20, 3); // integral in [m^3]
                C_im[NC*m + i] = 4 * M_PI / (hyd_par::k_B * TT) * Cim_int * 1E5;
                R1_index += n_shells[m];
            }
        }
    }
    
    // Fractional occupancy of cage m by component i
    theta_im = std::vector<double>(n_cages*NC);
    for (int m = 0; m < n_cages; m++) {
        double sum_cf{ 0 };
        for (int j = 0; j < NC; j++) {
            sum_cf += C_im[NC*m + j] * f1[j];
        }
        for (int i = 0; i < NC; i++) {
            theta_im[NC*m + i] = C_im[NC*m + i] * f1[i] / (1 + sum_cf);
        }
    }

    // Contribution of cage occupancy to chemical potential
    double dmu_H{ 0 };
    for (int m = 0; m < n_cages; m++) {
        double sumtheta{ 0 };
        for (int i = 0; i < NC; i++) {
            sumtheta += theta_im[NC*m + i];
        }
        dmu_H += Nm[m] / nH2O * log(1-sumtheta);
    }
    return dmu_H;
}

double Hydrate::dvH(double x) {
    // cage occupancy volume change

    // Function f(theta_im) [eq. 6.18/6.19]
    // fractional occupancy average molecular diameter
    double sum_theta_D{ 0 }, sum_theta{ 0 };
    for (int i = 0; i < NC; i++) {
        for (int m = 0; m < n_cages; m++) {
            sum_theta_D += theta_im[NC*m + i] * hyd_par::Di[ci[i]];
            sum_theta += theta_im[NC*m + i];
        }
    }
    double D_ = sum_theta_D / sum_theta;

    std::vector<double> f_im(NC*n_cages);
    for (int i = 0; i < NC; i++) {
        if (ci[i] != 0) {
            // For first cage (type 5^12), use
            f_im[0 + i] = (1 + zm[0] / nH2O) * theta_im[0 + i] / (1 + zm[0] / nH2O * theta_im[0 + i]) * exp(hyd_par::Di[ci[i]] - D_);
            // For other cages, use
            for (int m = 1; m < n_cages; m++) {
                f_im[NC*m + i] = (1 + zm[m] / nH2O) * theta_im[NC*m + i] / (1 + zm[m] / nH2O * theta_im[NC*m + i]);
            }
        }
    }

    // V_0(x) [eq. 4.42]
    double sum_m{ 0 };
    for (int m = 0; m < n_cages; m++) {
        double sum_i{ 0 };
        for (int i = 0; i < NC; i++) {
            sum_i += f_im[NC*m + i] * dr_im[hyd_par::nc*m + ci[i]];
        }
        sum_m += Nm[m] * sum_i;
    }

    // Kappa
    double kappa{ 0. };
    if (pi == 2) { kappa = 1E-7; }
    else {
        for (int i = 0; i < NC; i++) {
            kappa += hyd_par::kappa_iH[hyd_par::nc*pi + ci[i]] * theta_im[NC + i];
        }
    }
    kappa = hyd_par::v_a[4*pi + 3];


    // Molar volume of hydrate [eq. 4.41]
    double v0 = pow((hyd_par::a_0[pi] + sum_m), 3) * hyd_par::N_A / nH2O * 1E-24; // compositional dependence of molar volume [4.42]
    double vH = v0 * exp(hyd_par::v_a[4*pi + 0]*(TT-hyd_par::T_0) + hyd_par::v_a[4*pi + 1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[4*pi + 2]*pow((TT-hyd_par::T_0), 3) - kappa*(x-hyd_par::P_0)) * 1E-6; // eq. 4.41

    // Volume of empty hydrate lattice [eq. 4.41]
    double vB = hyd_par::v_0[pi] * exp(hyd_par::v_a[4*pi + 0]*(TT-hyd_par::T_0) + hyd_par::v_a[4*pi + 1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[4*pi + 2]*pow((TT-hyd_par::T_0), 3) - hyd_par::v_a[4*pi + 3]*(x-hyd_par::P_0)) * 1E-6;
    return vH - vB;
}

double Hydrate::h_io(double x) {
    // ideal gas enthalpy [eq. 3.3]
    double hi0 = hyd_par::hw_00 + hyd_par::hw_a[0]*(x-hyd_par::T_0) + 1/2*hyd_par::hw_a[1]*(pow(x, 2)-pow(hyd_par::T_0, 2)) + 1/3*hyd_par::hw_a[2]*(pow(x, 3)-pow(hyd_par::T_0, 3)) + 1/4*hyd_par::hw_a[3]*(pow(x, 4)-pow(hyd_par::T_0, 4));
    return hi0 / (hyd_par::R * pow(x, 2));
}

double Hydrate::h_beta(double x) {
    // molar enthalpy of empty hydrate lattice [eq. 4.40]
    double h_beta_ = hyd_par::h_beta0[pi] + hyd_par::h_beta_a[0]*(x-hyd_par::T_0) + 1/2.0 * hyd_par::h_beta_a[1]*(pow(x, 2)-pow(hyd_par::T_0, 2)) + 1/3.0 * hyd_par::h_beta_a[2]*(pow(x, 3)-pow(hyd_par::T_0, 3)) + 1/4.0 * hyd_par::h_beta_a[3]*(pow(x, 4)-pow(hyd_par::T_0, 4));
    return h_beta_ / (hyd_par::R * pow(x, 2));
}

double Hydrate::v_beta(double x) {
    // molar volume of empty hydrate lattice [eq. 4.41]
    double v_beta_ = hyd_par::v_0[pi] * exp(hyd_par::v_a[4*pi + 0]*(TT-hyd_par::T_0) + hyd_par::v_a[4*pi + 1]*pow((TT-hyd_par::T_0), 2) + hyd_par::v_a[4*pi + 2]*pow((TT-hyd_par::T_0), 3) - hyd_par::v_a[4*pi + 3]*(x-hyd_par::P_0)) * 1E-6;
    return v_beta_ / (hyd_par::R * 1E-5 * TT);
}

double Hydrate::pot(double x) {
    // hydrate cage cell potential omega(r) [eq. 3.43a]
    double om{ 0 };
    // hard core, soft core radius and potential well depth of guest molecule
    double ai = hyd_par::ai[comp_index]; double sigmai = hyd_par::sigma[comp_index]; double eik = hyd_par::eik[comp_index];
    // cout << comp_index << " " << ai << " " << sigmai << " " << eik << endl;

    // Loop over shell layers of cage m [eq. 4.45]
    for (int n = 0; n < n_shells[cage_index]; n++) {
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
    // cout << "Om/kT " << x << " " << comp_index << " " << cage_index << " " << eik/TT*om << endl;
    // cout << "Pot " << x << " " << comp_index << " " << cage_index << " " << pot << endl;

    return pot;
}

double Hydrate::Hyd_integrals(double a, double b, int steps, int function) {
    // integrals solved numerically with simpson's rule
    double s = 0;
    double h = (b-a)/steps;
    double x;

    switch (function)
    {
    case 0: // ideal gas enthalpy
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = T
            s += h*((h_io(x) + 4*h_io(x+h/2) + h_io(x+h))/6);
        };
        break;
    case 1: // molar enthalpy of empty hydrate lattice
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = T
            s += h*((h_beta(x) + 4*h_beta(x+h/2) + h_beta(x+h))/6);
        };
        break;
    case 2: // molar volume of empty hydrate lattice
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = p
            s += h*((v_beta(x) + 4*v_beta(x+h/2) + v_beta(x+h))/6);
        };
        break;
    case 3: // hydrate cage cell potential
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = r
            s += h*((pot(x) + 4*pot(x+h/2) + pot(x+h))/6);
            if (pot(x+h) < 1E-200) {break;} // otherwise integral blows up
        };
        break;
    case 4: // cage occupancy volume change
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = p
            s += h*((dvH(x) + 4*dvH(x+h/2) + dvH(x+h))/6) / (hyd_par::R*1E-5*TT);
        }
        break;
    }
    return s;
}

void Hydrate::xH() {
    double denominator{ 1. }; // denominator of eq. 3.50
    for (int m = 0; m < n_cages; m++) {
        for (int i = 0; i < NC; i++) {
            denominator += Nm[m]/nH2O * theta_im[NC*m + i];
        }
    }
    
    x_H[water_index] = 1.;
    for (int i = 0; i < NC; i++) {
        if (ci[i] != 0) {
            double numerator{ 0 }; // numerator of eq. 3.50
            for (int m = 0; m < n_cages; m++) { numerator += Nm[m]/nH2O * theta_im[NC*m + i]; }
            x_H[i] = numerator / denominator;
            x_H[water_index] -= numerator / denominator;
        }
    }
    return;
}

std::vector<double> Hydrate::fugacityCoefficient(double p, double T, std::vector<double> x) {
    f = x; // fugacity 
    hydrateFugacity(T, f);
    std::vector<double> phi_ij(NC);
    f[water_index] = f_wH;
    for (int i = 0; i < NC; i++) {
        phi_ij[i] = f[i]/(x_H[i] * p);
    }
    return phi_ij;
}