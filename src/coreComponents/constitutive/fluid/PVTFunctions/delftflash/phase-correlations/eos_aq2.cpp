#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "eos.h"
#include "../global/misc.h"

#define M_PI 3.14159265358979323846 /* pi */

using namespace std;

namespace aq2_par {
    double R = 8.3145; double T_0 = 298.15; double P_0 = 1.0;
    //components{                "H2O",        "CO2",        "N2",         "H2S",        "C1",         "C2",         "C3",         "nC4",        "nC5",        "nC6",       "nC7" };
    std::vector<double> gi_0 =  {-237129,      -385974,       18188,       -27920,       -34451,       -17000,       -7550,        -940,          9160,         18200,       27500 };     // gibbs energy of pure H2O or 1 molal solution at p0, T0
    std::vector<double> hi_0 =  {-285830,      -413798,      -10439,       -37660,       -87906,       -103136,      -131000,      -152000,      -173887,      -199200,     -225000 };    // molar enthalpy of pure H2O or 1 molal solution at p0, T0
    std::vector<double> gi_00 = {-228700,      -394600,       0,           -33100,       -50830,       -32900,       -23500,       -17200,       -8370,        -290,         8120 };      // gibbs energy of ideal gas at p0, T0
    std::vector<double> hi_00 = {-242000,      -393800,       0,           -20200,       -74900,       -84720,       -103900,      -126200,      -146500,      -167300,     -187900 };    // molar enthalpy of ideal gas at p0, T0
    std::vector<double> hi_a =  { 3.8747*R,     2.6751*R,     3.4736*R,     3.5577*R,     2.3902*R,     0.8293*R,    -0.4861*R,     0.4755*R,     0.8142*R,     0.8338*R,   -0.6184*R,    // a0
                                  0.0231E-2*R,  0.7188E-2*R, -0.0189E-2*R,  0.1574E-2*R,  0.6039E-2*R,  2.0752E-2*R,  3.6629E-2*R,  4.4650E-2*R,  5.4598E-2*R,  6.6373E-2*R, 8.1268E-2*R, // a1
                                  0.1269E-5*R, -0.4208E-5*R,  0.0971E-5*R,  0.0686E-5*R,  0.1525E-5*R, -0.7699E-5*R, -1.8895E-5*R, -2.2041E-5*R, -2.6997E-5*R, -3.444E-5*R, -4.388E-5*R,  // a2
                                 -0.4321E-9*R,  0.8977E-9*R, -0.3453E-9*R, -0.3959E-9*R, -1.3234E-9*R,  0.8756E-9*R,  3.8143E-9*R,  4.2068E-9*R,  5.0824E-9*R,  6.9342E-9*R, 9.2037E-9*R };  // a3
    std::vector<double> omega = { 0,            -8368,        -145101,      -41840,       -133009,      -169870,      -211418,      -253592,      -300955,      -335180,     -380158 };
    std::vector<double> cp =    { 0,            167.50,       149.75,       135.14,       176.12,       226.67,       277.52,       330.77,       373.24,       424.53,      472.37,      // partial molar heat capacity c1
                                  0,            5304066,      5046230,      2850801,      6310762,      9011737,      11749531,     14610096,     16955051,     19680558,    22283347 };  // partial molar heat capacity c2
    std::vector<double> vp =    { 0,            2.614,        2.596,        2.724,        2.829,        3.612,        4.503,        5.500,        6.282,        7.175,       8.064,       // partial molar volume v1
                                  0,            3125.9,       3083.0,       2833.6,       3651.8,       5565.2,       7738.2,       11014.4,      12082.2,      14264.4,     16435.2,     // partial molar volume v2
                                  0,            11.7721,      11.9407,      24.9559,      9.7119,       2.1778,      -6.3316,      -14.9298,     -23.4091,     -32.0202,    -40.5342,     // partial molar volume v3
                                  0,           -129198,      -129018,      -127989,      -131365,      -139277,      -148260,      -157256,      -166218,      -175238,     -184213 };    // partial molar volume v4

    // std::vector<char> type = {'w', 'm', 'm', 'm', 'm', 'i'}; // w(ater), m(olecular), i(onic)
    int z_a = 1; int z_c = -1;
    double e = 1.60218E-19; double eps0 = 8.85419E-12; double N_A = 6.022E23; double Mw = 18.015;
    int nc = omega.size();
}

AQ2::AQ2(std::vector<int> ci_) {
    ci = ci_; NC = ci.size();
    gi = std::vector<double>(NC);
    hi = std::vector<double>(NC);
    vi = std::vector<double>(NC);
    gi0 = std::vector<double>(NC);

    water_index = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 0));

    // salt = std::find(std::begin(ci), std::end(ci), 5) != std::end(ci); // check if salt is present
    salt = 0;

}

void AQ2::init(double p, double T) {
    pp = p; TT = T; // initialize for reference p & T in integrals
    for (int i = 0; i < NC; i++) {
        index = ci[i];
        
        // ideal gas Gibbs energy
        double hio = Aq2_integrals(aq2_par::T_0, T, 40, 0);
        double gio = aq2_par::gi_00[index]/(aq2_par::R * aq2_par::T_0);
        gi0[i] = gio - hio;

        // Gibbs energy of aqueous phase
        gi[i] = aq2_par::gi_0[index]/(aq2_par::R*aq2_par::T_0);
        hi[i] = Aq2_integrals(aq2_par::T_0, T, 20, 1);
        vi[i] = Aq2_integrals(aq2_par::P_0, p, 20, 2);        
    }

    if (salt) {
        salt_index = std::distance(ci.begin(), std::find(ci.begin(), ci.end(), 5)); //
        B_ca = -0.554860699 + 4.2795E-3*TT - 6.529E-6*pow(TT, 2);
        C_ca = -0.016131327 - 1.25089E-5*TT + 5.89E-8*pow(TT, 2);
        D_ca = -1.12161E-3 + 2.49474E-5*TT - 4.603E-8*pow(TT, 2);
        eps = (243.9576 + 0.039037 * pp - 1.01261E-5 * pow(pp, 2)) + (-0.7520846 - 2.12309E-4 * pp + 6.04961E-8 * pow(pp, 2)) * TT + (6.60648E-4 + 3.18021E-7 * pp - 9.33341E-11 * pow(pp, 2)) * pow(TT, 2);
        double rho_s = 1000;
        A_DH = pow(pow(aq2_par::e, 2)/(aq2_par::eps0*eps*aq2_par::R*TT), 1.5) * pow(aq2_par::N_A, 2)/(8*M_PI)*sqrt(2*rho_s); // kg^0.5/mol^0.5
    }
}

double AQ2::xt(double x) {
    const double eps1 = (243.9576 + 0.039037 * pp - 1.01261E-5 * pow(pp, 2)) + (-0.7520846 - 2.12309E-4 * pp + 6.04961E-8 * pow(pp, 2)) * x + (6.60648E-4 + 3.18021E-7 * pp - 9.33341E-11 * pow(pp, 2)) * pow(x, 2);
    const double d2ln = (eps*2*(6.60648E-4 + 3.18021E-7*pp - 9.33341E-11*pow(pp, 2)) - (-0.7520846 - 2.12309E-4*pp + 6.04961E-8*pow(pp, 2) + 2.0*x*(6.60648E-4 + 3.18021E-7*pp - 9.33341E-11*pow(pp, 2))))/(pow(eps, 2));
    const double dln2 = pow((-0.7520846 - 2.12309E-4*pp + 6.04961E-8*pow(pp, 2) + 2*x*(6.60648E-4 + 3.18021E-7*pp - 9.33341E-11*pow(pp, 2))), 2)/(pow(eps, 2));
    return 1/eps1 * (d2ln - dln2) * x;
}

double AQ2::h_io(double x) {
    double hi0 = aq2_par::hi_00[index] + aq2_par::hi_a[0*aq2_par::nc + index]*(x-aq2_par::T_0) + 1/2.0 * aq2_par::hi_a[1*aq2_par::nc + index]*(pow(x, 2)-pow(aq2_par::T_0, 2)) + 1/3.0 * aq2_par::hi_a[2*aq2_par::nc + index]*(pow(x, 3)-pow(aq2_par::T_0, 3)) + 1/4.0 * aq2_par::hi_a[3*aq2_par::nc + index]*(pow(x, 4)-pow(aq2_par::T_0, 4));
    return hi0 / (aq2_par::R * pow(x, 2));
}

double AQ2::h_i(double x) {
    if (index == 0) { // pure water enthalpy
        double hw = aq2_par::hi_0[index] + 8.712 * aq2_par::R * (x - aq2_par::T_0) + 1 / 2.0 * 0.125E-2 * aq2_par::R * (pow(x, 2) - pow(aq2_par::T_0, 2)) - 1 / 3.0 * 0.018E-5 * aq2_par::R * (pow(x, 3) - pow(aq2_par::T_0, 3));
        return hw / (aq2_par::R * pow(x, 2));
    }
    else { // solute molar enthalpy
        double xt = Aq2_integrals(aq2_par::T_0, x, 1, 3);
        double int_cp = aq2_par::cp[0*aq2_par::nc + index]*(x-aq2_par::T_0) - aq2_par::cp[1*aq2_par::nc + index]*(1/(x-228)-1/(aq2_par::T_0-228)) + aq2_par::omega[index]*xt;
        return (aq2_par::hi_0[index] + int_cp)/(aq2_par::R * pow(x, 2));
    }
}

double AQ2::v_i(double x) {
    if (index == 0) {
        double vw = ((31.1251 - 2.46176E-2 * x + 8.69425E-6 * pow(x, 2) - 6.03348E-10 * pow(x, 3)) +
                    (-1.14154E-1 + 2.15663E-4 * x - 7.96939E-8 * pow(x, 2) + 5.57791E-12 * pow(x, 3)) * TT +
                    (3.10034E-4 - 6.48160E-7 * x + 2.45391E-10 * pow(x, 2) - 1.72577E-14 * pow(x, 3)) * pow(TT, 2) +
                    (-2.48318E-7 + 6.47521E-10 * x - 2.51773E-13 * pow(x, 2) + 1.77978E-17 * pow(x, 3)) * pow(TT, 3)) * 1E-6;
        return vw / (aq2_par::R * 1E-5 * TT);
    }
    else {
        double eps1 = (243.9576 + 0.039037 * x - 1.01261E-5 * pow(x, 2)) + (-0.7520846 - 2.12309E-4 * x + 6.04961E-8 * pow(x, 2)) * TT + (6.60648E-4 + 3.18021E-7 * x - 9.33341E-11 * pow(x, 2)) * pow(TT, 2);
        double dedp = 0.039037 - 2.12309E-4*TT + 3.18021E-7*pow(TT, 2) + 2*x*(-1.01261E-5 + 6.04961E-8*TT - 9.33341E-11*pow(TT, 2));
        double vi_ = (aq2_par::vp[0*aq2_par::nc + index] + aq2_par::vp[1*aq2_par::nc + index]/(2600+x) + (aq2_par::vp[2*aq2_par::nc + index] + aq2_par::vp[3*aq2_par::nc + index]/(2600+x))*1/(TT-228) - aq2_par::omega[index]/(pow(eps1, 2)) * dedp);
        return vi_ / (aq2_par::R * TT);
    }
}

double AQ2::Aq2_integrals(double a, double b, int steps, int function) {
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
    case 1: // hi
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = T
            s += h*((h_i(x) + 4*h_i(x+h/2) + h_i(x+h))/6);
        };
        break;
    case 2: // vi
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = p
            s += h*((v_i(x) + 4*v_i(x+h/2) + v_i(x+h))/6);
        };
        break;
    case 3: // XT integral (for h_i)
        for (int i = 0; i < steps; i++) {
            x = a + h*i; // x = T
            s += h*((xt(x) + 4*xt(x+h/2) + xt(x+h))/6);
        }
        break;
    }
    return s;
}

std::vector<double> AQ2::ln_a(std::vector<double> m, std::vector<int> cii) {
    // Activity of water and solutes
    std::vector<double> lna(NC, 0);
    
    double I{ 0 };
    if (salt) {
        I = 0.5*(2*m[salt_index]);
    }
    for (int i = 0; i < NC; i++) {
        if (cii[i] == 0) { // H2O
            lna[i] = 0;
            if (salt) {
                // Ionic contribution to water activity
                double gamma_Bca = (0.13816 + 0.6*B_ca)*I*abs(aq2_par::z_c*aq2_par::z_a)/1.5 * ((1+3*I/abs(aq2_par::z_c*aq2_par::z_a))/pow(1+3*I/(2*abs(aq2_par::z_c*aq2_par::z_a)), 2) - log(1+3*I/(2*abs(aq2_par::z_c*aq2_par::z_a)))/(3*I/(2*abs(aq2_par::z_c*aq2_par::z_a)))) + 2/abs(aq2_par::z_c*aq2_par::z_a) * (1/2.0*B_ca*pow(I, 2) + 2/3.0*C_ca*pow(I, 3) + 3/4.0*D_ca*pow(I, 4));
                double gamma_DH = 2*A_DH * ((1-pow((1+sqrt(I)), 2))/(1+sqrt(I)) + 2*log(1+sqrt(I)));
                double gamma_Ica = -aq2_par::Mw*(2*I/abs(aq2_par::z_c*aq2_par::z_a) + gamma_DH + gamma_Bca);
                double ionic = m[salt_index]*m[salt_index]*gamma_Ica;
                ionic /= (m[salt_index]*pow(aq2_par::z_c, 2) * m[salt_index]*pow(aq2_par::z_a, 2));
                lna[i] += ionic;
            }
            // Contribution of molecular/ionic species to water activity
            double molecular{ 0 };
            for (int j = 0; j < NC; j++) { // all molecular/ionic species
                for (int k = 0; k < NC; k++) {
                    if (((cii[j] == 1) && (cii[k] == 5)) || ((cii[j] == 5) && (cii[k] == 1))) { // CH4-NaCl
                        double gamma_P1 = 0.025 - 5E-3/(2*I) * (1-(1+2*sqrt(I))*exp(-2*sqrt(I)));
                        double gamma_P2 = -5E-3*(1-(1+2*sqrt(I) + 2*I)*exp(-2*sqrt(I)));
                        molecular += 2*m[j]*m[k]*(gamma_P1-gamma_P2/(2*I));
                    }
                    else if ((cii[j] == 2) && (cii[k] == 2)) { // CO2-CO2
                        double gamma_P1 = 0.107 - 4.5E-4*TT;
                        molecular += m[j]*m[k]*gamma_P1;
                    }
                }
                if ((cii[j] != 0) && (cii[j] != 5)) { // all molecular species
                    molecular += m[j];
                }                
            }
            lna[i] -= aq2_par::Mw/1000*molecular;
        }
        else {
            double lnj;
            double gamma_P1{ 0 };
            if ((salt) && (cii[i] == 1)) {
                gamma_P1 = 0.025 - 5E-3/(2*I) * (1-(1+2*sqrt(I))*exp(-2*sqrt(I)));
            }
            else if (ci[i] == 2) {
                gamma_P1 = 0.107 - 4.5E-4*TT;
            }
            lnj = 2*m[i]*gamma_P1;
            lna[i] += log(m[i]*exp(lnj));
        }
    }
    
    return lna;
}

std::vector<double> AQ2::fugacityCoefficient(double p, double T, std::vector<double> x) {
    // Construct fugacity coefficients
    std::vector<double> phi_ij(NC);
    T = T;
    // Calculate molality
    std::vector<double> m(NC);
    for (int i = 0; i < NC; i++) {
        m[i] = 55.509 * x[i] / x[water_index];
    }

    std::vector<double> lna = ln_a(m, ci);
    
    for (int i = 0; i < NC; i++) {
        double mu = gi[i] - hi[i] + vi[i] + lna[i];
        phi_ij[i] = 1 * exp(mu - gi0[i])/(x[i]*p);
    }


    
    return phi_ij;
}