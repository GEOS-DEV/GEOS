#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>
#include <unordered_map>

#include "eos.h"
#include "../global/global.h"

using namespace std;

// namespace aq2_par {
//     double R = 8.3145; double T_0 = 298.15; double P_0 = 1.0;
//     // gibbs energy of pure H2O or 1 molal solution at p0, T0
//     std::unordered_map<std::string, double> gi_0 =  {{"H2O", -237129}, {"CO2", -385974}, {"N2", 18188}, {"H2S", -27920}, {"C1", -34451}, {"C2", -17000}, {"C3", -7550}, {"nC4", -940}, {"nC5", 9160}, {"nC6", 18200}, {"nC7", 27500} };
//     // molar enthalpy of pure H2O or 1 molal solution at p0, T0
//     std::unordered_map<std::string, double> hi_0 =  {{"H2O", -285830}, {"CO2", -413798}, {"N2", -10439}, {"H2S", -37660}, {"C1", -87906}, {"C2", -103136}, {"C3", -131000}, {"nC4", -152000}, {"nC5", -173887}, {"nC6", -199200}, {"nC7", -225000} };
//     // gibbs energy of ideal gas at p0, T0
//     std::unordered_map<std::string, double> gi_00 = {{"H2O", -228700}, {"CO2", -394600}, {"N2", 0}, {"H2S", -33100}, {"C1", -50830}, {"C2", -32900}, {"C3", -23500}, {"nC4", -17200}, {"nC5", -8370}, {"nC6", -290}, {"nC7", 8120} };
//     // molar enthalpy of ideal gas at p0, T0
//     std::unordered_map<std::string, double> hi_00 = {{"H2O", -242000}, {"CO2", -393800}, {"N2", 0}, {"H2S", -20200}, {"C1", -74900}, {"C2", -84720}, {"C3", -103900}, {"nC4", -126200}, {"nC5", -146500}, {"nC6", -167300}, {"nC7", -187900} };
//     // ideal gas heat capacity parameters [eq. 3.4]
//     std::unordered_map<std::string, std::vector<double>> hi_a =  {
//         {"H2O", {3.8747*R, 0.0231E-2*R, 0.1269E-5*R, -0.4321E-9*R}},
//         {"CO2", {2.6751*R, 0.7188E-2*R, -0.4208E-5*R, 0.8977E-9*R}},
//         {"N2", {3.4736*R, -0.0189E-2*R, 0.0971E-5*R, -0.3453E-9*R}},
//         {"H2S", {3.5577*R, 0.1574E-2*R, 0.0686E-5*R, -0.3959E-9*R}},
//         {"C1", {2.3902*R, 0.6039E-2*R, 0.1525E-5*R, -1.3234E-9*R}},
//         {"C2", {0.8293*R, 2.0752E-2*R, -0.7699E-5*R, 0.8756E-9*R}},
//         {"C3", {-0.4861*R, 3.6629E-2*R, -1.8895E-5*R, 3.8143E-9*R}},
//         {"nC4", {0.4755*R, 4.4650E-2*R, -2.2041E-5*R, 4.2068E-9*R}},
//         {"nC5", {0.8142*R, 5.4598E-2*R, -2.6997E-5*R, 5.0824E-9*R}},
//         {"nC6", {0.8338*R, 6.6373E-2*R, -3.444E-5*R, 6.9342E-9*R}},
//         {"nC7", {-0.6184*R, 8.1268E-2*R, -4.388E-5*R, 9.2037E-9*R}},
//     };
//     // Born constants of solutes [eq. 3.6-3.7]
//     std::unordered_map<std::string, double> omega = {{"CO2", -8368}, {"N2", -145101}, {"H2S", -41840}, {"C1", -133009}, {"C2", -169870}, {"C3", -211418}, {"nC4", -253592}, {"nC5", -300955}, {"nC6", -335180}, {"nC7", -380158} };
//     // Partial molar heat capacity terms [eq. 3.6]
//     std::unordered_map<std::string, double> cp1 = {{"CO2", 167.50}, {"N2", 149.75}, {"H2S", 135.14}, {"C1", 176.12}, {"C2", 226.67}, {"C3", 277.52}, {"nC4", 330.77}, {"nC5", 373.24}, {"nC6", 424.53}, {"nC7", 472.37}, };
//     std::unordered_map<std::string, double> cp2 = {{"CO2", 5304066}, {"N2", 5046230}, {"H2S", 2850801}, {"C1", 6310762}, {"C2", 9011737}, {"C3", 11749531}, {"nC4", 14610096}, {"nC5", 16955051}, {"nC6", 19680558}, {"nC7", 22283347}, };
//     // Partial molar volume terms [eq. 3.7]
//     std::unordered_map<std::string, std::vector<double>> vp = {
//         {"CO2", {2.614, 3125.9, 11.7721, -129198}},
//         {"N2", {2.596, 3083.0, 11.9407, -129018}},
//         {"H2S", {2.724, 2833.6, 24.9559, -127989}},
//         {"C1", {2.829, 3651.8, 9.7119, -131365}},
//         {"C2", {3.612, 5565.2, 2.1778, -139277}},
//         {"C3", {4.503, 7738.2, -6.3316, -148260}},
//         {"nC4", {5.500, 11014.4, -14.9298, -157256}},
//         {"nC5", {6.282, 12082.2, -23.4091, -166218}},
//         {"nC6", {7.175, 14264.4, -32.0202, -175238}},
//         {"nC7", {8.064, 16435.2, -40.5342, -184213}},
//     };
//     // Dielectric constant of water coefficients [eq. 3.9]
//     std::vector<double> eps = {243.9576, 0.039037, -1.01261E-5, 
//                             -0.7520846, -2.12309E-4, 6.04961E-8, 
//                             6.60648E-4, 3.18021E-7, -9.33341E-11};

//     // Solute interaction parameters eq. 3.18-3.20
//     std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Bca = {
//         {"Na+", {{"Cl-", {-0.554860699, 4.2795E-3, -6.529E-6}}}},
//         {"K+", {{"Cl-", {0.178544751, -9.55043E-4, 1.8208E-6}}}},
//         {"Ca2+", {{"Cl-", {0.549244833, -1.870735E-3, 3.3604E-6}}}},
//     };
//     std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Cca = {
//         {"Na+", {{"Cl-", {-0.016131327, -1.25089E-5, 5.89E-8}}}},
//         {"K+", {{"Cl-", {-5.546927E-3, 4.22294E-5, -9.038E-8}}}},
//         {"Ca2+", {{"Cl-", {-0.011031685, 7.49491E-5, -1.639E-7}}}},
//     };
//     std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> Dca = {
//         {"Na+", {{"Cl-", {-1.12161E-3, 2.49474E-5, -4.603E-8}}}},
//         {"K+", {{"Cl-", {7.12650E-5, -6.04659E-7, 1.327E-9}}}},
//         {"Ca2+", {{"Cl-", {1.08383E-4, -1.03524E-6, 2.3878E-9}}}},
//     };
//     std::unordered_map<std::string, std::unordered_map<std::string, std::vector<double>>> B = {
//         {"H2S", {{"H2S", {0.1, -4.7e-4, 0., 0.}}}},
//         {"CO2", {{"CO2", {0.107, -4.5e-4, 0., 0.}}}},
//         {"C1", {{"Na+", {0.025, 0., 0. -5e-3}}, {"Cl-", {0.025, 0., 0. -5e-3}}}},
//         {"C3", {{"Na+", {-0.09809, 4.19e-4, -6.2e-6, 0.}}, {"Cl-", {-0.09809, 4.19e-4, -6.2e-6, 0.}}}},
//     };

//     double e = 1.60218E-19; double eps0 = 8.85419E-12;
// }

void AQ2::parameters(double p_, double T_) {
    p = p_; pp = p_; 
    T = T_; TT = T_; // initialize for reference p & T in integrals

    gi = std::vector<double>(NC);
    hi = std::vector<double>(NC);
    vi = std::vector<double>(NC);
    gi0 = std::vector<double>(NC);

    water_index = std::distance(components.begin(), std::find(components.begin(), components.end(), "H2O"));

    for (int i = 0; i < NC; i++) 
    {
        comp = components[i];
        // ideal gas Gibbs energy
        double hio = Aq2_integrals(aq2_par::T_0, T, 10, 0);  // eq. 3.3
        double gio = aq2_par::gi_00[comp]/(aq2_par::R * aq2_par::T_0);
        gi0[i] = gio - hio;  // eq. 3.2

        // Gibbs energy of aqueous phase
        gi[i] = aq2_par::gi_0[comp]/(aq2_par::R*aq2_par::T_0);
        hi[i] = Aq2_integrals(aq2_par::T_0, T, 10, 1);
        vi[i] = Aq2_integrals(aq2_par::P_0, p, 10, 2);
    }

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
    
    return;
}

double AQ2::xt(double x) {
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

double AQ2::Aq2_functions(double x, int function) {
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

double AQ2::Aq2_integrals(double a, double b, int steps, int function) {
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

std::vector<double> AQ2::ln_a(std::vector<double> m_s) {
    // Activity of water and solutes
    std::vector<double> lna(NC, 0.);
    
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
    for (int i = 0; i < NC; i++)
    {
        if (components[i] == "H2O")
        {
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
                lna[i] += numerator/(c_denominator*a_denominator);
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
            lna[i] -= CompProp::prop["Mw"]["H2O"]/1000*molecular;
        }
        else
        {
            // m-m and m-i interactions
            double lnj = 0;
            for (int k = 0; k < NS; k++)
            {
                lnj += 2 * m_s[k] * gamma_P1[i*NS + k];
            }
            lna[i] += log(m_s[i] * exp(lnj));
        }
    }
    return lna;
}

std::vector<double> AQ2::fugacityCoefficient(std::vector<double> x) {
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

    std::vector<double> lna = ln_a(m_s);
    
    for (int i = 0; i < NC; i++) 
    {
        double mu = gi[i] - hi[i] + vi[i] + lna[i];  // eq. 3.5/3.10
        phi[i] = 1 * exp(mu - gi0[i])/(x[i]*p);
    }
    
    return phi;
}