#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "eos.h"
#include "../global/global.h"

using namespace std;

void TwoParCubic::parameters(double p_, double T_) {
    // calculate x - independent part of parameters
    p = p_; T = T_;

    Ai = std::vector<double>(NC);
    Bi = std::vector<double>(NC);
    for (int i = 0; i < NC; i++) 
    {
        double p_r = p / CompProp::prop["Pc"][components[i]];
        double T_r = T / CompProp::prop["Tc"][components[i]];
        Ai[i] = omegaA * p_r / pow(T_r, 2) * pow(1. + kappa[i] * (1 - sqrt(T_r)), 2);
        Bi[i] = omegaB * p_r / T_r;
    }

    Aij = std::vector<double>(NC*NC);
    for (int i = 0; i < NC; i++) 
    {
        for (int j = 0; j < NC; j++) 
        {
            double kij = CompProp::bic[components[i]][components[j]];
            Aij[i*NC + j] = sqrt(Ai[i] * Ai[j]) * (1 - kij);
        }
    }

    return;
}

std::vector<std::complex<double>> TwoParCubic::cubic_roots(double a, double b, double c, double d) {
    // Following Srikanth (2020) - A note on the solutions of cubic equations of state in low temperature region
    b /= a; c /= a; d /= a;
    double p = c - pow(b, 2)/3.;
    double q = d - b * c / 3. + 2 * pow(b, 3) / 27;
    double D = pow(p, 3) / 27 + pow(q, 2) / 4;

    std::vector<std::complex<double>> Z(3);

    if (D > 0.) 
    {
        std::vector<double> real(3);
        std::vector<double> imag(3);

        double u = std::cbrt(-q/2. + sqrt(D));
        double v = std::cbrt(-q/2. - sqrt(D));
        std::complex<double> w(-0.5, 0.5*sqrt(3));
        std::complex<double> w2 = std::conj(w);

        Z[0] = u + v - b/3.;
        Z[1] = w*u + w2*v - b/3.;
        Z[2] = w2*u + w*v - b/3.;
    } 
    else 
    {
        for (int k = 0; k < 3; k++) 
        {
            double theta = acos(-q/2 * sqrt(-27 / pow(p, 3)));
            Z[k] = -b / 3. + 2 * sqrt(-p / 3) * cos((theta + 2*M_PI*k)/3);
        }
    }

    return Z;
}

double TwoParCubic::calc_z(std::vector<double> x) {
    // Calculate x-dependent parameters a, b, A and B
    A = 0.;
    B = 0.;
    for (int i = 0; i < NC; i++) 
    {
        for (int j = 0; j < NC; j++) 
        {
            A += x[i] * x[j] * Aij[i*NC + j];
        }
        B += x[i] * Bi[i];
    }

    // Find roots of cubic EoS
    double a = 1;
    double b = (d1 + d2 - 1) * B - 1;
    double c = A + d1 * d2 * pow(B, 2) - (d1 + d2) * B * (B + 1);
    double d = - (A * B + d1 * d2 * pow(B, 2) * (B + 1));
    std::vector<std::complex<double>> Z = cubic_roots(a, b, c, d);

    // Check if three real roots
    bool real;
    if (Z[2].imag() != 0.) 
    { 
        real = false; 
    }
    else 
    { 
        real = true; 
    }

    double z;
    if (!real) 
    {
        // only one real root
        z = Z[0].real();
        
        double gI = 0;
        for (int i = 0; i < NC; i++) 
        { 
            gI += x[i] * log(x[i]); 
        }

        double gE = (z - 1) - log(z - B) - A / (B * (d1 - d2)) * log((z + d1 * B) / (z + d2 * B));
        
        g = gI + gE;
    } 
    else 
    {
        // Find zmin and zmax
        double z_min = Z[0].real();
        if (Z[1].real() < z_min) { z_min = Z[1].real(); }
        if (Z[2].real() < z_min) { z_min = Z[2].real(); }
        double z_max = Z[0].real();
        if (Z[1].real() > z_max) { z_max = Z[1].real(); }
        if (Z[2].real() > z_max) { z_max = Z[2].real(); }

        double gI = 0;
        for (int i = 0; i < NC; i++)
        { 
            gI += x[i] * log(x[i]); 
        }
        
        if (z_min <= B) 
        {
            z = z_max;

            double gI = 0;
            for (int i = 0; i < NC; i++)
            { 
                gI += x[i] * log(x[i]); 
            }
            
            double gE = (z - 1) - log(z - B) - A / (B * (d1 - d2)) * log((z + d1 * B) / (z + d2 * B));
            
            g = gI + gE;
        } else {
            double gI = 0;
            for (int i = 0; i < NC; i++)
            { 
                gI += x[i] * log(x[i]); 
            }
                    
            double gE_l = (z_min - 1) - log(z_min - B) - A / (B * (d1 - d2)) * log((z_min + d1 * B) / (z_min + d2 * B));
            double gE_v = (z_max - 1) - log(z_max - B) - A / (B * (d1 - d2)) * log((z_max + d1 * B) / (z_max + d2 * B));
            
            if (gE_v < gE_l) 
            {
                z = z_max;
                g = gI + gE_v;
            } 
            else 
            {
                z = z_min;
                g = gI + gE_l;
            }                        
        }
    }

    return z;
}

std::vector<double> TwoParCubic::fugacityCoefficient(std::vector<double> x) {
    double z = calc_z(x);

    // Calculate fugacity coefficient
    std::vector<double> phi(NC);
    for (int k = 0; k < NC; k++) 
    {
        double psi = 0;
        for (int j = 0; j < NC; j++) 
        {
            psi += x[j] * Aij[k*NC + j];
        }

        phi[k] = exp((z - 1) * Bi[k] / B - log(z - B) - A / ((d1 - d2) * B) * (2 * psi / A - Bi[k] / B) * log((z + d1 * B) / (z + d2 * B)));
    }
    return phi;
}

double TwoParCubic::calc_H(double p_, double T_, std::vector<double> x) {
    p = p_; T = T_;
    double Tf = T + 0.001;
    
    // calc da/dT with forward differentiation
    Ai = std::vector<double>(NC);
    Bi = std::vector<double>(NC);
    std::vector<double> Afi(NC);
    // std::vector<double> ai(NC);
    // std::vector<double> dadT(NC);
    for (int i = 0; i < NC; i++) 
    {
        double p_r = p / CompProp::prop["Pc"][components[i]];
        double T_r = T / CompProp::prop["Tc"][components[i]];
        double Tf_r = Tf / CompProp::prop["Tc"][components[i]];
        Ai[i] = omegaA * p_r / pow(T_r, 2) * pow(1 + kappa[i] * (1 - sqrt(T_r)), 2);
        Afi[i] = omegaA * p_r / pow(Tf_r, 2) * pow(1 + kappa[i] * (1 - sqrt(Tf_r)), 2);
        Bi[i] = omegaB * p_r / T_r;

        // ai[i] = omegaA * pow(R, 2) * pow(CompProp::Tc[components[i]], 2) / CompProp::Pc[components[i]];
        // double alpha = pow(1+kappa[i]*(1-sqrt(T_r)), 2);
        // dadT[i] = -kappa[i] * ai[i] * sqrt(alpha/(T*CompProp::Tc[components[i]]));
    }

    Aij = std::vector<double>(NC*NC);
    A = 0; B = 0;
    double Af = 0;
    for (int i = 0; i < NC; i++) 
    {
        for (int j = 0; j < NC; j++) 
        {
            double kij = CompProp::bic[components[i]][components[j]];
            Aij[i*NC + j] = sqrt(Ai[i] * Ai[j]) * (1 - kij);
            double Afij = sqrt(Afi[i] * Afi[j]) * (1 - kij);

            A += x[i] * x[j] * Aij[i*NC + j];
            Af += x[i] * x[j] * Afij;

            // dAdT += 0.5*x[i]*x[j]*(1-kij)*sqrt(ai[i]*ai[j]) * (dadT[i]/ai[i] + dadT[j]/ai[j]);
        }
        B += x[i] * Bi[i];
    }

    double a = A / p * pow(R, 2) * pow(T, 2);
    double af = Af / p * pow(R, 2) * pow(Tf, 2);
    double b = B / p * R * T;

    double dAdT = (af - a) / (Tf - T);

    double z = calc_z(x);
    
    double H_dev = R*T*(z-1) + (T*dAdT - a) / ((d1-d2)*b) * log((z + d1*B)/(z + d2*B));
    H_dev *= 1e5;

    return H_dev;
}

void PR::init() {
    // PR-EoS parameters
    d1 = 1 + sqrt(2);
    d2 = 1 - sqrt(2);

    omegaA = 0.45724;
    omegaB = 0.0778;

    kappa = std::vector<double>(NC);
    for (int i = 0; i < NC; i++) 
    {
        double ac = CompProp::prop["ac"][components[i]];
        if (ac <= 0.49) 
        {
            kappa[i] = 0.37464 + 1.54226 * ac - 0.26992 * pow(ac, 2);
        } 
        else 
        {
            kappa[i] = 0.379642 + 1.48503 * ac - 0.164423 * pow(ac, 2) + 0.016667 * pow(ac, 3);
        }
    }
    return;
}

void SRK::init() {   
    // SRK-EoS parameters
    d1 = 0;
    d2 = 1;

    omegaA = 0.42748;
    omegaB = 0.08664;

    kappa = std::vector<double>(NC);
    for (int i = 0; i < NC; i++) 
    {
        double ac = CompProp::prop["ac"][components[i]];
        kappa[i] = 0.48 + 1.574 * ac - 0.176 * pow(ac, 2);
    }
    return;
}