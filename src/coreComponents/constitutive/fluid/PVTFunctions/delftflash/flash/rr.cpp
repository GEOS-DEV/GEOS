#include <iostream>
#include <cmath>
#include <algorithm>

#include "rr.h"
#include "../global/misc.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

RR::RR(std::vector<double> z_, std::vector<double> K_) {
	z = z_; K = K_;
	NC = z.size(); NP = K.size()/NC + 1;

	v_min = std::vector<double>(NP-1);
	v_max = std::vector<double>(NP-1);

	// Check if domain is bounded: admissible K-values (gji*gki<0)
	admissible = bounded();
}

std::vector<double> RR::getx(std::vector<double> V) {
	std::vector<double> x(NP*NC);

	for (int i = 0; i < NC; i++) 
	{
		double m_i = 1;
		for (int k = 1; k < NP; k++) 
		{
			m_i += V[k] * (K[(k-1)*NC + i] - 1);
		}
		x[i] = z[i] / m_i;
	}

	for (int j = 1; j < NP; j++) 
	{
		for (int i = 0; i < NC; i++) 
		{
			x[j*NC + i] = K[(j-1)*NC + i] * x[i];
		}
	}

	return x;
}

bool RR::bounded() {
	// Calculate directional derivative of mi(V) along mj(V)
	// gij = Dm_i dotted with u_j
	// Dm_i[i, k] = K[i, k] - 1
	// u_j[j] = [1, 1, ..., sum((K[j, k] - 1)/(1-K[j, np-1]))]

	std::vector<double> Dm_i(NP-1, 0.); 
	std::vector<double> u_j(NP-1, 1.);
	std::vector<double> gij(NC*NC, 0.);
	for (int i = 0; i < NC; i++) 
	{
		for (int k = 0; k < NP-1; k++) 
		{
			Dm_i[k] = K[k*NC + i] - 1;
		}
		for (int j = 0; j < NC; j++) 
		{
			if (NP > 2) 
			{
				u_j[NP-2] = 0;
				for (int k = 0; k < NP-2; k++) 
				{
					u_j[NP-2] += (K[k*NC + j] - 1)/(1-K[(NP-2)*NC + j]);
				}
			}
			for (int k = 0; k < NP-1; k++) 
			{
				gij[i*NC + j] += Dm_i[k]*u_j[k];
			}
		}
	}

	// Check if there exist j, k such that gji*gki < 0 for i = 0, ..., NC-1
	std::vector<bool> gi(NC, false);
	for (int i = 0; i < NC; i++) 
	{
		for (int j = 0; j < NC; j++) 
		{
			for (int k = 0; k < NC; k++) 
			{
				if (gij[j*NC+i] * gij[k*NC+i] < 0) 
				{
					gi[i] = true;
					k = NC; j = NC; // go to next i in loop
				}
			}
		}
		if (gi[i] == false) 
		{ 
			return false;  // there exist no j, k such that gji*gki < 0 for this i, so region is not bounded -> exit loop
		}
	}

	return true;
}

std::vector<double> RR::V_limits(std::vector<double> V_j, int J) {
	// Calculate V_min and V_max from corner points of domain

	// This function calculates all intersections of planes m_i(v) and then checks which ones are actual corner points and range of the domain for V
	// First calculate total number of intersections between mi and mj: NC!/((NP-1)!(NC-(NP-1))!)
	int ni = factorial(NC)/(factorial(J)*factorial(NC-J));
	std::vector<double> intersections(J*ni);
	
	// Find all unique combinations of i of length NP-1
	Combinations c(ni, NC, J);  // contains vector with (NP-1)*(number of intersections) entries, stores combinations of i for unique intersections of m_i(v)
	
	// Solve linear system with Eigen library to find all V_j at intersection of m_i's: Ax = b
	// m_i = 1 + sum(K_ik-1)*v_k = 0 -> sum(K_ik-1)*v_k = -1 for each m_i
	// A_ik = [K_ik-1]; b_k = [-1, ..., -1]
	Eigen::MatrixXd A(J, J);
	Eigen::VectorXd b(J);
	Eigen::VectorXd v(J);
	for (int j = 0; j < J; j++) 
	{ 
		b(j) = -1; 
	}

	int index;  // index i from combination
	for (int i = 0; i < ni; i++) 
	{ 
		// Fill matrix A for combination of m_i's
		for (int j = 0; j < J; j++) 
		{
			index = c.getIndex(i, j); // index returns component j for i-th combination
			for (int k = 0; k < J; k++) 
			{  
				// loop over phase k to fill row j
				A(j, k) = K[k*NC + index] - 1; // A_jk = K_jk-1
			}
			for (int k = J; k < NP-1; k++) { 
				b(j) -= V_j[k]*(K[k*NC + index] - 1);  // if V_j+ are known
			}
		}
		v = A.partialPivLu().solve(b);  // with LU factorization
		for (int j = 0; j < J; j++) { 
			intersections[J*i + j] = v(j); 
			b(j) = -1; 
		} 
	}

	// If all but the three intersecting lines for m_i(v) > 0, then this intersection is a cornerpoint of the domain for V
	// Calculate m_i at intersections for all i = 0, ..., NC-1
	std::vector<double> Vlim(2*J);
	std::vector<double> m_i;
	int count_corners = 0;
	for (int i = 0; i < ni; i++) 
	{
		m_i = std::vector<double>(NC, 1.); // value of m_i's at the ni intersections = 1 + sum((K_kj-1)*v_j)
		for (int k = 0; k < NC; k++) 
		{
			for (int j = 0; j < J; j++) 
			{
				m_i[k] += intersections[J*i + j] * (K[j*NC + k]-1); // += v_j*(K_kj-1)
			}
			for (int j = J; j < NP-1; j++) 
			{
				m_i[k] += V_j[j] * (K[j*NC + k]-1);
			}
		}
		// Count the number of m_i's > 0
		int count = 0;
		for (int k = 0; k < NC; k++) 
		{ 
			if (m_i[k] > 1E-16) 
			{ 
				count++; 
			} 
		}
		if (count >= NC-J) 
		{ 
			count_corners++;
			for (int j = 0; j < J; j++) 
			{
				// v_j, min
				if (( intersections[J*i + j] < Vlim[2*j] ) || (count_corners == 1))
				{ 
					Vlim[2*j] = intersections[J*i + j]; 
				}
				// v_j, max
				if (( intersections[J*i + j] > Vlim[2*j+1] ) || (count_corners == 1))
				{ 
					Vlim[2*j+1] = intersections[J*i + j];
				}
			}
		}
	}
	if (count_corners - J < 1) 
	{ 
		cout << "Not enough corner points found: unable to calculate limits in J dimensions" << endl;
	}

	return Vlim;
}

std::vector<double> RR::f_df1(std::vector<double> V_j, int J) {
	std::vector<double> F(2, 0.);

	double m_i;
    for (int i = 0; i < NC; i++) 
	{
		m_i = 1;
		for (int k = 0; k < NP-1; k++) 
		{
			m_i += V_j[k] * (K[k*NC + i] - 1);
		}
		F[0] += z[i] * (1 - K[J*NC + i]) / m_i; // f_j
		F[1] -= z[i] * (1 - K[J*NC + i]) / pow(m_i, 2.) * (K[J*NC + i] - 1); // df_j/dv_j
	}
	
	return F;
}

// RRn::RRn(std::vector<double> z_, std::vector<double> K_) : RR(z_, K_) {
// 	Eigen::MatrixXd df = Eigen::MatrixXd::Zero(NP-1, NP-1);
// 	Eigen::VectorXd f = Eigen::VectorXd::Zero(NP-1);
// }

RRbn::RRbn(std::vector<double> z_, std::vector<double> K_) : RR(z_, K_) {
	// Initializing v and f
	v = std::vector<double>(NP-1);
	f = std::vector<double>(NP-1);
}

// RRb::RRb(std::vector<double> z_, std::vector<double> K_) : RR(z_, K_) {
// 	// Initializing v and f
// 	v_mid = std::vector<double>(NP-1);
// 	f_min = std::vector<double>(NP-1);
// 	f_mid = std::vector<double>(NP-1);
// }

/*
std::vector<double> RRn::solveRR(std::vector<double> K_) {
	// Solve Rachford-Rice equations using Newton's method to find set of phase fractions v_j
	// Equations to solve: f_j = (z_1*(1-K_1j) + ... + z_n*(1-K_nj))/(1 + sum(v_k*(K_ik-1))) for each j = 0, ..., NP-1
	K = K_;

	// Using Eigen library: v = df^-1*f
	Eigen::MatrixXd df = Eigen::MatrixXd::Zero(NP-1, NP-1);
	Eigen::VectorXd f = Eigen::VectorXd::Zero(NP-1);
	Eigen::VectorXd v(NP-1);
	Eigen::VectorXd dv(NP-1);
	for (int j = 0; j < NP-1; j++) { v(j) = v[j]; }  // Calculate v(j) = V_mid = (V_min + V_max)/2 for each phase
	
	// Residual of initial guess
	// double m_i;
	// for (int i = 0; i < NC; i++) {
	// 	m_i = 1.;
	// 	for (int k = 0; k < NP-1; k++) {
	// 		m_i += v(k)*(K[k*NC + i] - 1);
	// 	}
	// 	for (int j = 0; j < NP-1; j++) {
	// 		f(j) += z[i]*(1 - K[j*NC + i]) / m_i;
	// 		for (int k = 0; k < NP-1; k++) {
	// 			df(j, k) -= z[i]*(1 - K[j*NC + i]) / pow(m_i, 2.) * (K[k*NC + i] - 1);
	// 		}
	// 	}
	// }
	// Calculate residual of initial guess (f and df)
	f_df();
	
	// Convergence loop
	while (f.norm() > eps) {
		// iter++;
		// dv = df.colPivHouseholderQr().solve(f);
		dv = df.partialPivLu().solve(f);  // with LU factorization
		v -= dv;

		f_df();
		// for (int j = 0; j < NP-1; j++) { f(j) = 0; for (int k = 0; k < NP-1; k++) { df(j, k) = 0.; } }

		// for (int i = 0; i < NC; i++) {
		// 	m_i = 1.;
		// 	for (int k = 0; k < NP-1; k++) {
		// 		m_i += v(k)*(K[k*NC + i] - 1);
		// 	}
		// 	for (int j = 0; j < NP-1; j++) {
		// 		f(j) += z[i]*(1 - K[j*NC + i]) / m_i;
		// 		for (int k = 0; k < NP-1; k++) {
		// 			df(j, k) -= z[i]*(1 - K[j*NC + i]) / pow(m_i, 2.) * (K[k*NC + i] - 1);
		// 		}
		// 	}
		// }
	}

	// Update V
	V[0] = 1.;
	for (int j = 1; j < NP; j++) { V[j] = v(j-1); V[0] -= V[j]; }
	// cout << "Iterations: " << iter << endl;

	return V;
}

void RRn::f_df() {
	for (int j = 0; j < NP-1; j++) { f(j) = 0; for (int k = 0; k < NP-1; k++) { df(j, k) = 0.; } }
	
	double m_i;
	for (int i = 0; i < NC; i++) {
		m_i = 1.;
		for (int k = 0; k < NP-1; k++) {
			m_i += v(k)*(K[k*NC + i] - 1);
		}
		for (int j = 0; j < NP-1; j++) {
			f(j) += z[i]*(1 - K[j*NC + i]) / m_i;
			for (int k = 0; k < NP-1; k++) {
				df(j, k) -= z[i]*(1 - K[j*NC + i]) / pow(m_i, 2.) * (K[k*NC + i] - 1);
			}
		}
	}
}
*/

std::vector<double> RRbn::solveRR(std::vector<double> K_) {
	converged = std::vector<bool>(NP-1, false);
	
	// Latest K-values
	K = K_;
	for (int j = 1; j < NP; j++) 
	{ 
		for (int i = 0; i < NC; i++) 
		{ 
			if (K[(j-1)*NC + i] < std::numeric_limits<double>::epsilon()) 
			{ 
				cout << "K-value below epsilon" << endl; 
				K[(j-1)*NC + i] = std::numeric_limits<double>::epsilon(); 
			}
		}
	}

	// Output vector, contains all NP values of V
	std::vector<double> V(NP);
	
	// Check if domain is bounded: admissible K-values (gji*gki<0)
	// If not admissible, give error ##
	admissible = bounded();
	if (!admissible) 
	{ 
		cout << "Inadmissible K-values: region is not bounded" << endl;
		return V;
	}
	
	// Calculate limits of V domain
    int J = NP-1;
    std::vector<double> Vlim = V_limits(v, J); // Calculate v_min and v_max of j - based on values for v_mid[j+]
    v_min[J-1] = Vlim[2*(J-1)];
    v_max[J-1] = Vlim[2*(J-1) + 1];
    v[J-1] = (v_min[J-1] + v_max[J-1]) / 2;

    // Run recursive RR loop starting from NP-1 (so index=NP-2)
    rrLoop(J-1);

    // Update V
    V[0] = 1.;
	for (int j = 1; j < NP; j++) 
	{ 
		V[j] = v[j - 1]; V[0] -= V[j]; 
	}

	return V;
}

void RRbn::rrLoop(int J) {
	// This function recursively calculates the value of v[J] for which f_j[v] = 0
	// It is repeated until |f_j[v]| < eps
	// If f_j for j-1 has not yet converged, the function moves down one level to find [v_j-]
	while (!converged[J]) 
	{
		// If j > 0, then move down one level (j - 1)
		if (J > 0) 
		{
			// Guess V_mid for all j- phases
			std::vector<double> Vlim = V_limits(v, J); // Calculate Vmin and Vmax of j - based on values for V_mid[j+]
			v_min[J-1] = Vlim[2*(J-1)];
			v_max[J-1] = Vlim[2*(J-1) + 1];

			// if previous v[j] is out of bounds, recalculate v_mid[j]. Will save iterations if not necessary
			if (!(v_min[J-1] < v[J-1]) || !(v[J-1] < v_max[J-1])) 
			{ 
				v[J-1] = (v_min[J-1] + v_max[J-1]) / 2; 
			}
			
			converged[J-1] = false;

			rrLoop(J-1);
		}

		std::vector<double> F = f_df1(v, J);
		f[J] = F[0]; df = F[1];
		double fdf = f[J] / df;

		// if f[J] < 0 -> f[J] = 0 is above v[J], update v_min[J] to v[J]
		if (f[J] < 0) 
		{
			v_min[J] = v[J];
			if (v[J] - fdf > v_max[J]) 
			{
				// cout << "Correction applied, upper limit" << endl;
				v[J] = (v[J] + v_max[J]) / 2;
			}
			else 
			{
				v[J] -= fdf;
			}
		} 
		else 
		{ // else, f[J] > 0 -> f[J] = 0 is below v[J], update v_max[J] to v[J]
			v_max[J] = v[J];
			if (v[J] - fdf < v_min[J]) 
			{
				// cout << "Correction applied, lower limit" << endl;
				v[J] = (v[J] + v_min[J]) / 2;
			} 
			else 
			{
				v[J] -= fdf;
			}
		}
		if (abs(f[J]) < eps)
		{
			converged[J] = true;
			return;
		}
	}
}

/*
std::vector<double> RRb::solveRR() {
	// Output vector, contains all NP values of V
	V = std::vector<double>(NP);

	// If not admissible, give error ##
	if (!admissible) { 
		cout << "Inadmissible K-values: region is not bounded" << endl;
		return V;
	}
	
	// Calculate limits of V domain
    int J = NP-1;
    std::vector<double> Vlim = V_limits(v_mid, J); // Calculate v_min and v_max of j - based on values for v_mid[j+]
    v_min[J-1] = Vlim[2*(J-1)];
    v_max[J-1] = Vlim[2*(J-1) + 1];
    v_mid[J-1] = (v_min[J-1] + v_max[J-1]) / 2;

    // Run recursive RR loop starting from NP-1 (so index=NP-2)
    rrLoop(v_min[J-1], v_mid[J-1], v_max[J-1], J-1);

    // Update V
    V[0] = 1.;
	for (int j = 1; j < NP; j++) { V[j] = v_mid[j - 1]; V[0] -= V[j]; }

	return V;
}

void RRb::rrLoop(double V_j_min, double V_j_mid, double V_j_max, int J) {


	return;
}
*/
