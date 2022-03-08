#pragma once
#include <vector>

class RR {
protected:
    int NP, NC;
    bool admissible;
    std::vector<double> z, K, v_min, v_max;
    std::vector<bool> converged;

public:
    RR(std::vector<double> z_, std::vector<double> K_);
    virtual std::vector<double> solveRR(std::vector<double> K_) = 0;
    std::vector<double> getx(std::vector<double> V);

protected:
    bool bounded();
    std::vector<double> V_limits(std::vector<double> V_j, int j);
    std::vector<double> f_df1(std::vector<double> V_j, int J); // calculates f_j and df_j/dv_j (used only in B and BN)
};


// class RRn :public RR 
// {
// private:
//     Eigen::MatrixXd df;
//     Eigen::VectorXd f, v;
//     std::vector<double> v_min, v_max;
//     double eps = 1E-8;

// public:
//     RRn(std::vector<double> z_, std::vector<double> K_);
//     virtual std::vector<double> solveRR(std::vector<double> K_);
    
// private:
//     void f_df();  // calculate NP-1x1 vector containing f_j and NP-1xNP-1 matrix containing df_i/dv_j
// };


class RRbn :public RR
{
private:
    std::vector<double> f, v;
    double df;
    double eps = 1E-6; // find out why smaller eps doesn't converge for small composition
    int iter = 0;

public:
    RRbn(std::vector<double> z_, std::vector<double> K_);
    virtual std::vector<double> solveRR(std::vector<double> K_);

private:
    void rrLoop(int J);
};


// class RRb :public RR
// {
// private:
//     std::vector<double> v_mid, f_min, f_mid;
//     std::vector<bool> converged;
//     double eps = 1E-5;

// public:
//     RRb(std::vector<double> z_, std::vector<double> K_);
//     virtual std::vector<double> solveRR(std::vector<double> K_);

// private:
//     void rrLoop(double V_j_min, double V_j_mid, double V_j_max, int J);
// };
