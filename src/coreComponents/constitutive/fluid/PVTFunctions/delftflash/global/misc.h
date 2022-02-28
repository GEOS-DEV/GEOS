#pragma once
#include <complex>
#include <vector>
#include <iostream>

class Combinations {
private:
    std::vector<int> m_arr;
    std::vector<int> m_data;
    int Ne; // number of unique elements
    int L; // length of combinations
    int index0 = 0;
public:
    Combinations(int ni, int ne, int l);
    int getIndex(int i, int k) { return m_data[i*L + k]; }
private:
    void unique_combinations(std::vector<int> temp, int start, int end, int index);
};

std::vector<std::complex<double>> cubic_roots(double a, double b, double c, double d);
int searchString(std::vector<std::string> A, int size, std::string target);

int factorial(int n);