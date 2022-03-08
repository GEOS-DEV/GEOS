#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "misc.h"

using namespace std;

Combinations::Combinations(int ni, int ne, int l) {
	// https://www.geeksforgeeks.org/print-all-possible-combinations-of-r-elements-in-a-given-array-of-size-n/
	Ne = ne; L = l;
	m_arr = std::vector<int>(Ne);
	for (int i = 0; i < Ne; i++) 
	{ 
		m_arr[i] = i;   // array with component numbers for each m_i(v)
	}
	m_data = std::vector<int>(ni*L);
	std::vector<int> temp(L, 0);
	unique_combinations(temp, 0, Ne-1, 0);
}

void Combinations::unique_combinations(std::vector<int> temp, int start, int end, int index) {
	/* std::vector<int> arr ---> Input Array
	std::vector<int> data ---> Temporary array to store current combination
	start & end ---> Starting and ending indexes in arr
	index ---> Current index in data
	r ---> Size of a combination to be printed */
	// This code is contributed by rathbhupendra
	if (index == L)	
	{ 
		for (int j = 0; j < L; j++) 
		{ 
			m_data[index0 + j] = temp[j]; 
		}
		index0 += L;
		return; 
	}

	// replace index with all possible elements. The condition "end-i+1 >= r-index" makes sure
	// that including one element at index will make a combination with remaining elements at remaining positions
	for (int i = start; i <= end && end - i + 1 >= L - index; i++)
	{
		temp[index] = m_arr[i];
		unique_combinations(temp, i+1, end, index+1);
	}
}

int factorial(int n) {
    if (n > 1) 
	{ 
		return n * factorial(n - 1); 
	}
    else 
	{ 
		return 1; 
	}
}

int searchString(std::vector<std::string> A, int size, std::string target)
{
	for (int j = 0; j < size; j++) 
	{
		if (A[j] == target) 
		{
			return j;
		}
	}
	cout << "Component or phase not found!" << endl;
	return 100;
}