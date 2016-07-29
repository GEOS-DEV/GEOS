/*
 * StatisticalDistributionBaseT.cpp
 *
 *  Created on: Jun 7, 2012
 *      Author: johnson346
 */

#include "StatisticalDistributionBaseT.h"
#include "Utilities/FindRoots.h"

StatisticalDistributionBaseT::StatisticalDistributionBaseT()//: parameters()
{
  for(int i = 0; i < N_PARAMS; i++)
  {
    parameters[i] = std::numeric_limits<realT>::max();
  }
}

StatisticalDistributionBaseT::~StatisticalDistributionBaseT()
{
  // TODO Auto-generated destructor stub
}

void
StatisticalDistributionBaseT::AddWeibullParameters(const realT mean, const realT stdev)
{
  realT shape = -1;

  //find the zeros of the function
  {
    const realT stdev2 = stdev*stdev;
    const realT params[] = {mean, stdev2};
//    realT r = -1.9;
//    for(int i = 0; i < 500; i++, r+=0.01)
//    {
//      const realT val = pow(10.0, r);
//      const realT fval = WeibullShapeFunction(val, params);
//      std::cout << "wfunc " << val << " " << fval << std::endl;
//    }
    shape = FindWeibullShape(params, 1.2e-2, 1e3, 1e-6);
  }

  //scaling = mean / gamma(1+1/shape)
  const realT scaling = mean / Gamma(1.0+1.0/shape);

  AddParameter(WEIBULL_SHAPE, shape);
  AddParameter(WEIBULL_SCALE, scaling);
}

realT
StatisticalDistributionBaseT::FindWeibullShape(const realT* const params,
                                               const realT lower,
                                               const realT upper,
                                               const realT tol)
{
  return FindRoots::FindRoot(WeibullShapeFunction, params, lower, upper, tol);
}

realT
StatisticalDistributionBaseT::ErrorFunctionInverse(const realT erf_value,
                                                   const realT zlower,
                                                   const realT zupper,
                                                   const realT tol)
{
  const realT params[] = {erf_value};
  return FindRoots::FindRoot(ErfZero, params, zlower, zupper, tol);
}


