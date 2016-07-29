/*
 * FractalBaseT.cpp
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#include "FractalBaseT.h"
#include <limits.h>

/**
 * @brief Default constructor for the aperture generator class
 * @author Scott Johnson
 */
FractalBaseT::FractalBaseT()
{
}

void FractalBaseT::FillFractalParameters(const realT mean, const realT stdev,
                                         const realT hurst, const localIndex nlevels,
                                         Array2dT<realT>& parameters)
{
  const realT fct = pow(2.0,-hurst);

  realT meani = mean;
  realT stdevi = stdev * fct;
  realT sum = 0.0;
  for(localIndex i = 0; i < nlevels; i++, stdevi *= fct, meani *= 0)
  {
    parameters(i, 0) = meani;
    parameters(i, 1) = stdevi;
    sum += stdevi * stdevi;
  }

  //NOTE: do the following to renormalize the standard deviation to yield the originally requested one
  //the sum of N independent normal distributions is mu = sum(mu_i) and sigma^2 = sum(sigma_i^2)
  if(sum > 0)
  {
    const realT fct2 = stdev / sqrt(sum);
    for(localIndex i = 0; i < nlevels; i++)
    {
      parameters(i, 1) *= fct2;
    }
  }
}
