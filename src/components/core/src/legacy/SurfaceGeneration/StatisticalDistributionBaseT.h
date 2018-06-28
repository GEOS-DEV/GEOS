/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * StatisticalDistributionBaseT.h
 *
 *  Created on: Jun 7, 2012
 *      Author: johnson346
 */

#ifndef STATISTICALDISTRIBUTIONBASET_H_
#define STATISTICALDISTRIBUTIONBASET_H_

#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"
#include "Utilities/Utilities.h"
#include <cmath>
#include "math/TensorT/R1TensorT.h"

class StatisticalDistributionBaseT
{
public:
  StatisticalDistributionBaseT();
  virtual ~StatisticalDistributionBaseT();

  /**
   * @brief Randomly initializes the RNG seed
   * @author Scott Johnson
   */
  inline static unsigned InitializeRandom(const unsigned seed = 0)
  {
    const unsigned seedReturn = seed > 0 ? seed : (unsigned)time(0);
    srand(seedReturn);
    return seedReturn;
  };

  enum StatisticalDistributionParameters
  {
    STANDARD_DEVIATION = 0,
    MEAN = 1,
    HURST_EXPONENT = 2,
    MINIMUM_VALUE = 3,
    MAXIMUM_VALUE = 4,
    WEIBULL_SHAPE = 5,
    WEIBULL_SCALE = 6,
    POWERLAW_EXPONENT = 7,
    N_PARAMS = 8
  };

//  ///From "Handbook of Mathematical Functions" fml 7.1.26 (O(1e-7) error)
//  static inline realT ErrorFunctionA(const realT xx)
//  {
//    const realT sign = xx >= 0 ? 1.0 : -1.0;
//    const realT x = sign * xx;
//
//    // constants
//    const realT a1 = 0.254829592;
//    const realT a2 = -0.284496736;
//    const realT a3 = 1.421413741;
//    const realT a4 = -1.453152027;
//    const realT a5 = 1.061405429;
//    const realT p = 0.3275911;
//
//    //A&S formula 7.1.26
//    const realT t = 1.0 / (1.0 + p * x);
//    const realT y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t
// * exp(-x * x);
//    return sign * y; //std::erf(-x) = -std::erf(x)
//  }

  ///From "Numerical Recipes in Fortran" 1992, p. 214 (O(1e-7) error)
  static inline realT ErrorFunction(const realT xx)
  {
//    const realT sign = xx >= 0 ? 1.0 : -1.0;
//    const realT x = sign * xx;
//    const realT t = 1.0 / (1 + 0.5 * x);
//
//    // constants
//    const realT a0 = -1.26551223;
//    const realT a1 = 1.00002368;
//    const realT a2 = 0.37409196;
//    const realT a3 = 0.09678418;
//    const realT a4 = -0.18628806;
//    const realT a5 = 0.27886807;
//    const realT a6 = -1.13520398;
//    const realT a7 = 1.48851587;
//    const realT a8 = -0.82215223;
//    const realT a9 = 0.17087277;
//
//    const realT t2 = t * t;
//    const realT t4 = t2 * t2;
//    const realT t6 = t2 * t4;
//    const realT t8 = t4 * t4;
//    const realT tau = t * exp(-xx * xx + a0 + a1 * t + a2 * t2 +
//                              a3 * t * t2 + a4 * t4 + a5 * t4 * t +
//                              a6 * t6 + a7 * t6 * t + a8 * t8 + a9 * t8 * t);
//    return tau;

    //NOW THAT WE'RE USING C++11, WE CAN USE THE STD IMPLEMENTATION
    return erf(xx);
  }

  //Cumulative distribution function of the normalized Gaussian
  static inline realT PhiNormal(const realT z)
  {
    return 0.5 * (1.0 + erf(z / sqrt(2.0)));
  }

  //Find the probability that a normal variate is on the interval (z0, z1)
  static inline realT PhiNormalDifference(const realT z0, const realT z1)
  {
    return 0.5 * (erf(z1 / sqrt(2)) - erf(z0 / sqrt(2)));
  }

  //Find the probability that a normal variate is on the interval (z0, inf)
  static inline realT PhiNormalDifferenceInf(const realT z0)
  {
    //return 0.5 * (1 - std::erf(z0 / sqrt(2)));
    return 0.5 * erfc(z0 / sqrt(2));
  }

  static realT ErrorFunctionInverse(const realT erf_value,
                                    const realT zlower = -1e6,
                                    const realT zupper = 1e6,
                                    const realT tol = 1.0e-6);



  //Find the probability that a normal variate is on the interval
  static inline realT ProbabilityOfExceedence(const realT z)
  {
    //const realT erfc = 1.0 - ErrorFunction( z / sqrt(2.0));

    //find the probability that a normal variate is on the interval z <= x < inf
    const realT pz = 0.5*erfc(z / sqrt(2.0));

    //find the probability that a normal variate is on the interval -inf < x < z
    const realT probabilityOfExceedence = 1.0-pz;

    return probabilityOfExceedence;
  }


  /**
   * @brief Provides a normally distributed random sample
   * @author Scott Johnson
   * @param mean Mean of the normal distribution
   * @param stdev Standard deviation of teh normal distribution
   * @return Random sample
   */
  static inline realT NormalSample(const realT mean, const realT stdev)
  {
    realT z = BoxMuller();
    realT x = ScaleSample (z, mean, stdev);
    return x;
  };

  inline realT NormalSample() const
  {
    const realT* mean = GetParameter(MEAN);
    if(!mean)
      throw GPException("Cannot find mean for a normal sample");
    const realT* stdev = GetParameter(STANDARD_DEVIATION);
    if(!stdev)
      throw GPException("Cannot find mean for a normal sample");
    const realT* min = GetParameter(MINIMUM_VALUE);
    const realT* max = GetParameter(MAXIMUM_VALUE);

    realT ret = NormalSample(*mean, *stdev);
    if(min)
      ret = ret > *min ? ret : *min;
    if(max)
      ret = ret < *max ? ret : *max;
    return ret;
  };

  static inline realT UniformSample(const realT min, const realT max)
  {
    const float fct = (realT)rand()/(realT)RAND_MAX;
    const realT ret = min + (max - min) * fct;
    return ret;
  };

  static inline realT FrequencyToNumber(const realT frequency)
  {
    //floor operation on frequency + probability of one more given by the
    // modulus
    return ((int) frequency) + (UniformSample(0.0, 1.0) < fmod(frequency, 1.0) ? 1 : 0);
  };

  inline realT UniformSample() const
  {
    const realT* min = GetParameter(MINIMUM_VALUE);
    const realT* max = GetParameter(MAXIMUM_VALUE);
    if(!min)
      throw GPException("Cannot find minimum for a uniform sample");
    if(!min)
      throw GPException("Cannot find maximum for a uniform sample");

    return UniformSample(*min, *max);
  };

  inline realT WeibullSample(const realT tol = 1e-30) const
  {

    const realT* shape = GetParameter(WEIBULL_SHAPE);//k
    if(!shape || (isZero(*shape)))
      throw GPException("Cannot sample Weibull distribution without defining the shape parameter!");

    const realT* scaling = GetParameter(WEIBULL_SCALE);//lambda
    if(!scaling || (isZero(*scaling)))
      throw GPException("Cannot sample Weibull distribution without defining the scaling parameter!");

    const realT ishape = 1.0 / (*shape);
    const realT ialpha = pow((*scaling), (*shape));

    //    realT gamma;
    //    {
    //      realT* gammaPtr = GetParameter(WEIBULL_GAMMA);
    //      if(!gammaPtr)
    //      {
    //        gamma = Gamma(ishape + 1.0);
    //        AddParameter(WEIBULL_GAMMA, gamma);
    //      }
    //      else
    //      {
    //        gamma = *gammaPtr;
    //      }
    //    }

    const realT x = UniformSample(std::numeric_limits<realT>::min() > tol ? std::numeric_limits<realT>::min() : tol, 1.0);
    const realT ret = pow(-ialpha * log(x), ishape);
    return ret;
  }

  inline realT PowerLawSample() const
  {
    const realT* exponent = GetParameter(POWERLAW_EXPONENT);
    if (!exponent)
      throw GPException("Cannot find exponent for a power law sample");
    const realT* min = GetParameter(MINIMUM_VALUE);
    if (!min)
      throw GPException("Cannot find minimum for a power law sample");
    const realT* max = GetParameter(MAXIMUM_VALUE);
    if (!max)
      throw GPException("Cannot find maximum for a power law sample");

    const realT Np1 = (*exponent) + 1;
    if(isZero(Np1))
      throw GPException("Cannot transform a power law with exponent -1");
    const realT xmin_np1 = pow(*min, Np1);
    const realT xmax_np1 = pow(*max, Np1);
    const realT y = UniformSample(0.0, 1.0);
    const realT ret = pow((xmax_np1 - xmin_np1) * y + xmin_np1, 1.0 / Np1);
    return ret;
  }

  inline void AddParameter(const StatisticalDistributionParameters type, const realT value) { parameters[(int)type] = value; }

  inline const realT* GetParameter(const StatisticalDistributionParameters type) const
  {
    return parameters[(int)type] < std::numeric_limits<realT>::max() ?
           &parameters[(int)type] : 0;
    //return stlMapLookupPointer( this->parameters, (int)type);
  }

  inline localIndex NumberOfParameters() const
  {
    int n = 0;
    for(int i = 0 ; i < N_PARAMS ; i++, n+= ((parameters[i] < std::numeric_limits<realT>::max()) ? 1 : 0)) {}
    return n;
    //return this->parameters.size();
  }

  void AddWeibullParameters(const realT mean, const realT stdev);

private:

  //std::map<int, realT> parameters;
  realT parameters[N_PARAMS];

  static realT
  FindWeibullShape(const realT* const params,
                   const realT lower,
                   const realT upper,
                   const realT tol);

  static inline realT
  WeibullShapeFunction(const realT shape, const realT* const params)
  {
    //note: for Weibull, the scaling and shape parameters are related to the
    // mean and standard deviation by the following:
    //mean = scaling * gamma(1+1/shape)
    //scaling = mean / gamma(1+1/shape)
    //stdev^2 = scaling * scaling * gamma(1+2/shape) - mean * mean
    //Therefore,
    //We need to iterate to evaluate gamma, so find the zero of ...
    //mean * mean * gamma(1+2/shape) / (gamma(1+1/shape) * gamma(1+1/shape)) -
    // mean*mean - stdev*stdev = 0
    //what should we use for an initial value of "shape"?

    const realT mean = params[0];
    const realT stdev2 = params[1];
    const realT mean2 = mean*mean;

    const realT mean_div_g1 = mean / Gamma(1.0+1.0/shape);
    const realT g2 = Gamma(1.0+2.0/shape);

    return mean_div_g1*mean_div_g1*g2 - mean2 - stdev2;
  }

  static inline realT ErfZero(const realT x, const realT* params)
  {
    return erf(x) - params[0];
  }

  /**
   * @brief Scales an unscaled, normally distributed sample to the given mean
   * and stdev
   * @author Scott Johnson
   * @param z Unscaled sample
   * @param mean Mean of the normal distribution
   * @param stdev Standard deviation of the normal distribution
   * @return Scaled sample value
   */
  static inline realT ScaleSample(const realT z, const realT mean, const realT stdev)
  {
    realT x = stdev * z + mean;
    return x;
  };

  /**
   * @brief Generates a normally distributed sample with mean 0 and stdev 1
   * @author Scott Johnson
   * Generates a normally distributed sample using the Box-Muller algorithm
   * @return Normally distributed random sample (z = 1)
   */
  static inline realT BoxMuller()
  {
    float u = 0;
    float v = 0;

    while (u == 0)
    {
      u = (float)rand()/(float)RAND_MAX;
      v = (float)rand()/(float)RAND_MAX;
    }

    u = log(u);
    u *= -2.0;
    u = sqrt(u);
    v *= 2.0 * 3.141592653589793238462;
    v = cos(v);
    u *= v;
    return u;
  };

  ///Calculate the gamma function
  /**
     The gamma function is defined for integers
     as \f$\Gamma(n+1)=n!\f$; however, for non-integer
     numbers, \f$\Gamma(z)\f$ is often approximated using
     -Stirling's approximation
     -Lanczos' approximation
     .
     Here, the approximation used is from a 3rd party class:
     gamma.cpp -- computation of gamma function.
        Algorithms and coefficient values from "Computation of Special
        Functions", Zhang and Jin, John Wiley and Sons, 1996.
        (C) 2003, C. Bond. All rights reserved.

     If the input is an integer, the problem is well-defined,
     and an exact solution is given, though, overflow occurs at
     values < 0. Also, overflow occurs at values higher than 171
     for both real and integer values. Overflow is defined by a
     return value of 1e308.

     \param x Value being evaluated by the gamma operator
     \return Approximation of the gamma function
   */
  static inline realT Gamma(const realT x)
  {
    int i = 0, k = 0, m = 0;
    double ga = 0, gr = 0, r = 0, z = 0;

    const double g[] =
    { 1.0, 0.5772156649015329, -0.6558780715202538, -0.420026350340952e-1, 0.1665386113822915,
      -0.421977345555443e-1, -0.9621971527877e-2, 0.7218943246663e-2, -0.11651675918591e-2,
      -0.2152416741149e-3, 0.1280502823882e-3, -0.201348547807e-4, -0.12504934821e-5,
      0.1133027232e-5, -0.2056338417e-6, 0.6116095e-8, 0.50020075e-8, -0.11812746e-8, 0.1043427e-9,
      0.77823e-11, -0.36968e-11, 0.51e-12, -0.206e-13, -0.54e-14, 0.14e-14 };

    if (x > 171.0)
      return 1e308; // This value is an overflow flag.
    if (isEqual(x, (int) x))
    {
      if (x > 0.0)
      {
        ga = 1.0; // use factorial
        for (i = 2 ; i < x ; i++)
        {
          ga *= i;
        }
      }
      else
        ga = 1e308;
    }
    else
    {
      if (fabs(x) > 1.0)
      {
        z = fabs(x);
        m = (int) z;
        r = 1.0;
        for (k = 1 ; k <= m ; k++)
        {
          r *= (z - k);
        }
        z -= m;
      }
      else
        z = x;
      gr = g[24];
      for (k = 23 ; k >= 0 ; k--)
      {
        gr = gr * z + g[k];
      }
      ga = 1.0 / (gr * z);
      if (fabs(x) > 1.0)
      {
        ga *= r;
        if (x < 0.0)
        {
          ga = -M_PI / (x * ga * sin(M_PI * x));
        }
      }
    }
    return ga;
  }
};

#endif /* STATISTICALDISTRIBUTIONBASET_H_ */
