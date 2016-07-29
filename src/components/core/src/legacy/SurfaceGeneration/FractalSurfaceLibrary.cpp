/*
 * FractalSurfaceLib.cpp
 *
 *  Created on: Jun 5, 2012
 *      Author: johnson346
 */

#include "FractalSurfaceLibrary.h"
#include "FractalSurface.h"

static FractalSurface fractalSurface;

int GenerateFractalSurface(double loweru, double lowerv,
			   double upperu, double upperv,
			   double hurst, double hfct,
			   double mean, double stdev,
			   int nlevels, int n0, int n1)
{
  //u0, v0, u1, v1, 1.3, 1.2, mu, sig, 5, 5, 5

  //Initialize
  R1TensorT<2> lower;
  lower[0] = loweru;
  lower[0] = lowerv;

  R1TensorT<2> upper;
  upper[0] = upperu;
  upper[1] = upperv;

  fractalSurface.InitializeFractal(mean, stdev,
				lower, upper,
				hurst, nlevels,
				n0, n1, hfct);
  return 0;
}

double ValueFractalSurface(double u, double v)
{
  R1TensorT<2> point;
  point(0) = u;
  point(1) = v;
  return fractalSurface.Value(point);
}
