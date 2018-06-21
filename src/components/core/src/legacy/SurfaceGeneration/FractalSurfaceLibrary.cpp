// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.
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
