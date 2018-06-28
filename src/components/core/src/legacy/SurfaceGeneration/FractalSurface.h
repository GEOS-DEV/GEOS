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

/**
 * @file FractalSurface.h
 *
 *  Created on: Sep 27, 2011
 *  Author: scottjohnson
 */

#ifndef FRACTALSURFACE_H_
#define FRACTALSURFACE_H_

#include "../TensorT.old/R1TensorT.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"
#include "StatisticalDistributionBaseT.h"
#include "Interpolation/SurfaceKernel.h"
#include "FractalBaseT.h"

/**
 * @author Scott Johnson
 * @brief FractalSurface creates self-similar distributions of apertures
 */
class FractalSurface : public FractalBaseT
{
public:
  FractalSurface();

  unsigned InitializeFractal(const realT mean, const realT stdev,
                             const R1TensorT<2>& lower, const R1TensorT<2>& upper,
                             const realT hurst = 1.3, const localIndex nlevels = 6,
                             const localIndex n0 = 1, const localIndex n1 = 1,
                             const realT hfct = 1.2,
                             const unsigned seed = 0);

  unsigned Initialize(const Array2dT<realT>& parameters,
                      const R1TensorT<2>& lower, const R1TensorT<2>& upper,
                      const localIndex n0 = 1, const localIndex n1 = 1,
                      const realT hfct = 1.2,
                      const unsigned seed = 0);

  realT Value(const R1TensorT<2>& position) const;

protected:
  void InitializeSumMW(const int nj, const int nk, const int ioffset,
                       const localIndex ilevel, Array2dT<SurfaceKernel*>& curr);

  void InitializeLevel(const int nj, const int nk, const int ioffset,
                       const realT dx0, const realT dx1,
                       const realT mean, const realT stdev, const realT h,
                       const localIndex ilevel,
                       Array2dT<SurfaceKernel*>& curr, localIndex& icurr);

  void FillValues(const int ioffset, const localIndex nlevels,
                  array< Array2dT<SurfaceKernel*> >& vals);

  ///lower 3D coordinates of the points that will be queried
  R1TensorT<2> m_lower;

  ///upper 3D coordinates of the points that will be queried
  R1TensorT<2> m_upper;

  ///For each cell at the finest level, holds, for each level, the kernels
  // applicable to the cell
  Array2dT<array<array<SurfaceKernel*> > > m_values;

  ///Kernels referenced in values
  array<SurfaceKernel> m_kernels;
};

#endif /* FRACTALSURFACE_H_ */
