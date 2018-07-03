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
 * @file FractalVolume.h
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#ifndef FRACTALVOLUME_H_
#define FRACTALVOLUME_H_

#include "../TensorT.old/R1TensorT.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"
#include "StatisticalDistributionBaseT.h"
#include "Interpolation/VolumeKernel.h"
#include "FractalBaseT.h"

/**
 * @author Scott Johnson
 * @brief FractalVolume creates self-similar distributions of apertures
 */
class FractalVolume : public FractalBaseT
{
public:
  FractalVolume();

  unsigned InitializeFractal(const realT mean, const realT stdev,
                             const R1Tensor& lower, const R1Tensor& upper,
                             const realT hurst = 1.3, const localIndex nlevels = 6,
                             const localIndex n0 = 1, const localIndex n1 = 1,
                             const localIndex n2 = 1, const realT hfct = 1.2,
                             const unsigned seed = 0);

  unsigned Initialize(const Array2dT<realT>& parameters,
                      const R1Tensor& lower, const R1Tensor& upper,
                      const localIndex n0 = 1, const localIndex n1 = 1,
                      const localIndex n2 = 1, const realT hfct = 1.2,
                      const unsigned seed = 0);

  realT Value(const R1Tensor& position) const;

  localIndex Positions(const realT dx,
                       array<R1Tensor>& positions,
                       const realT weight = 1.0) const;

  localIndex Positions(const R1Tensor& min,
                       const R1Tensor& max,
                       const realT dx,
                       array<R1Tensor>& positions,
                       const realT weight = 1.0) const;

private:
  void InitializeSumMW(const int nj, const int nk, const int nl, const int ioffset,
                       const localIndex ilevel, Array3dT<VolumeKernel*>& curr);

  void InitializeLevel(const int nj, const int nk, const int nl, const int ioffset,
                       const realT dx0, const realT dx1, const realT dx2,
                       const realT mean, const realT stdev, const realT h,
                       Array3dT<VolumeKernel*>& curr, localIndex& icurr);

  void FillValues(const int ioffset, const localIndex nlevels,
                  array< Array3dT<VolumeKernel*>* >& vals);

  ///lower 3D coordinates of the points that will be queried
  R1Tensor m_lower;

  ///upper 3D coordinates of the points that will be queried
  R1Tensor m_upper;

  ///Initial value of n2
  localIndex m_n2;

  ///For each cell at the finest level, holds, for each level, the kernels
  // applicable to the cell
  Array3dT<array<array<VolumeKernel*> > > m_values;

  ///Kernels referenced in values
  array<VolumeKernel> m_kernels;
};

#endif /* FRACTALVOLUME_H_ */
