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
 * @file VolumetricInterpolant.h
 *
 *  Created on: September 4, 2013
 *  Author: scottjohnson
 */

#ifndef VOLUMEKERNEL_H_
#define VOLUMEKERNEL_H_

#include "../TensorT.old/R1TensorT.h"
#include "Interpolation/InterpolantBaseT.h"

/**
 * @author Scott Johnson
 * @brief VolumetricInterpolant holds the definition of a 3D interpolation
 * kernel with compact support
 */
class VolumeKernel : public InterpolantBaseT
{
public:
  VolumeKernel(): m_x(0.0)
  {};

  void Initialize(const R1Tensor& x, const realT h, const realT mean, const realT stdev);

  virtual realT SetInverseSmoothingLengthFactor(const realT h);

  void Evaluate(const R1Tensor& x, realT& fmW4, realT& mW4) const;

  inline realT mW(const realT dxd2) const
  {
    return m_value * W(dxd2);
  }

private:
  inline realT W(const realT dxd2) const
  {
    return M4(dxd2 * m_invh2, Alpha()) * m_invhd;
  }

  static inline realT Alpha()
  {
    static realT alpha = InterpolantBaseT::Alpha(3);
    return alpha;
  }

public:
  R1Tensor m_x;

};

#endif /* VOLUMEKERNEL_H_ */
