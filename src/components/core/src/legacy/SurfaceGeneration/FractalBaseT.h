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
 * @file FractalBaseT.h
 *
 *  Created on: June 23, 2012
 *  Author: scottjohnson
 */

#ifndef FRACTALBASET_H_
#define FRACTALBASET_H_

#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "ArrayT/ArrayT.h"

/**
 * @author Scott Johnson
 * @brief FractalBaseT creates self-similar distributions of properties
 */
class FractalBaseT
{
public:
  FractalBaseT();

protected:
  void FillFractalParameters(const realT mean, const realT stdev,
                             const realT hurst, const localIndex nlevels,
                             Array2dT<realT>& parameters);

  ///Smoothing length multiplication factor
  realT m_hfct;

  ///Initial values of n0, n1, and n2
  localIndex m_n0, m_n1;
};

#endif /* FRACTALBASET_H_ */
