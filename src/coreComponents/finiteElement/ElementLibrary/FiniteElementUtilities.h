/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file FiniteElementUtilities.h
 */

#ifndef FINITEELEMENTUTILITIES_H_
#define FINITEELEMENTUTILITIES_H_

#include <cassert>
#include "Common/DataTypes.hpp"
//#include "ElementLibrary/Basis.h"

namespace FiniteElementUtilities
{
void Integrate( const R2SymTensor & fieldvar,
                const R1Tensor * const dNdX,
                const realT & detJ,
                const realT & detF,
                const R2Tensor & Finv,
                const int numPoints,
                R1Tensor * const result );


//  void Interp(const R1Tensor &globalCoord,
//	      const array1d<R1Tensor> &nodeCoords,
//              const array1d<real64> &nodeValues,
//              BasisBase *basis,
//	      const unsigned int &ndofs,
//	      const unsigned int &ndim,
//              realT &result);
//
//
//  void InterpdNdX(const R1Tensor &globalCoord,
//		  const array1d<R1Tensor> &nodeCoords,
//		  BasisBase *basis,
//		  const unsigned int &ndofs,
//                  const unsigned int &ndim,
//		  array1d<R1Tensor> &result);
//
//
//  void FindLocalCoord(const R1Tensor &globalCoord,
//		      const array1d<R1Tensor> &nodeCoords,
//		      BasisBase *basis,
//		      const unsigned int &ndofs,
//		      const unsigned int &ndim,
//		      R1Tensor &localCoord);

}

#endif /* FINITEELEMENTUTILITIES_H_ */
