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
 * @file FiniteElementUtilities.cpp
 */

#include "FiniteElementUtilities.h"
#include "Utilities/TrilinosUtilities.h"

namespace FiniteElementUtilities
{
void Integrate( const R2SymTensor & fieldvar,
                const R1Tensor * const dNdX,
                const realT & detJ,
                const realT & detF,
                const R2Tensor & Finv,
                const int numPoints,
                R1Tensor * const result )
{
  const realT integrationFactor = detJ * detF;

  R2Tensor P;
  P.AijBkj( fieldvar, Finv );
  P *= integrationFactor;

  for( int a=0; a<numPoints; ++a )      // loop through all shape functions in
                                        // element
  {
    result[a].minusAijBj( P, dNdX[a] );
  }


}

//
//  void Interp(const R1Tensor &globalCoord,
//	      const array1d<R1Tensor> &nodeCoords,
//              const array1d<real64> &nodeValues,
//              BasisBase *basis,
//	      const unsigned int &ndofs,
//	      const unsigned int &ndim,
//              realT &result)
//  {
//
//    assert(nodeValues.size() == ndofs);
//    assert(nodeCoords.size() == ndofs);
//
//    R1Tensor localCoord;
//
//    FindLocalCoord(globalCoord, nodeCoords, basis, ndofs, ndim, localCoord);
//
//    result = 0.0;
//
//    for(unsigned int i =0; i < ndofs; i++)
//    {
//
//      result += nodeValues[i] * basis->value(i, localCoord);
//
//    }
//
//  }
//
//  // Calcuate dN/dX at a specific position
//
//  void InterpdNdX(const R1Tensor &globalCoord,
//		  const array1d<R1Tensor> &nodeCoords,
//		  BasisBase *basis,
//		  const unsigned int &ndofs,
//                  const unsigned int &ndim,
//		  array1d<R1Tensor> &result)
//  {
//
//    assert(nodeCoords.size() == ndofs);
//
//    R1Tensor localCoord;
//
//    FindLocalCoord(globalCoord, nodeCoords, basis, ndofs, ndim, localCoord);
//
//    R2TensorT<3> jacobian;
//    R2TensorT<3> inv_jacobian;
//
//    result.resize(ndofs);
//
//    jacobian = 0;
//
//    for(unsigned int i = 0; i < ndofs; i++)
//    {
//
//      jacobian.plus_dyadic_ab(nodeCoords[i], basis->gradient(i, localCoord));
//
//    }
//
//    if(ndim == 2)
//    {
//      jacobian(2, 2) = 1;
//
//    }
//
//    inv_jacobian.Inverse(jacobian);
//
//    for(unsigned int i = 0; i < ndofs; i++)
//    {
//      result[i].AijBi(inv_jacobian, basis->gradient(i, localCoord));
//    }
//
//}
//
//
//  void FindLocalCoord(const R1Tensor &globalCoord,
//		      const array1d<R1Tensor> &nodeCoords,
//		      BasisBase *basis,
//		      const unsigned int &ndofs,
//		      const unsigned int &ndim,
//		      R1Tensor &localCoord)
//  {
//
//    assert(nodeCoords.size() == ndofs);
//    array1d<real64> resid(ndim);
//    rArray2d A(ndim, ndim);
//    array1d<real64> x(ndim);
//
//    localCoord = 0.0;
//
//    array1d<real64> values(ndofs);
//    array1d<R1Tensor>  gradients(ndofs);
//
//    realT nrTol = 1.0e-6;
//    int iterMax = 4;
//
//    realT maxV;
//
//    bool converged = 0;
//
//    int iter = 0;
//
//    while(iter < iterMax) {
//
//      for(unsigned int i = 0; i < ndofs; i++) {
//
//	values[i] = basis->value(i, localCoord);
//	gradients[i] = basis->gradient(i, localCoord);
//
//      }
//
//      for(unsigned int i = 0; i < ndim; i++) {
//
//	for(unsigned int j = 0; j < ndim; j++) {
//
//	  A[i][j] = 0.0;
//
//	  for(unsigned int k = 0; k < ndofs; k++) {
//
//	    A[i][j] += nodeCoords[k][i] * gradients[k][j];
//
//	  }
//
//	}
//
//	resid[i] = -globalCoord[i];
//
//	for(unsigned int j = 0; j < ndofs; j++) {
//
//	  resid[i] += nodeCoords[j][i] * values[j];
//
//	}
//
//      }
//
//      maxV = -1e30;
//
//      for(unsigned int i = 0; i < ndim; i++) {
//
//	if(fabs(resid[i]) > maxV)
//	  maxV = fabs(resid[i]);
//
//      }
//
//      if(maxV < nrTol) {
//	converged = 1;
//	break;
//      }
//
//      LinSolve_Local(A, x, resid);
//
//      for(unsigned int i = 0; i < ndim; i++) {
//
//	localCoord[i] -= x[i];
//
//      }
//
//      iter++;
//
//    }
//
//    //sanity check
//
//    for(unsigned int i = 0; i < ndim; i++) {
//
//      if(localCoord[i] < 0.0 || localCoord[i] > 1.0) {
//	converged = 0;
//	break;
//      }
//
//    }
//
//    if(!converged) {
//      throw GPException("ERROR: Failed to find the local coord");
//    }
//
//  }

}
