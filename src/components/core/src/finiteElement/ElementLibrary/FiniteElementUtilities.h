/*
 * FiniteElementUtilities.h
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#ifndef FINITEELEMENTUTILITIES_H_
#define FINITEELEMENTUTILITIES_H_

#include <cassert>
#include "Common/DataTypes.hpp"
//#include "ElementLibrary/Basis.h"

namespace FiniteElementUtilities
{
void Integrate( const R2SymTensor& fieldvar,
                const R1Tensor* const dNdX,
                const realT& detJ,
                const realT& detF,
                const R2Tensor& Finv,
                const int numPoints,
                R1Tensor* const result);


//  void Interp(const R1Tensor &globalCoord,
//	      const array<R1Tensor> &nodeCoords,
//              const array<real64> &nodeValues,
//              BasisBase *basis,
//	      const unsigned int &ndofs,
//	      const unsigned int &ndim,
//              realT &result);
//
//
//  void InterpdNdX(const R1Tensor &globalCoord,
//		  const array<R1Tensor> &nodeCoords,
//		  BasisBase *basis,
//		  const unsigned int &ndofs,
//                  const unsigned int &ndim,
//		  array<R1Tensor> &result);
//
//
//  void FindLocalCoord(const R1Tensor &globalCoord,
//		      const array<R1Tensor> &nodeCoords,
//		      BasisBase *basis,
//		      const unsigned int &ndofs,
//		      const unsigned int &ndim,
//		      R1Tensor &localCoord);

}

#endif /* FINITEELEMENTUTILITIES_H_ */
