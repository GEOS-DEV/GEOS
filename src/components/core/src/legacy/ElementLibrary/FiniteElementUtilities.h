/*
 * FiniteElementUtilities.h
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#ifndef FINITEELEMENTUTILITIES_H_
#define FINITEELEMENTUTILITIES_H_

#include <cassert>
#include "Common/Common.h"
#include "ElementLibrary/Basis.h"

namespace FiniteElementUtilities
{
  void Integrate( const R2SymTensor& fieldvar,
                  const R1Tensor* const dNdX,
                  const realT& detJ,
                  const realT& detF,
                  const R2Tensor& Finv,
                  const int numPoints,
                  R1Tensor* const result);


  void Interp(const R1Tensor &globalCoord,
	      const Array1dT<R1Tensor> &nodeCoords,
              const rArray1d &nodeValues,
              Basis *basis, 
	      const unsigned int &ndofs,
	      const unsigned int &ndim,
              realT &result);


  void InterpdNdX(const R1Tensor &globalCoord,
		  const Array1dT<R1Tensor> &nodeCoords,
		  Basis *basis,
		  const unsigned int &ndofs,
                  const unsigned int &ndim, 
		  Array1dT<R1Tensor> &result);


  void FindLocalCoord(const R1Tensor &globalCoord,
		      const Array1dT<R1Tensor> &nodeCoords,
		      Basis *basis,
		      const unsigned int &ndofs,
		      const unsigned int &ndim, 
		      R1Tensor &localCoord);

}

#endif /* FINITEELEMENTUTILITIES_H_ */
