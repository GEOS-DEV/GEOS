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
 * @file LinearTriangleShapeFunctionKernel.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_LINEARTRIANGLEFACE
#define GEOSX_CORE_FINITEELEMENT_LINEARTRIANGLEFACE

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{

class LinearTriangleFaceShapeFunctionKernel : public FiniteElementBase
{
public:

  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 3;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 1;

  constexpr static real64 parentArea = 0.5;
  constexpr static real64 weight = parentArea / numQuadraturePoints;
  constexpr static real64 quadratureFactor = 1.0 / 3.0;


  virtual ~LinearTriangleFaceShapeFunctionKernel() override final
  {}

  virtual localIndex getNumQuadraturePoints() const override final
  {
    return numQuadraturePoints;
  }

  virtual localIndex getNumSupportPoints() const override final
  {
    return numNodes;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void shapeFunctionValues( localIndex const q,
                                   real64 (& N)[numNodes] )
  {
    GEOSX_UNUSED_VAR( q );
    N[0] = N[1] = N[2] = 1.0 / 3.0;
  }

  GEOSX_HOST_DEVICE
  static real64 JxW( localIndex const q,
                     real64 const (&X)[numNodes][3] )
  {
    GEOSX_UNUSED_VAR( q );
    real64 n[3] = { ( X[1][1] - X[0][1] ) * ( X[2][2] - X[0][2] ) - ( X[2][1] - X[0][1] ) * ( X[1][2] - X[0][2] ),
                    ( X[2][0] - X[0][0] ) * ( X[1][2] - X[0][2] ) - ( X[1][0] - X[0][0] ) * ( X[2][2] - X[0][2] ),
                    ( X[1][0] - X[0][0] ) * ( X[2][1] - X[0][1] ) - ( X[2][0] - X[0][0] ) * ( X[1][1] - X[0][1] )};
    return sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] ) * weight;
  }

private:

};

}
}
#endif //GEOSX_CORE_FINITEELEMENT_LINEARTRIANGLEFACE
