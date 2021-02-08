/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VirtualElementBase.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
#define GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_

#include "mesh/CellBlock.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/NodeManager.hpp"

namespace geosx
{
namespace virtualElement
{
/**
 * @brief Base class for VEM implementations.
 */
class VirtualElementBase
{
public:

  /// Default constructor.
  VirtualElementBase() = default;

  /// Default destructor.
  virtual ~VirtualElementBase() = default;

  // /**
  //  * @brief Virtual getter for the number of quadrature points per element.
  //  * @return The number of quadrature points per element.
  //  */
  // virtual localIndex getNumQuadraturePoints() const = 0;

  // /**
  //  * @brief Virtual getter for the number of support points per element.
  //  * @return The number of support points per element.
  //  */
  // virtual localIndex getNumSupportPoints() const = 0;

  /**
   * @brief Virtual getter for a value of the stabilization matrix.
   * @param iBasisFunction The row index.
   * @param jBasisFunction The column index.
   * @return The requested value.
   */
  // virtual real64 getStabilizationValue( localIndex const iBasisFunction,
  //                                       localIndex const jBasisFunction
  //                                       ) const = 0;

  /**
   * @brief Get the shape function projected derivatives at a given quadrature point.
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projected derivatives. Size will be @ref
   * getNumSupportPoints() x 3.
   * @note See documentation of derived classes to know the projection used.
   */
  // virtual void getGradN( localIndex const q, arrayView2d< real64 const > & gradN ) const = 0;

  /**
   * @brief Get the shape function projections at a given quadrature point.
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projections. Size will be @ref getNumSupportPoints().
   * @note See documentation of derived classes to know the projection used.
   */
  // virtual void getN( localIndex const q, arrayView1d< real64 const > & N ) const = 0;

  /**
   * @brief Get the integration weight for a quadrature point.
   * @param q Index of the quadrature point.
   * @return The weight.
   */
  // virtual real64 transformedQuadratureWeight( localIndex const q ) const = 0;
};
}
}

#endif // GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
