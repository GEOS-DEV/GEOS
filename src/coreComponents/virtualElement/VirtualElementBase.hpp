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
 * @file VirtualElementBase.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
#define GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_

#include "mesh/MeshLevel.hpp"

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

      /// @brief Compute VEM projectors on the geometry.
      virtual void ComputeProjectors( MeshLevel const &, localIndex const &,
                                      localIndex const &, localIndex const &) = 0;

      /**
       * @brief Virtual getter for the number of quadrature points per element.
       * @return The number of quadrature points per element.
       */
      virtual localIndex getNumQuadraturePoints() const = 0;

      /**
       * @brief Virtual getter for the number of support points per element.
       * @return The number of support points per element.
       */
      virtual localIndex getNumSupportPoints() const = 0;
    };
  }
}

#endif // GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
