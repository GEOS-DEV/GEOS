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
 * @file ConformingVirtualElement_1.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
#define GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_

#include "VirtualElementBase.hpp"

namespace geosx
{
  namespace virtualElement
  {
    class ConformingVirtualElementOrder1 final : public VirtualElementBase
    {
      protected:

      localIndex numQuadraturePoints;
      localIndex numSupportPoints;

      public:
      ConformingVirtualElementOrder1() {}
      ~ConformingVirtualElementOrder1() {}

      void ComputeProjectors( MeshLevel const & mesh,
                              localIndex const & regionIndex,
                              localIndex const & subRegionIndex,
                              localIndex const & cellIndex) override;

      void ComputeFaceProjectors( MeshLevel const & mesh,
                                  localIndex const & faceId,
                                  array1d<real64> & basisProjections ) override;

      localIndex getNumQuadraturePoints() const override { return numQuadraturePoints; }
      localIndex getNumSupportPoints() const override { return numSupportPoints; }
    };
  }
}

#endif // GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
