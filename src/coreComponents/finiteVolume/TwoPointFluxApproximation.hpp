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
 * @file TwoPointFluxApproximation.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geosx
{

class TwoPointFluxApproximation : public FluxApproximationBase
{
public:

  static std::string CatalogName() { return "TwoPointFluxApproximation"; }

  TwoPointFluxApproximation() = delete;

  TwoPointFluxApproximation(std::string const & name, dataRepository::Group * const parent);

protected:

  virtual void computeCellStencil( DomainPartition const & domain ) override;

  virtual void addToFractureStencil( DomainPartition const & domain,
                                     string const & faceElementRegionName,
                                     bool const initFlag ) override;

  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       SortedArray<localIndex> const & faceSet,
                                       BoundaryStencil & stencil ) override;

};

}


#endif //GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
