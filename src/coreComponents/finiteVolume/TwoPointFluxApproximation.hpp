/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file TwoPointFluxApproximation.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geosx
{

class TwoPointFluxApproximation : public FluxApproximationBase
{
public:

  static std::string CatalogName() { return "TwoPointFluxApproximation"; }

  TwoPointFluxApproximation() = delete;

  TwoPointFluxApproximation(std::string const & name, dataRepository::ManagedGroup * const parent);

protected:

  virtual void computeCellStencil( DomainPartition const & domain ) override;

  virtual void addToFractureStencil( DomainPartition const & domain,
                                     string const & faceElementRegionName ) override;

  virtual void computeBoundaryStencil( DomainPartition const & domain,
                                       set<localIndex> const & faceSet,
                                       BoundaryStencil & stencil ) override;

};

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
