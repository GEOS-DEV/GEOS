/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
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

  void computeCellStencil(DomainPartition const * domain,
                          CellStencil & stencil) const override;

  void computeFaceStencil(DomainPartition const * domain,
                          set<localIndex> const & faceSet,
                          FaceStencil & stencil) const override;

  virtual void computeWellStencil( DomainPartition const * domain,
                                   WellBase const * well,
                                   WellStencil & stencil ) const override;

};

}


#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
