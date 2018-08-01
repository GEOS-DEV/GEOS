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
 * @file FluxApproximationBase.hpp
 *
 */

#ifndef SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
#define SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_

#include <dataRepository/ManagedGroup.hpp>
#include <finiteVolume/StencilCollection.hpp>
#include <managers/DomainPartition.hpp>

namespace geosx
{

/**
 * @struct A structure containing a single cell (element) identifier triplet
 */
struct CellDescriptor
{
  localIndex region;
  localIndex subRegion;
  localIndex index;
};

/**
 * @class FluxApproximationBase
 *
 * Base class for various flux approximation classes.
 * Stores the stencil, construction is implemented in derived classes
 */
class FluxApproximationBase : public dataRepository::ManagedGroup
{
public:

  using CatalogInterface = cxx_utilities::CatalogInterface<FluxApproximationBase, string const &, ManagedGroup * const >;
  static typename CatalogInterface::CatalogType& GetCatalog();

  using StencilType = StencilCollection<CellDescriptor, real64>;

  void FillDocumentationNode() override;

  FluxApproximationBase() = delete;

  FluxApproximationBase(string const & name, dataRepository::ManagedGroup * const parent);

  void ReadXML_PostProcess() override;

  /// provides access to the stencil collection
  StencilType const & getStencil() const { return m_stencilCellToCell; }

  /// actual computation of the stencil, to be overridden by implementations
  virtual void compute(DomainPartition * domain) = 0;

  struct viewKeyStruct
  {
    static constexpr auto fieldNameString = "fieldName";
  };

protected:

  StencilType m_stencilCellToCell;

  /// name of the coefficient field
  string m_fieldName;

};

}

#endif //SRC_COMPONENTS_CORE_SRC_FINITEVOLUME_FLUXAPPROXIMATIONBASE_HPP_
