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
 * @file SimpleWell.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP

#include "WellBase.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
static constexpr auto simpleWell = "SimpleWell";
}
}

class DomainPartition;

class SimpleWell : public WellBase
{
public:
  explicit SimpleWell( string const & name, dataRepository::ManagedGroup * const parent );
  ~SimpleWell() override;

  SimpleWell() = delete;
  SimpleWell( SimpleWell const &) = delete;
  SimpleWell( SimpleWell && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return dataRepository::keys::simpleWell; }

  const string getCatalogName() const override;

  void FillDocumentationNode() override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;

  virtual void FinalInitialization( ManagedGroup * const group ) override;

  localIndex numConnectionsGlobal() const { return numSubGroups(); }
  localIndex numConnectionsLocal() const { return size(); }

  struct viewKeyStruct : public WellBase::viewKeyStruct
  {

    static constexpr auto connectionElementRegionString = "connectionElementRegion";
    static constexpr auto connectionElementSubregionString = "connectionElementSubregion";
    static constexpr auto connectionElementIndexString = "connectionElementIndex";

    static constexpr auto pressureString = "pressure";
    static constexpr auto transmissibilityString = "transmissibility";
    static constexpr auto gravityDepthString = "gravityDepth";

    dataRepository::ViewKey connectionElementRegion    = { connectionElementRegionString    };
    dataRepository::ViewKey connectionElementSubregion = { connectionElementSubregionString };
    dataRepository::ViewKey connectionElementIndex     = { connectionElementIndexString     };

    dataRepository::ViewKey pressure         = { pressureString         };
    dataRepository::ViewKey transmissibility = { transmissibilityString };
    dataRepository::ViewKey gravityDepth     = { gravityDepthString     };

  } viewKeys;

private:

  void ConnectToCells( DomainPartition const * domain );
  void PrecomputeData( DomainPartition const * domain );

  array1d<localIndex> m_connectionElementRegion;
  array1d<localIndex> m_connectionElementSubregion;
  array1d<localIndex> m_connectionElementIndex;

};

} //namespace geosx


#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_SIMPLEWELL_HPP
