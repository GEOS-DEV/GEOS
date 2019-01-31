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
 * @file WellBase.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_

#include "managers/ObjectManagerBase.hpp"
#include "SegmentManager.hpp"
#include "ConnectionManager.hpp"
#include "PerforationManager.hpp"


namespace geosx
{

class WellBase : public ObjectManagerBase
{
public:

  enum class Type { PRODUCER, INJECTOR };

  explicit WellBase( string const & name, dataRepository::ManagedGroup * const parent );
  ~WellBase() override;

  WellBase() = delete;
  WellBase( WellBase const &) = delete;
  WellBase( WellBase && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return "WellBase"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< WellBase, std::string const &, ManagedGroup * const >;

  static CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual const string getCatalogName() const override { return CatalogName(); }

  virtual dataRepository::ManagedGroup * CreateChild(string const & childKey, string const & childName) override;

  virtual void PostProcessInput() override;

  localIndex numSegmentsGlobal() const { return m_segmentManager.numSegmentsGlobal(); }
  localIndex numSegmentsLocal()  const { return m_segmentManager.numSegmentsLocal(); }

  localIndex numConnectionsGlobal() const { return m_connectionManager.numConnectionsGlobal(); }
  localIndex numConnectionsLocal()  const { return m_connectionManager.numConnectionsLocal();  }

  localIndex numPerforationsGlobal() const { return m_perforationManager.numPerforationsGlobal(); }
  localIndex numPerforationsLocal()  const { return m_perforationManager.numPerforationsLocal();  }
 
  SegmentManager * geSegments()              { return &m_segmentManager; }
  SegmentManager const * getSegments() const { return &m_segmentManager; }

  ConnectionManager * getConnections()             { return &m_connectionManager; }
  ConnectionManager const * getConnections() const { return &m_connectionManager; }
  
  PerforationManager * getPerforations()             { return &m_perforationManager; }
  PerforationManager const * getPerforations() const { return &m_perforationManager; }
  
  R1Tensor const & getGravityVector() const;

  real64 getReferenceDepth() const { return m_referenceDepth; }
  Type getType() const { return m_type; }

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto referenceDepthString = "referenceDepth";
    static constexpr auto typeString           = "type";

    dataRepository::ViewKey referenceDepth             = { referenceDepthString };
    dataRepository::ViewKey type                       = { typeString };

  } viewKeysWellBase;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {

    static constexpr auto segmentsString     = dataRepository::keys::segments;
    static constexpr auto connectionsString  = dataRepository::keys::connections;
    static constexpr auto perforationsString = dataRepository::keys::perforations;

    dataRepository::GroupKey segments     = { segmentsString     };
    dataRepository::GroupKey connections  = { connectionsString  };
    dataRepository::GroupKey perforations = { perforationsString };

  } groupKeysWellBase;

protected:

  SegmentManager     m_segmentManager;
  ConnectionManager  m_connectionManager;
  PerforationManager m_perforationManager;

  real64 m_referenceDepth;
  string m_typeString;
  Type   m_type;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_WELLBASE_HPP_
