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
 * @file Well.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_
#define GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_

#include "managers/ObjectManagerBase.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "WellElementManager.hpp"
#include "ConnectionData.hpp"
#include "PerforationData.hpp"
#include "PerforationManager.hpp"

namespace geosx
{

class Well : public ObjectManagerBase
{
public:

  enum class Type { PRODUCER, INJECTOR };

  explicit Well( string const & name, dataRepository::ManagedGroup * const parent );
  ~Well() override;

  Well() = delete;
  Well( Well const & ) = delete;
  Well( Well && ) = delete;

  /// Catalog name interface
  static string CatalogName() { return "Well"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< Well, std::string const &, ManagedGroup * const >;

  static CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual const string getCatalogName() const override { return CatalogName(); }

  virtual dataRepository::ManagedGroup * CreateChild(string const & childKey, string const & childName) override;

  virtual void InitializePostSubGroups( ManagedGroup * const rootGroup ) override;

  virtual void PostProcessInput() override;

  localIndex numWellElementsLocal() const { return m_wellElementSubRegion.size(); }
  
  localIndex numConnectionsLocal()  const { return m_connectionData.numConnectionsLocal(); }

  localIndex numPerforationsLocal()  const { return m_perforationData.numPerforationsLocal(); }

  WellElementSubRegion * getWellElements()             { return &m_wellElementSubRegion; }
  WellElementSubRegion const * getWellElements() const { return &m_wellElementSubRegion; }
  
  ConnectionData * getConnections()             { return &m_connectionData; }
  ConnectionData const * getConnections() const { return &m_connectionData; }
  
  PerforationData * getPerforations()             { return &m_perforationData; }
  PerforationData const * getPerforations() const { return &m_perforationData; }
  
  R1Tensor const & getGravityVector() const;

  real64 getReferenceDepth() const { return m_referenceDepth; }
  Type getType() const { return m_type; }

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto referenceDepthString = "referenceDepth";
    static constexpr auto typeString           = "type";

    dataRepository::ViewKey referenceDepth             = { referenceDepthString };
    dataRepository::ViewKey type                       = { typeString };

  } viewKeysWell;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto wellElementsString    = dataRepository::keys::wellElements;
    static constexpr auto wellElementDataString = dataRepository::keys::wellElementData;
    static constexpr auto connectionDataString  = dataRepository::keys::connectionData;
    static constexpr auto perforationsString    = dataRepository::keys::perforations;
    static constexpr auto perforationDataString = dataRepository::keys::perforationData;

    dataRepository::GroupKey wellElements    = { wellElementsString };
    dataRepository::GroupKey wellElementData = { wellElementDataString };
    dataRepository::GroupKey connectionData  = { connectionDataString };
    dataRepository::GroupKey perforations    = { perforationsString };
    dataRepository::GroupKey perforationData = { perforationDataString };

  } groupKeysWell;

protected:

  WellElementSubRegion m_wellElementSubRegion;
  WellElementManager m_wellElementManager;
  
  ConnectionData  m_connectionData;

  PerforationData    m_perforationData;
  PerforationManager m_perforationManager;
  
  real64 m_referenceDepth;
  string m_typeString;
  Type   m_type;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_WELL_HPP_
