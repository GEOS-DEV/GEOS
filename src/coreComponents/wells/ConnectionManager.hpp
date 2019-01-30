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
 * @file ConnectionManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTIONMANAGER_HPP
#define GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTIONMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto segments = "Segments";
}
}

class DomainPartition;
class Connection;

class ConnectionManager : public ObjectManagerBase
{
public:

  explicit ConnectionManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~SegmentManager() override;

  ConnectionManager() = delete;
  ConnectionManager( ConnectionManager const &) = delete;
  ConnectionManager( ConnectionManager && ) = delete;

  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual const string getCatalogName() const override;

  localIndex numConnectionsGlobal() const { return integer_conversion<localIndex>(m_connList.size()); }
  localIndex numConnectionsLocal()  const { return integer_conversion<localIndex>(size());         }

  Connection const * getConnection( localIndex iseg ) const;
  Connection *       getConnection( localIndex iseg );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto connIndexString    = "connIndex";

    static constexpr auto gravityDepthString = "gravityDepth";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey connIndex    = { connIndexString };

    ViewKey gravityDepth = { gravityDepthString };

  } viewKeysPerfManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto connectionString = "Connection";

    dataRepository::GroupKey connection = { connectionString };

  } groupKeysSegmentManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  void PrecomputeData( MeshLevel const * domain );

  array1d<localIndex> m_connIndex;

  array1d<real64> m_gravityDepth;

  string_array m_connectionList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_MANAGERS_WELLS_CONNECTIONMANAGER_HPP
