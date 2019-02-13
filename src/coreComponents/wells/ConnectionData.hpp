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
 * @file ConnectionData.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_CONNECTIONDATA_HPP
#define GEOSX_CORECOMPONENTS_WELLS_CONNECTIONDATA_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto connections = "Connections";
}
}

class DomainPartition;
class Connection;

class ConnectionData : public ObjectManagerBase
{
public:

  explicit ConnectionData( string const & name, dataRepository::ManagedGroup * const parent );
  ~ConnectionData() override;

  ConnectionData() = delete;
  ConnectionData( ConnectionData const &) = delete;
  ConnectionData( ConnectionData && ) = delete;

  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  virtual const string getCatalogName() const override;

  localIndex numConnectionsGlobal() const { return integer_conversion<localIndex>(m_connectionList.size()); }
  localIndex numConnectionsLocal()  const { return integer_conversion<localIndex>(size());         }

  Connection const * getConnection( localIndex iconn ) const;
  Connection *       getConnection( localIndex iconn );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {

    static constexpr auto connectionIndexString = "connectionIndex";

    using ViewKey = dataRepository::ViewKey;
    
    ViewKey connectionIndex = { connectionIndexString };

  } viewKeysConnectionData;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto connectionString = "Connection";

    dataRepository::GroupKey connection = { connectionString };

  } groupKeysConnectionData;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  void PrecomputeData( MeshLevel const * domain );

  array1d<localIndex> m_connectionIndex;

  string_array m_connectionList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_CONNECTIONDATA_HPP
