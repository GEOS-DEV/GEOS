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
 * @file PerforationManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP
#define GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto perforations = "Perforations";
}
}

class DomainPartition;
class Perforation;

class PerforationManager : public ObjectManagerBase
{
public:

  explicit PerforationManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~PerforationManager() override;

  PerforationManager() = delete;
  PerforationManager( PerforationManager const &) = delete;
  PerforationManager( PerforationManager && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  globalIndex numPerforationsGlobal()  const
  { return integer_conversion<globalIndex>(m_globalPerforationList.size()); }

  Perforation const * getPerforation( globalIndex iperf ) const;
  Perforation *       getPerforation( globalIndex iperf );
  
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
  } viewKeysPerforationManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto perforationString = "Perforation";

    dataRepository::GroupKey perforation = { perforationString };

  } groupKeysPerforationManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  string_array m_globalPerforationList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_PERFORATIONMANAGER_HPP
