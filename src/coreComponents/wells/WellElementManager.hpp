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
 * @file WellElementManager.hpp
 *
 */

#ifndef GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP
#define GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP

#include "dataRepository/ManagedGroup.hpp"
#include "managers/ObjectManagerBase.hpp"

namespace geosx
{

class MeshLevel;

namespace dataRepository
{
namespace keys
{
static constexpr auto wellElements = "Segments";
}
}

class DomainPartition;
class WellElement;

class WellElementManager : public ObjectManagerBase
{
public:

  explicit WellElementManager( string const & name, dataRepository::ManagedGroup * const parent );
  ~WellElementManager() override;

  WellElementManager() = delete;
  WellElementManager( WellElementManager const &) = delete;
  WellElementManager( WellElementManager && ) = delete;

  virtual const string getCatalogName() const override;
  
  dataRepository::ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  globalIndex numWellElementsGlobal()  const
  { return integer_conversion<globalIndex>(m_globalWellElementList.size()); }

  WellElement const * getWellElement( globalIndex iwelem ) const;
  WellElement *       getWellElement( globalIndex iwelem );

  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
  } viewKeysWellElementManager;

  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    static constexpr auto wellElementString = "Segment";

    dataRepository::GroupKey wellElement = { wellElementString };

  } groupKeysWellElementManager;

protected:

  virtual void InitializePreSubGroups( ManagedGroup * const problemManager ) override;

  virtual void InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager ) override;

private:

  string_array m_globalWellElementList;

};

} //namespace geosx

#endif //GEOSX_CORECOMPONENTS_WELLS_WELLELEMENTMANAGER_HPP
