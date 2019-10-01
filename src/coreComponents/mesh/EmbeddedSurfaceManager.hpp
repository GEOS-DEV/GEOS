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
 * @file ElementManagerT.h
 */

#ifndef EMBEDDEDSURFACEMANAGER_H_
#define EMBEDDEDSURFACEMANAGER_H_

//#include "Common.h"
//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "EmbeddedSurfaceSubRegion.hpp"
#include "managers/ObjectManagerBase.hpp"

//#include "legacy/ArrayT/bufvector.h"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const cellBlocks = "cellBlocks";
}
}


/**
 * Class to manage the data stored for each Embedded surface element.
 */
class EmbeddedSurfaceManager : public ObjectManagerBase
{
public:
  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static string CatalogName()
  {
    return "EmbeddedSurfaceManager";
  }

  virtual const string getCatalogName() const override final
  { return EmbeddedSurfaceManager::CatalogName(); }



  ///@}

  EmbeddedSurfaceManager( string const &, ManagedGroup * const parent );
  virtual ~EmbeddedSurfaceManager() override;


//  void Initialize(  ){}

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  using ManagedGroup::resize;

  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

//  CellBlock & CreateRegion( string const & regionName,
//                               string const & elementType,
//                               integer const & numElements );

  CellBlock * GetRegion( string const & regionName )
  {
    return this->GetGroup(dataRepository::keys::cellBlocks)->GetGroup<CellBlock>(regionName);
  }

  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    ManagedGroup * elementRegions = this->GetGroup(dataRepository::keys::cellBlocks);
    elementRegions->forSubGroups<CellBlock>( lambda );
  }
private:
  EmbeddedSurfaceManager( const EmbeddedSurfaceManager& );
  EmbeddedSurfaceManager& operator=( const EmbeddedSurfaceManager&);


};
}
#endif /* ELEMENTMANAGERT_H_ */
