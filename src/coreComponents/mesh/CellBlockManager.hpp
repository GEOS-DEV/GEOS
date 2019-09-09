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
 * @file CellBlockManager.hpp
 */

#ifndef CELLBLOCKMANAGER_H_
#define CELLBLOCKMANAGER_H_

#include "CellBlock.hpp"

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
 * Class to manage the data stored at the element level.
 */
class CellBlockManager : public ObjectManagerBase
{
public:
  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  static string CatalogName()
  {
    return "CellBlockManager";
  }

  virtual const string getCatalogName() const override final
  { return CellBlockManager::CatalogName(); }



  ///@}

  CellBlockManager( string const &, Group * const parent );
  virtual ~CellBlockManager() override;


//  void Initialize(  ){}

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  using Group::resize;

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
    Group * elementRegions = this->GetGroup(dataRepository::keys::cellBlocks);
    elementRegions->forSubGroups<CellBlock>( lambda );
  }
private:
  CellBlockManager( const CellBlockManager& );
  CellBlockManager& operator=( const CellBlockManager&);


};
}
#endif /* ELEMENTMANAGERT_H_ */
