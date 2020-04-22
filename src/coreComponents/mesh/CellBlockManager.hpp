/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellBlockManager.hpp
 */

#ifndef GEOSX_MESH_CELLBLOCKMANAGER_H_
#define GEOSX_MESH_CELLBLOCKMANAGER_H_

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
    return this->GetGroup( dataRepository::keys::cellBlocks )->GetGroup< CellBlock >( regionName );
  }

  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    Group * elementRegions = this->GetGroup( dataRepository::keys::cellBlocks );
    elementRegions->forSubGroups< CellBlock >( lambda );
  }
private:
  CellBlockManager( const CellBlockManager & );
  CellBlockManager & operator=( const CellBlockManager & );


};
}
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
