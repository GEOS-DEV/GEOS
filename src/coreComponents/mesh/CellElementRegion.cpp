/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "CellElementRegion.hpp"
#include "CellElementSubRegion.hpp"
#include "mesh/generators/CellBlockABC.hpp"

namespace geos
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceCellBlockNamesString(), &m_cellBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL );
}

CellElementRegion::~CellElementRegion()
{}


bool isSelectingAllCellsInOneRegion( string_array const & cellBlockNames, CellElementRegion const & region )
{
  if( cellBlockNames.size() == 1 &&
      cellBlockNames[0] == CellElementRegion::viewKeyStruct::selectAllCellBlocksString())
  {
    // check that we have only one cellBlock included
    int n = 0;
    region.getParent().forSubGroups< CellElementRegion >( [&] ( CellElementRegion const & elemRegion )
    {
      n += elemRegion.getCellBlockNames().size();
    } );
    return n == 1;
  }
  else
  {
    return false;
  }
}

void registerSubRegion( CellElementRegion const & region, Group & elementSubRegions,
                        CellBlockABC const & cellBlock )
{
  GEOS_THROW_IF( cellBlock.getName() == CellElementRegion::viewKeyStruct::selectAllCellBlocksString() ||
                 cellBlock.getRegionName() == CellElementRegion::viewKeyStruct::selectAllCellBlocksString(),
                 GEOS_FMT( "{}: cellBlock named '{}' uses the reserved keyword '{}'.",
                           region.getWrapperDataContext( CellElementRegion::viewKeyStruct::sourceCellBlockNamesString() ),
                           cellBlock.getName(),
                           CellElementRegion::viewKeyStruct::selectAllCellBlocksString() ),
                 InputError );

  // For now, subRegion name must be the same as the cellBlock.
  CellElementSubRegion & subRegion =
    elementSubRegions.registerGroup< CellElementSubRegion >( cellBlock.getName() );

  subRegion.copyFromCellBlock( cellBlock );
}

void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );

  // if we would like to select all cellBlocks on an unique CellElementRegion
  if( isSelectingAllCellsInOneRegion( m_cellBlockNames, *this ) )
  {
    cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
    {
      registerSubRegion( *this, elementSubRegions, cellBlock );
    } );
  }
  else
  {
    for( string const & cellBlockName : this->m_cellBlockNames )
    {
      GEOS_THROW_IF( cellBlockName == viewKeyStruct::selectAllCellBlocksString(),
                     GEOS_FMT( "{0}: The keyword '{1}' is useful to include all existing cells in a {2},"
                               " it should be used alone in a unique {2}.",
                               getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                               viewKeyStruct::selectAllCellBlocksString(),
                               CellElementRegion::catalogName() ),
                     InputError );

      // if we want to add a specific cellBlock ("1_tetrahedra" form typically)
      CellBlockABC const * exactCellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cellBlockName );
      if( exactCellBlock != nullptr )
      {
        registerSubRegion( *this, elementSubRegions, *exactCellBlock );
      }
      else
      {
        bool foundOneCellBlock = false;

        // if we want to add all region
        cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & foundCellBlock )
        {
          if( foundCellBlock.getRegionName() == cellBlockName )
          {
            registerSubRegion( *this, elementSubRegions, foundCellBlock );
            foundOneCellBlock = true;
          }
        } );

        if( !foundOneCellBlock )
        {
          std::set< string > regions;
          std::vector< string > subRegions;
          cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cb )
          {
            regions.insert( string( cb.getRegionName() ) );
            subRegions.push_back( string( cb.getName() ) );
          } );
          GEOS_THROW( GEOS_FMT( "{}: region or sub-region '{}' not found.\n"
                                "Available region list: {}\n"
                                "Available sub-region list: {}",
                                getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                                cellBlockName,
                                stringutilities::join( regions, ", " ),
                                stringutilities::join( subRegions, ", " ) ),
                      InputError );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
