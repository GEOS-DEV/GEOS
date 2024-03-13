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
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( GEOS_FMT( "The list of cell-blocks this {} contains.\n"
                              "If the source mesh has a data array named as the VTKMesh::regionAttribute value, the list can contain:\n"
                              "  - a list of sub-region(s) (= cell-blocks), each refered as \"regionAttribute_elementType\" (ie: \"1_tetrahedra\"),\n"
                              "  - a list of region(s), each refered as \"regionAttribute\",\n"
                              "  - all ellement of the mesh at once, refered with \"all\".\n"
                              "If the source mesh has no \"regionAttribute\" data array, the list can contain:\n"
                              "  - a list of sub-region(s) (= cell-blocks), each refered as \"elementType\" (ie: \"tetrahedra\"),\n"
                              "  - all ellement of the mesh at once, refered with \"all\".\n",
                              catalogName() ) );

  //TODO: add documentation ?
  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 );
}

CellElementRegion::~CellElementRegion()
{}


bool CellElementRegion::isSelectingAllCells()
{
  if( m_cellBlockNames.size() == 1 &&
      m_cellBlockNames[0] == CellElementRegion::viewKeyStruct::selectAllCellBlocksString() )
  {
    return true;
  }
  else
  {
    // error if the all keyword is mixed with other cellBlocks names
    for( string const & cellBlockName : m_cellBlockNames )
    {
      GEOS_THROW_IF( cellBlockName == viewKeyStruct::selectAllCellBlocksString(),
                     GEOS_FMT( "{0}: The keyword '{1}' is useful to include all existing cells in a {2},"
                               " it should be used alone in a unique {2}.",
                               getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                               viewKeyStruct::selectAllCellBlocksString(),
                               CellElementRegion::catalogName() ),
                     InputError );
    }
    return false;
  }
}

void CellElementRegion::registerSubRegion( CellBlockABC const & cellBlock )
{
  GEOS_THROW_IF( cellBlock.getName() == CellElementRegion::viewKeyStruct::selectAllCellBlocksString(),
                 GEOS_FMT( "{}: cellBlock named '{}' uses the reserved keyword '{}'.",
                           getWrapperDataContext( CellElementRegion::viewKeyStruct::sourceCellBlockNamesString() ),
                           cellBlock.getName(),
                           CellElementRegion::viewKeyStruct::selectAllCellBlocksString() ),
                 InputError );

  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );
  // For now, subRegion name must be the same as the cellBlock.
  CellElementSubRegion & subRegion =
    elementSubRegions.registerGroup< CellElementSubRegion >( cellBlock.getName() );

  subRegion.copyFromCellBlock( cellBlock );
}

void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  // if we would like to select all cellBlocks on an unique CellElementRegion
  if( isSelectingAllCells() )
  {
    cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
    {
      registerSubRegion( cellBlock );
    } );
  }
  else
  {
    for( string const & cellBlockName : m_cellBlockNames )
    {
      // if we want to add a specific cellBlock ("1_tetrahedra" form typically)
      CellBlockABC const * exactCellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cellBlockName );
      if( exactCellBlock != nullptr )
      {
        registerSubRegion( *exactCellBlock );
      }
      else
      {
        // TODO: replace commented code by new cellBlockMatch & cellBlockAttributeValues wrappers
        // bool foundOneCellBlock = false;

        // if we want to add all region
        // cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & foundCellBlock )
        // {
        //   if( foundCellBlock.getRegionName() == cellBlockName )
        //   {
        //     registerSubRegion( foundCellBlock );
        //     foundOneCellBlock = true;
        //   }
        // } );

        // if( !foundOneCellBlock )
        // {
        // std::set< string > regions;
        std::vector< string > subRegions;
        cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cb )
        {
          // regions.insert( string( cb.getRegionName() ) );
          subRegions.push_back( string( cb.getName() ) );
        } );
        GEOS_THROW( GEOS_FMT( "{}: region or sub-region '{}' not found.\n"
                              // "Available region list: {}\n"
                              "Available sub-region list: {}",
                              getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                              cellBlockName,
                              // stringutilities::join( regions, ", " ),
                              stringutilities::join( subRegions, ", " ) ),
                    InputError );
        // }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
