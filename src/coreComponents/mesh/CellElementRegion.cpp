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

#include <fnmatch.h>

namespace geos
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  std::vector< string > elementNames;
  for( int i = 0; i < numElementTypes(); ++i )
  {
    if( getElementDim( (ElementType)i ) == 3 )
    {
      elementNames.push_back( getElementTypeName( (ElementType)i ) );
    }
  }

  registerWrapper( viewKeyStruct::sourceCellBlockNamesString(), &m_cellBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( GEOS_FMT( "List of the desired cell-blocks from the mesh to contain in this {}.\n"
                              "The form of this attribute is of \"regionAttribute_elementType\", so \"1_tetrahedra\" select the "
                              "cellBlock that contains the tetrahedric elements for which the regionAttribute is 1.\n"
                              "The element types are: {}.\n"
                              "This setting cannot be used simultaneously with {} nor {}.",
                              catalogName(), stringutilities::join( elementNames, ", " ),
                              viewKeyStruct::cellBlockAttributeValuesString(),
                              viewKeyStruct::cellBlockMatchPatternsString() ) );

  registerWrapper( viewKeyStruct::cellBlockAttributeValuesString(), &m_cellBlockAttributeValues ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( GEOS_FMT( "List of regionAttribute values for which we want to add the cells in this {}.\n"
                              "I.e. {{ 1 }} selects the {{ 1_tetrahedra, 1_pyramid, 1_hexahedra... }} cellBlocks.\n"
                              "This setting cannot be used simultaneously with {} nor {}.",
                              catalogName(),
                              viewKeyStruct::sourceCellBlockNamesString(),
                              viewKeyStruct::cellBlockMatchPatternsString() ) );

  registerWrapper( viewKeyStruct::cellBlockMatchPatternsString(), &m_cellBlockMatchPatterns ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( GEOS_FMT( "List of fnmatch pattern to match cellBlock names to add them in this {}.\n"
                              "I.e. {{ 1_* }} selects the {{ 1_tetrahedra, 1_pyramid, 1_hexahedra... }} cellBlocks.\n"
                              "This setting cannot be used simultaneously with {} nor {}.",
                              catalogName(),
                              viewKeyStruct::sourceCellBlockNamesString(),
                              viewKeyStruct::cellBlockAttributeValuesString() ) );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 );
}

CellElementRegion::~CellElementRegion()
{}


void CellElementRegion::registerSubRegion( CellBlockABC const & cellBlock )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );
  // For now, subRegion name must be the same as the cellBlock (so we can match them and reference them in errors).
  CellElementSubRegion & subRegion =
    elementSubRegions.registerGroup< CellElementSubRegion >( cellBlock.getName() );

  subRegion.copyFromCellBlock( cellBlock );
}

string CellElementRegion::getCellBlockAttributeValue( string_view cellBlockName )
{
  return string( stringutilities::removeStringAndFollowingContent( cellBlockName, "_" ) );
}

std::set< string > getAvailableAttributeValues( Group const & cellBlocks )
{
  std::set< string > attributeValues;
  cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
  {
    attributeValues.insert( CellElementRegion::getCellBlockAttributeValue( cellBlock.getName() ) );
  } );
  return attributeValues;
}

std::set< string > getCellBlockNamesSet( Group const & cellBlocks )
{
  std::set< string > cellBlockNames;
  cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
  {
    cellBlockNames.insert( cellBlock.getName() );
  } );
  return cellBlockNames;
}

std::set< string > CellElementRegion::computeSelectedCellBlocks( Group const & cellBlocks ) const
{
  std::set< string > selectedCellBlocks;

  GEOS_THROW_IF( m_cellBlockAttributeValues.empty() && m_cellBlockMatchPatterns.empty() && m_cellBlockNames.empty(),
                 GEOS_FMT( "{}: cellBlocks must be selected to fill this region (using {}, {} or {}).",
                           getDataContext(),
                           viewKeyStruct::sourceCellBlockNamesString(),
                           viewKeyStruct::cellBlockAttributeValuesString(),
                           viewKeyStruct::cellBlockMatchPatternsString() ),
                 InputError );
  GEOS_THROW_IF( ( !m_cellBlockAttributeValues.empty() + !m_cellBlockMatchPatterns.empty() + !m_cellBlockNames.empty() ) != 1,
                 GEOS_FMT( "{}: Only one setting must be used to select cellBlocks ({}, {} or {}).",
                           getDataContext(),
                           viewKeyStruct::sourceCellBlockNamesString(),
                           viewKeyStruct::cellBlockAttributeValuesString(),
                           viewKeyStruct::cellBlockMatchPatternsString() ),
                 InputError );

  {
    std::set< string > cellBlockMatchPatterns;

    // for each desired attribute values, we add the following the match patterns:
    // Attribute value, followed by an underscore, followed by one or more characters.
    static string_view constexpr attributeValuePatternSuffix = "_?*";
    for( integer attributeValue : m_cellBlockAttributeValues )
    {
      cellBlockMatchPatterns.insert( std::to_string( attributeValue ) + string( attributeValuePatternSuffix ) );
    }

    for( string const & matchPattern : m_cellBlockMatchPatterns )
    {
      cellBlockMatchPatterns.insert( matchPattern );
    }

    for( string const & matchPattern : cellBlockMatchPatterns )
    {
      bool matching = false;
      cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
      {
        // if the pattern matches the tested cellBlock name
        if( fnmatch( matchPattern.c_str(), cellBlock.getName().c_str(), 0 ) == 0 )
        {
          matching = true;
          selectedCellBlocks.insert( cellBlock.getName() );
        }
      } );

      if( !matching )
      {
        // a) check for the case where the user is requesting attribute values
        if( !m_cellBlockAttributeValues.empty() )
        {
          auto const attributeValueLength = matchPattern.size() - attributeValuePatternSuffix.size();
          auto const attributeValue = matchPattern.substr( 0, attributeValueLength );
          GEOS_THROW( GEOS_FMT( "{}: Attribute value '{}' not found.\nAvailable attribute list: {{ {} }}",
                                getWrapperDataContext( viewKeyStruct::cellBlockAttributeValuesString() ),
                                attributeValue,
                                stringutilities::join( getAvailableAttributeValues( cellBlocks ), ", " ) ),
                      InputError );
        }
        // b) check for the case where the user is requesting match pattern
        GEOS_THROW_IF( !m_cellBlockMatchPatterns.empty(),
                       GEOS_FMT( "{}: No cellBlock name is satisfying the pattern '{}'.\nAvailable cellBlock list: {{ {} }}",
                                 getWrapperDataContext( viewKeyStruct::cellBlockMatchPatternsString() ),
                                 matchPattern,
                                 stringutilities::join( getCellBlockNamesSet( cellBlocks ), ", " ) ),
                       InputError );
        // this error exists just in case, we should not fall there.
        GEOS_ERROR( GEOS_FMT( "{}: Unknown error, no cellBlock matching the pattern '{}'.\nAvailable cellBlock list: {{ {} }}",
                              matchPattern, stringutilities::join( getCellBlockNamesSet( cellBlocks ), ", " ) ) );
      }
    }
  }

  // classic one-by-one cellblock selecting
  {
    for( string const & requestedCellBlockName : m_cellBlockNames )
    {
      bool found = false;
      cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
      {
        if( cellBlock.getName() == requestedCellBlockName )
        {
          selectedCellBlocks.insert( requestedCellBlockName );
          found = true;
        }
      } );
      GEOS_THROW_IF( !found,
                     GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}",
                               getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                               requestedCellBlockName,
                               stringutilities::join( getCellBlockNamesSet( cellBlocks ), ", " ) ),
                     InputError );
    }
  }

  return selectedCellBlocks;
}


void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  {
    std::set< string > selectedCellBlocks = computeSelectedCellBlocks( cellBlocks );
    m_cellBlockNames.clear();
    m_cellBlockNames.insert( 0, selectedCellBlocks.cbegin(), selectedCellBlocks.cend() );
  }

  for( string const & cellBlockName : m_cellBlockNames )
  {
    CellBlockABC const * cellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cellBlockName );
    // should not be thrown as computeSelectedCellBlocks() should check every cellBlock names.
    GEOS_ERROR_IF( cellBlock == nullptr,
                   GEOS_FMT( "{}: Unexpected error with {} named '{}'.\nAvailable cellBlock list: {{ {} }}",
                             getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                             cellBlockName,
                             stringutilities::join( getCellBlockNamesSet( cellBlocks ), ", " ) ) );

    registerSubRegion( *cellBlock );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
