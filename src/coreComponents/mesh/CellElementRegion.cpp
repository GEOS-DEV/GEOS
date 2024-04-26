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
    setDescription( GEOS_FMT( "List of the desired cellBlocks names to contain in this {}.\n"
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
                              "I.e. {{ 1, 2 }} selects the {{ 1_tetrahedra, 1_pyramid, 2_tetrahedra, 2_pyramid... }} cellBlocks.\n"
                              "This setting cannot be used simultaneously with {} nor {}.",
                              catalogName(),
                              viewKeyStruct::sourceCellBlockNamesString(),
                              viewKeyStruct::cellBlockMatchPatternsString() ) );

  registerWrapper( viewKeyStruct::cellBlockMatchPatternsString(), &m_cellBlockMatchPatterns ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( GEOS_FMT( "List of fnmatch pattern to match cellBlock names to add them in this {}.\n"
                              "I.e. \"{{ * }}\" selects every elements, {{ [1-5]_* }} selects the "
                              "{{ 1_tetrahedra, 2_tetrahedra, ..., 5_tetrahedra, 1_pyramid... }} cellBlocks.\n"
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

std::set< string > getAvailableAttributeValues( std::set< string > const & cellBlocksNames )
{
  std::set< string > attributeValues;
  for( auto const & cellBlockName : cellBlocksNames )
  {
    attributeValues.insert( CellElementRegion::getCellBlockAttributeValue( cellBlockName ) );
  }
  return attributeValues;
}

std::set< string >
CellElementRegion::computeSelectedCellBlocks( std::set< string > const & cellBlocksNames ) const
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
      for( auto const & cellBlockName : cellBlocksNames )
      {
        // if the pattern matches the tested cellBlock name
        if( fnmatch( matchPattern.c_str(), cellBlockName.c_str(), 0 ) == 0 )
        {
          matching = true;
          selectedCellBlocks.insert( cellBlockName );
        }
      }

      if( !matching )
      {
        // a) in the case where the user is requesting attribute values
        if( !m_cellBlockAttributeValues.empty() )
        {
          auto const attributeValueLength = matchPattern.size() - attributeValuePatternSuffix.size();
          auto const attributeValue = matchPattern.substr( 0, attributeValueLength );
          GEOS_THROW( GEOS_FMT( "{}: Attribute value '{}' not found.\nAvailable attribute list: {{ {} }}",
                                getWrapperDataContext( viewKeyStruct::cellBlockAttributeValuesString() ),
                                attributeValue,
                                stringutilities::join( getAvailableAttributeValues( cellBlocksNames ), ", " ) ),
                      InputError );
        }
        else
        {
          // b) in the case where the user is requesting match pattern(s)
          GEOS_THROW_IF( !m_cellBlockMatchPatterns.empty(),
                         GEOS_FMT( "{}: No cellBlock name is satisfying the pattern '{}'.\nAvailable cellBlock list: {{ {} }}",
                                   getWrapperDataContext( viewKeyStruct::cellBlockMatchPatternsString() ),
                                   matchPattern,
                                   stringutilities::join( cellBlocksNames, ", " ) ),
                         InputError );
        }
      }
    }
  }

  // classic one-by-one cellblock selecting
  {
    for( string const & requestedCellBlockName : m_cellBlockNames )
    {
      GEOS_THROW_IF( cellBlocksNames.count( requestedCellBlockName ) == 0,
                     GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}",
                               getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                               requestedCellBlockName,
                               stringutilities::join( cellBlocksNames, ", " ) ),
                     InputError );
      selectedCellBlocks.insert( requestedCellBlockName );
    }
  }

  return selectedCellBlocks;
}


std::set< string >
CellElementRegion::generateMesh( Group const & cellBlocks,
                                 std::set< string > cellBlocksNames,
                                 std::map< string, CellElementRegion const * > & cellBlocksMatchers )
{
  std::set< string > selectedCellBlocks = computeSelectedCellBlocks( cellBlocksNames );

  m_cellBlockNames.clear();
  m_cellBlockNames.insert( 0, selectedCellBlocks.cbegin(), selectedCellBlocks.cend() );

  for( string const & cellBlockName : m_cellBlockNames )
  {
    // Duplication detection :
    // For now, multiple regions per cell is not supported (by ElementRegionManager::getCellBlockToSubRegionMap())
    // TODO: refactor the CellElementRegion & Mesh classes so regions are mapped to cellblocks IN the mesh (and potencially
    // to multiple regions per cell). So, for external meshes, the cellblocks would no longer be exposed to the final user.
    auto const [existingMatcherIt, isNewMatcher] = cellBlocksMatchers.emplace( cellBlockName, this );
    GEOS_THROW_IF( !isNewMatcher,
                   GEOS_FMT( "The cellBlock '{}' has been referenced in multiple {}:\n- {}\n- {}",
                             cellBlockName, catalogName(),
                             existingMatcherIt->second->getDataContext(), getDataContext() ),
                   InputError );

    CellBlockABC const * cellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cellBlockName );
    // should not be thrown as computeSelectedCellBlocks() should check every cellBlock names.
    GEOS_ERROR_IF( cellBlock == nullptr,
                   GEOS_FMT( "{}: Unexpected error, {} not found.\nAvailable cellBlock list: {{ {} }}",
                             getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ), cellBlockName,
                             stringutilities::join( cellBlocksNames, ", " ) ) );
    registerSubRegion( *cellBlock );
  }

  return selectedCellBlocks;
}
void CellElementRegion::generateMesh( Group const & )
{
  GEOS_ERROR( "Not implemented, use the overloading of this method." );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
