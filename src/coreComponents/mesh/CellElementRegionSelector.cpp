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

#include "mesh/CellElementRegionSelector.hpp"

#include <fnmatch.h>
#include <optional>


namespace geos
{
using namespace dataRepository;
using ViewKeys = CellElementRegion::viewKeyStruct;


std::optional< string > getCellBlockAttributeValue( string_view cellBlockName )
{
  static constexpr string_view separator = "_";
  return cellBlockName.find( separator ) != string_view::npos ?
         stringutilities::removeStringAndFollowingContent( cellBlockName, separator ) :
         std::optional< string >( std::nullopt );
}


CellElementRegionSelector::CellElementRegionSelector( Group const & cellBlocks ):
  m_cellBlocks( cellBlocks )
{
  cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
  {
    string const name = cellBlock.getName();
    m_cellBlockNames.emplace( name );

    std::optional< string > regionAttributeValue = getCellBlockAttributeValue( name );
    if( regionAttributeValue )
    {
      m_regionAttributeValues.emplace( *regionAttributeValue );
    }
  } );
  m_orphanCellBlockNames = m_cellBlockNames;
  m_orphanRegionAttributes = m_regionAttributeValues;
}


std::set< string >
CellElementRegionSelector::getFNMatchPatterns( CellElementRegion const & region,
                                               std::set< integer > const & requestedAttributeValues,
                                               std::set< string > const & requestedMatchPatterns ) const
{
  std::set< string > matchPatterns;

  // user region attribute matching patterns creation
  for( integer attributeValue : requestedAttributeValues )
  {
    string attributeValueStr = std::to_string( attributeValue );

    GEOS_THROW_IF( m_regionAttributeValues.count( attributeValueStr ) == 0,
                   GEOS_FMT( "{}: Region attribute value '{}' not found.\nAvailable region attribute list: {{ {} }}",
                             region.getWrapperDataContext( ViewKeys::cellBlockAttributeValuesString() ),
                             attributeValueStr,
                             stringutilities::join( m_regionAttributeValues, ", " ) ),
                   InputError );

    // for each desired attribute values, we add the following the match patterns:
    // Attribute value, followed by an underscore, followed by one or more characters.
    matchPatterns.emplace( GEOS_FMT( "{}_?*", attributeValueStr ) );
  }

  // user fnMatch patterns
  for( string const & matchPattern : requestedMatchPatterns )
  {
    matchPatterns.emplace( matchPattern );
  }

  return matchPatterns;
}

std::set< string >
CellElementRegionSelector::getFNMatchSelection( CellElementRegion const & region,
                                                std::set< string > const & requestedMatchPatterns ) const
{
  std::set< string > matchedCellBlocks;
  for( string const & matchPattern : requestedMatchPatterns )
  {
    bool matching = false;
    for( auto const & cellBlockName : m_cellBlockNames )
    {
      // if the pattern matches the tested cellBlock name
      if( fnmatch( matchPattern.c_str(), cellBlockName.c_str(), 0 ) == 0 )
      {
        matching = true;
        matchedCellBlocks.emplace( cellBlockName );
      }
    }

    GEOS_THROW_IF( !matching,
                   GEOS_FMT( "{}: No cellBlock name is satisfying the pattern '{}'.\nAvailable cellBlock list: {{ {} }}",
                             region.getWrapperDataContext( ViewKeys::cellBlockMatchPatternsString() ),
                             matchPattern,
                             stringutilities::join( m_cellBlockNames, ", " ) ),
                   InputError );
  }
  return matchedCellBlocks;
}

std::set< string >
CellElementRegionSelector::getOneByOneSelection( CellElementRegion const & region,
                                                 std::set< string > const & requestedCellBlocks ) const
{
  std::set< string > oneByOneCellBlocks;
  for( string const & requestedCellBlockName : requestedCellBlocks )
  {
    GEOS_THROW_IF( m_cellBlockNames.count( requestedCellBlockName ) == 0,
                   GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}",
                             region.getWrapperDataContext( ViewKeys::sourceCellBlockNamesString() ),
                             requestedCellBlockName,
                             stringutilities::join( m_cellBlockNames, ", " ) ),
                   InputError );
    oneByOneCellBlocks.insert( requestedCellBlockName );
  }
  return oneByOneCellBlocks;
}


template< typename T >
std::set< T > toStdSet( arrayView1d< T const > const & array )
{
  return std::set< T >( array.begin(), array.end() );
}


std::set< string > CellElementRegionSelector::selectRegionCellBlocks( CellElementRegion const & region )
{
  std::set< integer > const requestedAttributeValues = toStdSet(
    region.getReference< integer_array >( ViewKeys::cellBlockAttributeValuesString() ) );
  std::set< string > const requestedMatchPatterns = toStdSet(
    region.getReference< string_array >( ViewKeys::cellBlockMatchPatternsString() ) );
  std::set< string > const requestedCellBlock = toStdSet(
    region.getReference< string_array >( ViewKeys::sourceCellBlockNamesString() ) );

  // at least one selection method is needed
  GEOS_THROW_IF( requestedAttributeValues.empty() && requestedMatchPatterns.empty() && requestedCellBlock.empty(),
                 GEOS_FMT( "{}: cellBlocks must be selected to fill this region (using {}, {} or {}).",
                           region.getDataContext(),
                           ViewKeys::sourceCellBlockNamesString(),
                           ViewKeys::cellBlockAttributeValuesString(),
                           ViewKeys::cellBlockMatchPatternsString() ),
                 InputError );

  // only one selection method allowed for now
  GEOS_THROW_IF( ( !requestedAttributeValues.empty() + !requestedMatchPatterns.empty() + !requestedCellBlock.empty() ) != 1,
                 GEOS_FMT( "{}: Only one setting must be used to select cellBlocks ({}, {} or {}).",
                           region.getDataContext(),
                           ViewKeys::sourceCellBlockNamesString(),
                           ViewKeys::cellBlockAttributeValuesString(),
                           ViewKeys::cellBlockMatchPatternsString() ),
                 InputError );

  std::set< string > const matchPatterns = getFNMatchPatterns( region, requestedAttributeValues, requestedMatchPatterns );
  std::set< string > const matchedCellBlocks = getFNMatchSelection( region, matchPatterns );
  std::set< string > cellBlocksSelection = getOneByOneSelection( region, requestedCellBlock );
  cellBlocksSelection.insert( matchedCellBlocks.begin(), matchedCellBlocks.end() );

  selectRegionCellBlocks( region, requestedAttributeValues, cellBlocksSelection );

  return cellBlocksSelection;
}

void CellElementRegionSelector::selectRegionCellBlocks( CellElementRegion const & region,
                                                        std::set< integer > const & attributeValues,
                                                        std::set< string > const & cellBlockNames )
{
  for( integer attributeValue : attributeValues )
  {
    string attributeValueStr = std::to_string( attributeValue );

    auto const [existingOwnerIt, isNewOwner] = m_regionAttributeOwner.emplace( attributeValueStr, &region );
    GEOS_THROW_IF( !isNewOwner && ( existingOwnerIt->second != &region ),
                   GEOS_FMT( "The region attribute '{}' has been referenced in multiple {}:\n- {}\n- {}",
                             attributeValueStr, CellElementRegion::catalogName(),
                             existingOwnerIt->second->getDataContext(), region.getDataContext() ),
                   InputError );

    // the regionAttribute is now selected by the region
    m_orphanRegionAttributes.erase( attributeValueStr );
  }

  for( string const & cellBlockName : cellBlockNames )
  {
    auto const [existingOwnerIt, isNewOwner] = m_cellBlocksOwner.emplace( cellBlockName, &region );
    GEOS_THROW_IF( !isNewOwner && ( existingOwnerIt->second != &region ),
                   GEOS_FMT( "The cellBlock '{}' has been referenced in multiple {}:\n- {}\n- {}",
                             cellBlockName, CellElementRegion::catalogName(),
                             existingOwnerIt->second->getDataContext(), region.getDataContext() ),
                   InputError );

    // the cellBlocks is now selected by the region
    m_orphanCellBlockNames.erase( cellBlockName );
  }
}


CellBlockABC const & CellElementRegionSelector::getCellBlock( CellElementRegion const & region,
                                                              string const & cellBlockName ) const
{
  CellBlockABC const * cellBlock = m_cellBlocks.getGroupPointer< CellBlockABC >( string( cellBlockName ) );

  // Error should never be called as every cell-block names should already be checked.
  GEOS_ERROR_IF( cellBlock == nullptr,
                 GEOS_FMT( "{}: Unexpected error, {} not found.\nAvailable cell-block list: {{ {} }}",
                           region.getWrapperDataContext( ViewKeys::sourceCellBlockNamesString() ),
                           cellBlockName, stringutilities::join( m_cellBlockNames, ", " ) ) );

  return *cellBlock;
}


void CellElementRegionSelector::checkForNoOrphanCellBlocks() const
{
  if( !m_orphanCellBlockNames.empty() )
  {
    std::ostringstream oss;
    if( !m_orphanRegionAttributes.empty() )
    {
      oss << GEOS_FMT( "The {} {{ {} }} has not been referenced in any region.\n",
                       ViewKeys::cellBlockAttributeValuesString(),
                       stringutilities::join( m_orphanRegionAttributes, ", " ));
    }
    oss << GEOS_FMT( "The following cell-blocks has not been referenced in any region: {{ {} }}.\n",
                     stringutilities::join( m_orphanCellBlockNames, ", " ));
    oss << GEOS_FMT( "Please add it in an existing {} (through {}, {} or {}), or consider creating a new one to describe your model.",
                     CellElementRegion::catalogName(),
                     ViewKeys::cellBlockAttributeValuesString(),
                     ViewKeys::sourceCellBlockNamesString(),
                     ViewKeys::cellBlockMatchPatternsString() );
    GEOS_THROW( oss.str(), InputError );
  }
}


} /* namespace geos */
