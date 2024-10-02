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


// std::optional< string > CellElementRegionSelector::getCellBlockAttributeValue( string_view cellBlockName )
// {
//   if( cellBlockName.find( cellBlockTypeSeparator ) != string_view::npos )
//   {
//     for( )
//     {
//       if( stringutilities::endsWith( cellBlockName, primitiveSuffix ))
//       {
//         return stringutilities::removeStringAndFollowingContent( cellBlockName, cellBlockTypeSeparator );
//       }
//     }
//   }
//   return {};
// }

// bool CellElementRegionSelector::isRegionCellBlock( string_view cellBlockName )
// {
//   return cellBlockName.find( cellBlockTypeSeparator ) != string_view::npos;
// }


CellElementRegionSelector::CellElementRegionSelector( Group const & cellBlocks,
                                                      RegionAttributesCellBlocksMap const & regionsCellBlocks ):
  m_regionAttributesCellBlocks( regionsCellBlocks )
{
  cellBlocks.forSubGroups< CellBlockABC >( [&] ( CellBlockABC const & cellBlock )
  {
    string const name = cellBlock.getName();
    m_cellBlocksOwners.emplace( name, std::vector< CellElementRegion const * >() );//possible de supprimer en gardant un set des noms de
                                                                                   // cellblocks ?
  } );

  for( auto const & regionAttribute : regionsCellBlocks )
  {
    m_regionAttributesOwners.emplace( regionAttribute.first, std::vector< CellElementRegion const * >() );//possible de supprimer
  }
}


// std::set< string >
// CellElementRegionSelector::buildMatchPatterns( CellElementRegion const & region,
//                                                std::set< integer > const & requestedAttributeValues,
//                                                std::set< string > const & requestedMatchPatterns ) const
// {
//   std::set< string > matchPatterns;

//   // user region attribute matching patterns creation
//   for( integer attributeValue : requestedAttributeValues )
//   {
//     string attributeValueStr = std::to_string( attributeValue );
//     // if attribute value does not exist in the mesh
//     GEOS_THROW_IF( m_regionAttributesOwners.count( attributeValueStr ) == 0,
//                    GEOS_FMT( "{}: Region attribute value '{}' not found.\nAvailable region attribute list: {{ {} }}",
//                              region.getWrapperDataContext( ViewKeys::cellBlockAttributeValuesString() ),
//                              attributeValueStr,
//                              stringutilities::joinLamda( m_regionAttributesOwners, ", ",
//                                                          []( auto pair ) { return pair->first; } ) ),
//                    InputError );

//     // for each desired attribute values, we add the following the match patterns:
//     // Attribute value, followed by an underscore, followed by one or more characters.
//     matchPatterns.emplace( GEOS_FMT( "{}_?*", attributeValueStr ) );
//   }

//   // user fnMatch patterns
//   for( string const & matchPattern : requestedMatchPatterns )
//   {
//     matchPatterns.emplace( matchPattern );
//   }

//   return matchPatterns;
// }

// std::set< string >
// CellElementRegionSelector::getMatchingCellblocks( CellElementRegion const & region,
//                                                   std::set< string > const & matchPatterns ) const
// {
//   std::set< string > matchedCellBlocks;
//   for( string const & matchPattern : matchPatterns )
//   {
//     bool matching = false;
//     for( auto const & [cellBlockName, owners] : m_cellBlocksOwners )
//     {
//       // if the pattern matches the tested cellBlock name
//       if( fnmatch( matchPattern.c_str(), cellBlockName.c_str(), 0 ) == 0 )
//       {
//         matching = true;
//         matchedCellBlocks.emplace( cellBlockName );
//       }
//     }

//     GEOS_THROW_IF( !matching,
//                    GEOS_FMT( "{}: No cellBlock name is satisfying the pattern '{}'.\nAvailable cellBlock list: {{ {} }}",
//                              region.getWrapperDataContext( ViewKeys::cellBlockMatchPatternsString() ),
//                              matchPattern,
//                              stringutilities::joinLamda( m_cellBlocksOwners, ", ",
//                                                          []( auto pair ) { return pair->first; } ) ),
//                    InputError );
//   }
//   return matchedCellBlocks;
// }
std::set< string >
CellElementRegionSelector::getMatchingCellblocks( CellElementRegion const & region,
                                                  string_view matchPattern ) const
{
  std::set< string > matchedCellBlocks;
  bool matching = false;
  for( auto const & [cellBlockName, owners] : m_cellBlocksOwners )
  {
    // if the pattern matches the tested cellBlock name
    if( fnmatch( matchPattern.data(), cellBlockName.c_str(), 0 ) == 0 )
    {
      matching = true;
      matchedCellBlocks.emplace( cellBlockName );
    }
  }

  GEOS_THROW_IF( !matching,
                 GEOS_FMT( "{}: No cellBlock name is satisfying the qualifier '{}'.\n"
                           "Available cellBlock list: {{ {} }}\nAvailable region attribute list: {{ {} }}",
                           region.getWrapperDataContext( ViewKeys::sourceCellBlockQualifiersString() ),
                           matchPattern,
                           stringutilities::joinLamda( m_regionAttributesOwners, ", ",
                                                       []( auto pair ) { return pair->first; } ),
                           stringutilities::joinLamda( m_cellBlocksOwners, ", ",
                                                       []( auto pair ) { return pair->first; } ) ),
                 InputError );
  return matchedCellBlocks;
}

void
CellElementRegionSelector::verifyRequestedCellBlocks( CellElementRegion const & region,
                                                      std::set< string > const & cellBlockNames ) const
{
  for( string const & requestedCellBlockName : cellBlockNames )
  {
    // if cell block does not exist in the mesh
    GEOS_THROW_IF( m_cellBlocksOwners.count( requestedCellBlockName ) == 0,
                   GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}",
                             region.getWrapperDataContext( ViewKeys::sourceCellBlockQualifiersString() ),
                             requestedCellBlockName,
                             stringutilities::joinLamda( m_cellBlocksOwners, ", ",
                                                         []( auto pair ) { return pair->first; } ) ),
                   InputError );
  }
}


void
CellElementRegionSelector::registerRegionSelection( CellElementRegion const & region,
                                                    std::set< string > const & cellBlockNames,
                                                    std::set< string > const & attributeValues )
{
  for( string const & attributeValue : attributeValues )
  {
    m_regionAttributesOwners[attributeValue].push_back( &region );
    GEOS_LOG( "ajout de attributeValue="<<attributeValue<<" dans la region "<<region.getName());
  }

  for( string const & cellBlockName : cellBlockNames )
  {
    m_cellBlocksOwners[cellBlockName].push_back( &region );
    GEOS_LOG( "ajout de cellBlockName="<<cellBlockName<<" dans la region "<<region.getName());
  }
}


std::set< string > CellElementRegionSelector::buildCellBlocksSelection( CellElementRegion const & region )
{
  std::set< string > cellBlocksSelection;
  std::set< string > regionAttributeSelection;

  auto const qualifiers = region.getCellBlockQualifiers();
  for( string const & qualifier : qualifiers )
  {
    auto const regionCellBlocks = m_regionAttributesCellBlocks.find( qualifier );
    if( regionCellBlocks != m_regionAttributesCellBlocks.end() )
    { // if the qualifier is a region attribute value, let's select it
      GEOS_LOG(GEOS_FMT("selection du regionAttribute {}", qualifier));
      regionAttributeSelection.emplace( regionCellBlocks->first );
      for( string const & cellBlock : regionCellBlocks->second )
      {
        cellBlocksSelection.emplace( cellBlock );
      }
    }
    else
    { // the qualifier is a match pattern, or a simple cellblock name, let's select all matching cell-blocks
      GEOS_LOG(GEOS_FMT("selection du matchPattern {}", qualifier));
      std::set< string > const matchedCellBlocks = getMatchingCellblocks( region, qualifier );
      cellBlocksSelection.insert( matchedCellBlocks.begin(), matchedCellBlocks.end() );
    }
  }

  verifyRequestedCellBlocks( region, cellBlocksSelection );
  registerRegionSelection( region, cellBlocksSelection, regionAttributeSelection );

  return cellBlocksSelection;
}

void CellElementRegionSelector::checkSelectionConsistency() const
{
  auto const getRegionStr = []( auto regionPtrIterator ) -> string {
    return GEOS_FMT( "- {}", (*regionPtrIterator)->getDataContext() );
  };

  auto const checkOwnerCount = [&]( string_view qualifierType,
                                    auto const & qualifiersOwners,
                                    auto & orphanList ) {
    // Search of never or multiple selected attribute values
    std::vector< string > multipleRefsErrors;
    for( auto const & [qualifier, owningRegions] : qualifiersOwners )
    {
      if( owningRegions.size() == 0 )
      {
        orphanList.insert( qualifier );
      }
      else if( owningRegions.size() > 1 )
      {
        GEOS_LOG(GEOS_FMT("WOW n°{} : {}",multipleRefsErrors.size(),qualifier));
        multipleRefsErrors.push_back(
          GEOS_FMT( "The {} '{}' has been referenced in multiple {}:\n{}",
                    qualifierType, qualifier, CellElementRegion::catalogName(),
                    stringutilities::joinLamda( owningRegions, '\n', getRegionStr ) ) );
      }
    }
    GEOS_THROW_IF( !multipleRefsErrors.empty(), stringutilities::join( multipleRefsErrors, "\n\n" ), InputError );
  };

  std::set< string > orphanRegionAttributes;
  std::set< string > orphanCellBlockNames;
  checkOwnerCount( "region attribute", m_regionAttributesOwners, orphanRegionAttributes );
  checkOwnerCount( "cell-block", m_cellBlocksOwners, orphanCellBlockNames );

  if( !orphanCellBlockNames.empty() )
  {
    std::ostringstream oss;
    if( !orphanRegionAttributes.empty() )
    {
      oss << GEOS_FMT( "The region attributes {{ {} }} has not been referenced in any {}.\n",
                       stringutilities::join( orphanRegionAttributes, ", " ),
                       CellElementRegion::catalogName() );
    }
    oss << GEOS_FMT( "The following cell-blocks has not been referenced in any region: {{ {} }}.\n",
                     stringutilities::join( orphanCellBlockNames, ", " ) );
    oss << GEOS_FMT( "Please add it in an existing {} (through the '{}' attribute), or consider creating a new one to describe your model.",
                     CellElementRegion::catalogName(),
                     ViewKeys::sourceCellBlockQualifiersString() );
    GEOS_THROW( oss.str(), InputError );
  }
}


} /* namespace geos */
