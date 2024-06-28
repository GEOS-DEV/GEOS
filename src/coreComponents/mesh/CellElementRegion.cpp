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
    setDescription( GEOS_FMT( "List of the desired cell-blocks names to contain in this {}.\n"
                              "The form of this attribute is of \"regionAttribute_elementType\", so \"1_tetrahedra\" select the "
                              "cell-block that contains the tetrahedric elements for which the regionAttribute is 1.\n"
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
    setDescription( GEOS_FMT( "List of fnmatch pattern to match cell-block names to add them in this {}.\n"
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


void CellElementRegion::generateMesh( Group const & cellBlocks )
{
  Group & subRegions = this->getGroup( viewKeyStruct::elementSubRegions() );
  for( string const & cellBlockName : m_cellBlockNames )
  {
    CellBlockABC const * cellBlock = cellBlocks.getGroupPointer< CellBlockABC >( cellBlockName );
    GEOS_THROW_IF( cellBlock == nullptr,
                   GEOS_FMT( "{}: No cellBlock named '{}'.\nAvailable cellBlock list: {{ {} }}\nNo CellElementRegionSelector has been used to verify the cellBlock selection.",
                             getWrapperDataContext( viewKeyStruct::sourceCellBlockNamesString() ),
                             cellBlockName,
                             stringutilities::join( m_cellBlockNames, ", " ) ),
                   InputError );

    // subRegion name must be the same as the cell-block (so we can match them and reference them in errors).
    CellElementSubRegion & subRegion = subRegions.registerGroup< CellElementSubRegion >( cellBlockName );
    subRegion.copyFromCellBlock( *cellBlock );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geos */
