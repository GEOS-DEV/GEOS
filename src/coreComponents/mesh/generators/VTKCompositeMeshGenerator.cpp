/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VTKCompositeMeshGenerator.cpp
 */

#include "VTKCompositeMeshGenerator.hpp"

#include "mesh/generators/VTKFaceBlockUtilities.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
{
using namespace dataRepository;


VTKCompositeMeshGenerator::VTKCompositeMeshGenerator( string const & name,
                                                      Group * const parent )
  : ExternalMeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::regionAttributeString(), &m_attributeName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Name of the VTK cell attribute to use as region marker" );

  registerWrapper( viewKeyStruct::nodesetNamesString(), &m_nodesetNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the VTK nodesets to import" );

  registerWrapper( viewKeyStruct::mainBlockNameString(), &m_mainBlockName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "main" ).
    setDescription( "For multi-block files, name of the 3d mesh block." );

  registerWrapper( viewKeyStruct::faceBlockNamesString(), &m_faceBlockNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "For multi-block files, names of the face mesh block." );

  registerWrapper( viewKeyStruct::useGlobalIdsString(), &m_useGlobalIds ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Controls the use of global IDs in the input file for cells and points."
                    " If set to 0 (default value), the GlobalId arrays in the input mesh are used if available, and generated otherwise."
                    " If set to a negative value, the GlobalId arrays in the input mesh are not used, and generated global Ids are automatically generated."
                    " If set to a positive value, the GlobalId arrays in the input mesh are used and required, and the simulation aborts if they are not available" );
}

void VTKCompositeMeshGenerator::fillCellBlockManager( CellBlockManager & GEOS_UNUSED_PARAM( cellBlockManager ), array1d< int > const & )
{}



void VTKCompositeMeshGenerator::importFieldOnArray( Block GEOS_UNUSED_PARAM( block ),
                                                    string const & GEOS_UNUSED_PARAM( blockName ),
                                                    string const & GEOS_UNUSED_PARAM( meshFieldName ),
                                                    bool GEOS_UNUSED_PARAM( isMaterialField ),
                                                    dataRepository::WrapperBase & GEOS_UNUSED_PARAM( wrapper ) ) const
{}

void VTKCompositeMeshGenerator::freeResources()
{
  m_cellMap.clear();
  m_faceBlockMeshes.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKCompositeMeshGenerator, string const &, Group * const )

} // namespace geos
