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

/**
 * @file PAMELAMeshGenerator.cpp
 */

#include "PAMELAMeshGenerator.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/TypeDispatch.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "CellBlockManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

// PAMELA includes
#include "Elements/Element.hpp"
#include "Mesh/MeshFactory.hpp"
#include "MeshDataWriters/Variable.hpp"
#include "MeshDataWriters/MeshParts.hpp"

#include <unordered_set>

namespace geosx
{
using namespace dataRepository;

PAMELAMeshGenerator::PAMELAMeshGenerator( string const & name, Group * const parent ):
  ExternalMeshGeneratorBase( name, parent )
{}

namespace
{

string makeRegionLabel( string const & regionName, string const & regionCellType )
{
  return regionName + "_" + regionCellType;
}

/*!
 * @brief Given the PAMELA Surface or Region label, return a simple unique name for GEOSX
 * @details Surface labels in PAMELA are composed of different parts such as the index, the type of cells, etc.
 * @param[in] pamelaLabel the surface or region label within PAMELA
 * @return the name of the surface or the region
 */
string retrieveSurfaceOrRegionName( string const & pamelaLabel )
{
  // The PAMELA label looks like: PART00002_POLYGON_POLYGON_GROUP_Ovbd1_Ovbd2_14
  // (Ovbd1_Ovbd2 may or may not be present), we keep just the part following GROUP
  // This can be done without splitting, but we want to keep flexibility to change.
  array1d< string > const splitLabel = stringutilities::tokenize( pamelaLabel, "_" );
  auto const it = std::find( splitLabel.begin(), splitLabel.end(), "GROUP" );
  GEOSX_THROW_IF( it == splitLabel.end(),
                  "GEOSX assumes that PAMELA places the word GROUP before the region/surface name",
                  InputError );
  // Handle both cases where region name is present and missing (keep the number in that case)
  auto const first = it + 1;
  auto const last = splitLabel.end() - 1;
  return stringutilities::join( first, (first == last) ? last + 1 : last, "_" );
}

string const & getElementLabel( PAMELA::ELEMENTS::TYPE const type )
{
  static std::map< PAMELA::ELEMENTS::TYPE, string > const typeToLabelMap =
  {
    { PAMELA::ELEMENTS::TYPE::VTK_VERTEX, "VERTEX" },
    { PAMELA::ELEMENTS::TYPE::VTK_LINE, "LINE" },
    { PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE, "TRIANGLE" },
    { PAMELA::ELEMENTS::TYPE::VTK_QUAD, "QUAD" },
    { PAMELA::ELEMENTS::TYPE::VTK_TETRA, "TETRA" },
    { PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON, "HEX" },
    { PAMELA::ELEMENTS::TYPE::VTK_WEDGE, "WEDGE" },
    { PAMELA::ELEMENTS::TYPE::VTK_PYRAMID, "PYRAMID" }
  };
  GEOSX_THROW_IF( typeToLabelMap.count( type ) == 0, "Unsupported PAMELA element type", std::runtime_error );
  return typeToLabelMap.at( type );
}

ElementType toGeosxElementType( PAMELA::ELEMENTS::TYPE const type )
{
  switch( type )
  {
    case PAMELA::ELEMENTS::TYPE::VTK_LINE: return ElementType::Line;
    case PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE: return ElementType::Triangle;
    case PAMELA::ELEMENTS::TYPE::VTK_QUAD: return ElementType::Quadrilateral;
    case PAMELA::ELEMENTS::TYPE::VTK_TETRA: return ElementType::Tetrahedron;
    case PAMELA::ELEMENTS::TYPE::VTK_PYRAMID: return ElementType::Pyramid;
    case PAMELA::ELEMENTS::TYPE::VTK_WEDGE: return ElementType::Wedge;
    case PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON: return ElementType::Hexahedron;
    default:
    {
      GEOSX_THROW( "Unsupported PAMELA element type", std::runtime_error );
    }
  }
}

PAMELA::ELEMENTS::TYPE toPamelaElementType( ElementType const type )
{
  switch( type )
  {
    case ElementType::Line: return PAMELA::ELEMENTS::TYPE::VTK_LINE;
    case ElementType::Triangle: return PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE;
    case ElementType::Quadrilateral: return PAMELA::ELEMENTS::TYPE::VTK_QUAD;
    case ElementType::Tetrahedron: return PAMELA::ELEMENTS::TYPE::VTK_TETRA;
    case ElementType::Pyramid: return PAMELA::ELEMENTS::TYPE::VTK_PYRAMID;
    case ElementType::Wedge: return PAMELA::ELEMENTS::TYPE::VTK_WEDGE;
    case ElementType::Hexahedron: return PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON;
    default:
    {
      GEOSX_THROW( "Unsupported PAMELA element type", std::runtime_error );
    }
  }
}

std::vector< int > getPamelaNodeOrder( PAMELA::ELEMENTS::TYPE const type,
                                       bool const isZReverse )
{
  switch( type )
  {
    case PAMELA::ELEMENTS::TYPE::VTK_LINE: return { 0, 1 };
    case PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE: return { 0, 1, 2 };
    case PAMELA::ELEMENTS::TYPE::VTK_QUAD: return { };
    case PAMELA::ELEMENTS::TYPE::VTK_TETRA: return { 0, 1, 2, 3 };
    case PAMELA::ELEMENTS::TYPE::VTK_PYRAMID: return { 0, 1, 3, 2, 4 };
    case PAMELA::ELEMENTS::TYPE::VTK_WEDGE: return { 0, 3, 1, 4, 2, 5 };
    case PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON:
    {
      // if the reverseZ option is on, we have to switch the node ordering
      // (top nodes become bottom nodes) to make sure that the hex volume is positive
      if( isZReverse )
      {
        return { 4, 5, 7, 6, 0, 1, 3, 2 };
      }
      else
      {
        return { 0, 1, 3, 2, 4, 5, 7, 6 };
      }
    }
    default:
    {
      GEOSX_THROW( "Unsupported PAMELA element type", std::runtime_error );
    }
  }
}

/// @return mesh length scale
real64 importNodes( PAMELA::Mesh & srcMesh, // PAMELA is not const-correct,
                    R1Tensor const & scale,
                    R1Tensor const & translate,
                    CellBlockManager & cellBlockManager )
{
  cellBlockManager.setNumNodes( srcMesh.get_PointCollection()->size_all() );
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > X = cellBlockManager.getNodePositions();

  arrayView1d< globalIndex > const nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();

  auto & nodeSets = cellBlockManager.getNodeSets();
  SortedArray< localIndex > & allNodes = nodeSets["all"];

  constexpr real64 minReal = LvArray::NumericLimits< real64 >::min;
  constexpr real64 maxReal = LvArray::NumericLimits< real64 >::max;
  real64 xMin[3] = { maxReal, maxReal, maxReal };
  real64 xMax[3] = { minReal, minReal, minReal };

  for( auto const & verticesIterator : *srcMesh.get_PointCollection() )
  {
    localIndex const vertexLocalIndex = verticesIterator->get_localIndex();
    globalIndex const vertexGlobalIndex = verticesIterator->get_globalIndex();
    X( vertexLocalIndex, 0 ) = ( verticesIterator->get_coordinates().x + translate[0] ) * scale[0];
    X( vertexLocalIndex, 1 ) = ( verticesIterator->get_coordinates().y + translate[1] ) * scale[1];
    X( vertexLocalIndex, 2 ) = ( verticesIterator->get_coordinates().z + translate[2] ) * scale[2];
    allNodes.insert( vertexLocalIndex );
    nodeLocalToGlobal[vertexLocalIndex] = vertexGlobalIndex;
    for( int i = 0; i < 3; ++i )
    {
      xMin[i] = std::min( xMin[i], X( vertexLocalIndex, i ) );
      xMax[i] = std::max( xMax[i], X( vertexLocalIndex, i ) );
    }
  }

  MpiWrapper::allReduce( xMin, xMin, 3, MPI_MIN, MPI_COMM_GEOSX );
  MpiWrapper::allReduce( xMax, xMax, 3, MPI_MAX, MPI_COMM_GEOSX );
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

void importCellBlock( PAMELA::SubPart< PAMELA::Polyhedron * > * const cellBlockPtr,
                      string const & cellBlockName,
                      bool const isZReverse,
                      CellBlockManager & cellBlockManager )
{
  GEOSX_ASSERT( cellBlockPtr != nullptr );

  CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );
  cellBlock.setElementType( toGeosxElementType( cellBlockPtr->ElementType ) );
  cellBlock.resize( LvArray::integerConversion< localIndex >( cellBlockPtr->SubCollection.size_owned() ) );

  std::vector< int > const nodeOrder = getPamelaNodeOrder( cellBlockPtr->ElementType, isZReverse );
  arrayView2d< localIndex, cells::NODE_MAP_USD > const cellToVertex = cellBlock.getElemToNode().toView();
  arrayView1d< globalIndex > const localToGlobal = cellBlock.localToGlobalMap();

  // Iterate on cells
  auto & subCollection = cellBlockPtr->SubCollection;
  for( auto cellItr = subCollection.begin_owned(); cellItr != subCollection.end_owned(); ++cellItr )
  {
    PAMELA::Polyhedron const * const cellPtr = *cellItr;
    GEOSX_ASSERT( cellPtr != nullptr );
    localIndex const cellLocalIndex = cellPtr->get_localIndex();
    globalIndex const cellGlobalIndex = cellPtr->get_globalIndex();

    std::vector< PAMELA::Point * > const & cornerList = cellPtr->get_vertexList();
    GEOSX_ASSERT_EQ( cornerList.size(), nodeOrder.size() );

    for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( cornerList.size() ); ++i )
    {
      cellToVertex[cellLocalIndex][i] = cornerList[nodeOrder[i]]->get_localIndex();
    }
    localToGlobal[cellLocalIndex] = cellGlobalIndex;
  }
}

void importSurface( PAMELA::Part< PAMELA::Polygon * > * const surfacePtr,
                    string const & surfaceName,
                    CellBlockManager & cellBlockManager )
{
  GEOSX_ASSERT( surfacePtr != nullptr );
  SortedArray< localIndex > & curNodeSet = cellBlockManager.getNodeSets()[ surfaceName ];

  for( auto const & subPart : surfacePtr->SubParts )
  {
    PAMELA::SubPart< PAMELA::Polygon * > * const cellBlockPtr = subPart.second;
    GEOSX_ASSERT( cellBlockPtr != nullptr );
    if( cellBlockPtr->ElementType == PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE ||
        cellBlockPtr->ElementType == PAMELA::ELEMENTS::TYPE::VTK_QUAD )
    {
      auto & subCollection = cellBlockPtr->SubCollection;
      for( auto cellItr = subCollection.begin_owned(); cellItr != subCollection.end_owned(); cellItr++ )
      {
        PAMELA::Polygon const * const cellPtr = *cellItr;
        GEOSX_ASSERT( cellPtr != nullptr );
        std::vector< PAMELA::Point * > const & cornerList = cellPtr->get_vertexList();
        for( auto corner : cornerList )
        {
          curNodeSet.insert( corner->get_localIndex() );
        }
      }
    }
  }
}

} // namespace

void PAMELAMeshGenerator::generateMesh( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
  m_pamelaMesh = std::unique_ptr< PAMELA::Mesh >( PAMELA::MeshFactory::makeMesh( m_filePath ) );
  m_pamelaMesh->CreateFacesFromCells();
  m_pamelaMesh->PerformPolyhedronPartitioning( PAMELA::ELEMENTS::FAMILY::POLYGON,
                                               PAMELA::ELEMENTS::FAMILY::POLYGON );
  m_pamelaMesh->CreateLineGroupWithAdjacency( "TopologicalC2C",
                                              m_pamelaMesh->getAdjacencySet()->get_TopologicalAdjacency( PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
                                                                                                         PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
                                                                                                         PAMELA::ELEMENTS::FAMILY::POLYGON ) );

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': generating GEOSX mesh data structure", catalogName(), getName() ) );

  domain.getMetisNeighborList() = m_pamelaMesh->getNeighborList();
  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
  //TODO for the moment we only consider on mesh level "Level0"
  meshBody.createMeshLevel( 0 );

  CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

  // Dealing with nodes
  real64 const lengthScale = importNodes( *m_pamelaMesh, m_scale, m_translate, cellBlockManager );
  meshBody.setGlobalLengthScale( lengthScale );

  // Use the PartMap of PAMELA to get the mesh
  PAMELA::PartMap< PAMELA::Polyhedron * > const polyhedronPartMap =
    std::get< 0 >( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get(), 0 ) );

  // Import cell blocks
  for( auto const & polyhedronPart : polyhedronPartMap )
  {
    PAMELA::Part< PAMELA::Polyhedron * > * const regionPtr = polyhedronPart.second;
    GEOSX_ASSERT( regionPtr != nullptr );
    string const regionName = retrieveSurfaceOrRegionName( regionPtr->Label );

    // Iterate on cell types
    for( auto const & subPart : regionPtr->SubParts )
    {
      // Ignore non-polyhedrons elements (somehow they exist within a polyhedron Part in PAMELA data model).
      auto const elementFamily = PAMELA::ELEMENTS::TypeToFamily.at( static_cast< int >( subPart.second->ElementType ) );
      if( elementFamily == PAMELA::ELEMENTS::FAMILY::POLYHEDRON )
      {
        size_t const numGlobalCells = MpiWrapper::sum( subPart.second->SubCollection.size_owned() );
        if( numGlobalCells > 0 )
        {
          string const cellBlockName = makeRegionLabel( regionName, getElementLabel( subPart.second->ElementType ) );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Importing cell block " << cellBlockName );
          importCellBlock( subPart.second, cellBlockName, m_scale[2] < 0, cellBlockManager );
          m_cellBlockRegions.emplace( cellBlockName, polyhedronPart.first );
        }
      }
    }
  }

  // Import surfaces
  PAMELA::PartMap< PAMELA::Polygon * > const polygonPartMap = std::get< 0 >( PAMELA::getPolygonPartMap( m_pamelaMesh.get(), 0 ));
  for( auto const & polygonPart : polygonPartMap )
  {
    PAMELA::Part< PAMELA::Polygon * > * const surfacePtr = polygonPart.second;
    string const surfaceName = retrieveSurfaceOrRegionName( surfacePtr->Label );
    GEOSX_LOG_LEVEL_RANK_0( 1, "Importing surface " << surfaceName );
    importSurface( surfacePtr, surfaceName, cellBlockManager );
  }

  cellBlockManager.buildMaps();
}

void PAMELAMeshGenerator::freeResources()
{
  m_pamelaMesh.reset();
  m_cellBlockRegions.clear();
}

namespace
{

void importRegularField( PAMELA::VariableDouble & source,
                         std::vector< int > const & indexMap,
                         WrapperBase & wrapper )
{
  // Scalar material fields are stored as 1D arrays, vector/tensor are 2D
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 1, 2 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( source.offset );
    localIndex const numComponentsDst = wrapperT.numArrayComp();
    GEOSX_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                           "Mismatch in number of components for field " << source.Label );

    // the assumption here is that before this step, GEOSX has added ghost elements **at the end** of the view
    localIndex const numCells = LvArray::integerConversion< localIndex >( indexMap.size() );
    for( localIndex i = 0; i < numCells; ++i )
    {
      // get_data() currently returns a new vector for each cell, but auto const & makes sure
      // this code keeps working if/when PAMELA is fixed to return a more reasonable type (span/pointer)
      auto const & data = source.get_data( indexMap[i] );
      int comp = 0;
      LvArray::forValuesInSlice( view[i], [&]( auto & val )
      {
        val = data[comp++];
      } );
    }
  } );
}

void importMaterialField( PAMELA::VariableDouble & source,
                          std::vector< int > const & indexMap,
                          WrapperBase & wrapper )
{
  // Scalar material fields are stored as 2D arrays, vector/tensor are 3D
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 2, 3 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( source.offset );
    localIndex const numComponentsDst = wrapperT.numArrayComp() / view.size( 1 );
    GEOSX_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                           "Mismatch in number of components for field " << source.Label );

    // the assumption here is that before this step, GEOSX has added ghost elements **at the end** of the view
    localIndex const numCells = LvArray::integerConversion< localIndex >( indexMap.size() );
    for( localIndex i = 0; i < numCells; ++i )
    {
      // get_data() currently returns a new vector for each cell, but auto const & makes sure
      // this code keeps working if/when PAMELA is fixed to return a more reasonable type (span/pointer)
      auto const & data = source.get_data( indexMap[i] );

      // For material fields, copy identical values into each quadrature point
      for( int q = 0; q < view.size( 1 ); ++q )
      {
        int comp = 0;
        LvArray::forValuesInSlice( view[i][q], [&]( auto & val )
        {
          val = data[comp++];
        } );
      }
    }
  } );
}

std::unordered_set< string > getMaterialWrapperNames( ElementSubRegionBase const & subRegion )
{
  using namespace constitutive;
  std::unordered_set< string > materialWrapperNames;
  subRegion.getConstitutiveModels().forSubGroups< ConstitutiveBase >( [&]( ConstitutiveBase const & material )
  {
    material.forWrappers( [&]( WrapperBase const & wrapper )
    {
      if( wrapper.sizedFromParent() )
      {
        materialWrapperNames.insert( ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() ) );
      }
    } );
  } );
  return materialWrapperNames;
}

} // namespace

void PAMELAMeshGenerator::importFieldsOnSubRegion( PAMELA::Part< PAMELA::Polyhedron * > & srcRegion,
                                                   PAMELA::SubPart< PAMELA::Polyhedron * > & srcSubRegion,
                                                   ElementRegionBase const & dstRegion,
                                                   ElementSubRegionBase & dstSubRegion ) const
{
  // Make a list of material field wrapper names on current subregion
  std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( dstSubRegion );

  for( localIndex fieldIndex = 0; fieldIndex < m_fieldsToImport.size(); fieldIndex++ )
  {
    string const & sourceName = m_fieldsToImport[fieldIndex];
    string const & wrapperName = m_fieldNamesInGEOSX[fieldIndex];

    // Find destination
    if( !dstSubRegion.hasWrapper( wrapperName ) )
    {
      // Skip - the user may have not enabled a particular physics model/solver on this dstRegion.
      GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Skipping import of {} -> {} on {}/{} (field not found)",
                                            sourceName, wrapperName, dstRegion.getName(), dstSubRegion.getName() ) );
      continue;
    }
    WrapperBase & wrapper = dstSubRegion.getWrapperBase( wrapperName );

    GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Importing field {} -> {} on {}/{}",
                                          sourceName, wrapperName, dstRegion.getName(), dstSubRegion.getName() ) );

    // Find source - and hard error if it doesn't exist
    PAMELA::VariableDouble * const meshProperty = srcRegion.FindVariableByName( sourceName );
    GEOSX_THROW_IF( meshProperty == nullptr,
                    GEOSX_FMT( "Field not found in source dataset: {}", sourceName ),
                    InputError );

    // Decide if field is constitutive or not
    if( materialWrapperNames.count( wrapperName ) > 0 && wrapper.numArrayDims() > 1 )
    {
      importMaterialField( *meshProperty, srcSubRegion.IndexMapping, wrapper );
    }
    else
    {
      importRegularField( *meshProperty, srcSubRegion.IndexMapping, wrapper );
    }
  }
}

void PAMELAMeshGenerator::importFields( DomainPartition & domain ) const
{
  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': importing field data from mesh dataset", catalogName(), getName() ) );
  GEOSX_ASSERT_MSG( m_pamelaMesh, "Must call generateMesh() before importFields()" );

  ElementRegionManager & elemManager = domain.getMeshBody( this->getName() ).getBaseDiscretization().getElemManager();

  PAMELA::PartMap< PAMELA::Polyhedron * > const polyhedronPartMap =
    std::get< 0 >( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get(), 0 ) );

  std::vector< string > regionNames;

  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex,
                                                                         localIndex,
                                                                         ElementRegionBase const & region,
                                                                         CellElementSubRegion & subRegion )
  {
    // Use stored map of names to get access to PAMELA region again
    GEOSX_ERROR_IF_NE( m_cellBlockRegions.count( subRegion.getName() ), 1 );
    string const & srcRegionName = m_cellBlockRegions.at( subRegion.getName() );
    PAMELA::Part< PAMELA::Polyhedron * > * const regionPtr = polyhedronPartMap.at( srcRegionName );
    GEOSX_ERROR_IF( regionPtr == nullptr, "Internal logic error (PAMELA region not found)" );

    // Use element type to get the PAMELA cell block (sub-part) corresponding to subregion
    PAMELA::ELEMENTS::TYPE const elemType = toPamelaElementType( subRegion.getElementType() );
    PAMELA::SubPart< PAMELA::Polyhedron * > * const cellBlockPtr = regionPtr->SubParts.at( static_cast< int >( elemType ) );
    GEOSX_ERROR_IF( cellBlockPtr == nullptr, "Internal logic error (PAMELA subregion not found)" );

    importFieldsOnSubRegion( *regionPtr, *cellBlockPtr, region, subRegion );
    regionNames.emplace_back( region.getName() );
  } );

  std::vector< string > tmp( m_fieldNamesInGEOSX.begin(), m_fieldNamesInGEOSX.end() );
  FieldIdentifiers fieldsToBeSync;
  fieldsToBeSync.addElementFields( tmp, regionNames );

  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       domain.getMeshBody( this->getName() ).getBaseDiscretization(),
                                                       domain.getNeighbors(),
                                                       false );
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, PAMELAMeshGenerator, string const &, Group * const )
}
