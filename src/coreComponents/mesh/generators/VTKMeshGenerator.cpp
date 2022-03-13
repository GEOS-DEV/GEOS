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
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "CellBlockManager.hpp"
#include "HexCellBlockManager.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"

#include "common/MpiWrapper.hpp"
#include "common/TypeDispatch.hpp"

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"

#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#ifdef GEOSX_USE_MPI

#include <vtkMPIController.h>
#include <vtkMPI.h>

#else
  #include <vtkDummyController.h>
#endif

#include <numeric>
#include <unordered_set>

namespace geosx
{
using namespace dataRepository;


VTKMeshGenerator::VTKMeshGenerator( string const & name,
                                    Group * const parent )
  : MeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the mesh file" );
}

void VTKMeshGenerator::postProcessInput()
{ }

Group * VTKMeshGenerator::createChild( string const &
                                       GEOSX_UNUSED_PARAM( childKey ),
                                       string const &
                                       GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

template< typename VTK_ARRAY >
VTK_ARRAY const * getDataArrayOptional( vtkFieldData * source,
                                        string const & name )
{
  // returns nullptr if either lookup or cast fail
  return VTK_ARRAY::FastDownCast( source->GetAbstractArray( name.c_str() ) );
}

template< typename VTK_ARRAY >
VTK_ARRAY const & getDataArray( vtkFieldData * source,
                                string const & name )
{
  VTK_ARRAY const * const array = getDataArrayOptional< VTK_ARRAY >( source, name );
  GEOSX_ERROR_IF( array == nullptr,
                  "VTK array '" << name << "' not found or has unexpected data type" );
  return *array;
}

/**
 * @brief Returns a string describing the VTK element.
 * @param[in] vtkCellType The VTK element type.
 * @return The name.
 * @warning This information will be visible in the input file... Consider refactoring with great care.
 */
string getCellName( VTKCellType vtkCellType )
{
  switch( vtkCellType )
  {
    case VTK_HEXAHEDRON:
      return "hexahedra";
    case VTK_TETRA:
      return "tetrahedra";
    case VTK_WEDGE:
      return "wedges";
    case VTK_PYRAMID:
      return "pyramids";
    default:
      GEOSX_ERROR( "VTK cell type \"" << vtkCellType << "\" is not a recognized cell type" );
      // Dummy return for compilation error.
      return "";
  }
}

/**
 * @brief Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getVTKController()
{
#ifdef GEOSX_USE_MPI
  vtkNew< vtkMPIController > controller;
  vtkMPICommunicatorOpaqueComm vtkGeosxComm( &MPI_COMM_GEOSX );
  vtkNew< vtkMPICommunicator > communicator;
  communicator->InitializeExternal( &vtkGeosxComm );
  controller->SetCommunicator( communicator );
#else
  vtkNew< vtkDummyController > controller;
#endif
  return controller;
}

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 */
vtkSmartPointer< vtkUnstructuredGrid > loadVTKMesh( Path const & filePath )
{
  string const extension = stringutilities::tokenize( filePath, "." ).back(); // TODO maybe code a method in Path to get the file extension?
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh;

  if( extension == "pvtu" )
  {
    vtkSmartPointer< vtkXMLPUnstructuredGridReader > vtkUgReader = vtkSmartPointer< vtkXMLPUnstructuredGridReader >::New();
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->UpdateInformation();
    int const numberOfPieces = vtkUgReader->GetNumberOfPieces();
    if( MpiWrapper::commSize() == 1 )
    {
      vtkUgReader->Update();
    }
    else if( MpiWrapper::commRank() < numberOfPieces )
    {
      vtkUgReader->UpdatePiece( MpiWrapper::commRank(), 2, 0 ); // TODO What's that numPieces = 2?
    }
    loadedMesh = vtkUgReader->GetOutput();
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
      auto read = [&]( auto vtkUgReader ) // auto can be multiple types in the same function
      {
        vtkUgReader->SetFileName( filePath.c_str() );
        vtkUgReader->Update();
        return vtkUgReader->GetOutput();
      };

      if( extension == "vtk" )
      {
        loadedMesh = read( vtkSmartPointer< vtkUnstructuredGridReader >::New() );
      }
      else if( extension == "vtu" )
      {
        loadedMesh = read( vtkSmartPointer< vtkXMLUnstructuredGridReader >::New() );
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    }
    else
    {
      loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
    }
  }

  return loadedMesh;
}

/**
 * @brief Redistribute the mesh among the available MPI ranks
 * @details this method will also generate global ids for points and cells in the VTK Mesh
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 * @param[out] cuts the bounding boxes used by the VTK partitioner needed to compute neighboring ranks
 */
vtkSmartPointer< vtkUnstructuredGrid > redistributeMesh( vtkUnstructuredGrid & loadedMesh,
                                                         std::vector< vtkBoundingBox > & cuts )
{
  // Redistribute data all over the available ranks
  vtkNew< vtkRedistributeDataSetFilter > rdsf;
  rdsf->SetInputDataObject( &loadedMesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->GenerateGlobalCellIdsOn();
  rdsf->Update();

  cuts = rdsf->GetCuts();

  // Generate global IDs for vertices and cells
  vtkNew< vtkGenerateGlobalIds > generator;
  generator->SetInputDataObject( vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject( 0 ) ) );
  generator->Update();
  vtkSmartPointer< vtkUnstructuredGrid > mesh = vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject( 0 ) );
  return mesh;
}

/**
 * @brief Copy the VTK mesh nodes into the nodeManager of GEOSX
 * @param[in] nodeManager the NodeManager of the domain in which the points will be copied.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @return the global length of the mesh (diagonal of the bounding box)
 */
double writeMeshNodes( CellBlockManagerBase & cellBlockManager,
                       vtkSmartPointer< vtkUnstructuredGrid > mesh )
{
  cellBlockManager.setNumNodes( mesh->GetNumberOfPoints() );
  arrayView1d< globalIndex > const & nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();


  std::cout<< " VTK nb points " << mesh->GetNumberOfPoints() << std::endl;

  // Writing the points
  GEOSX_ERROR_IF( mesh->GetNumberOfPoints() == 0, "Mesh is empty, aborting" );

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = cellBlockManager.getNodesPositions();
  std::map< string, SortedArray< localIndex > > & nodeSets = cellBlockManager.getNodeSets();
  SortedArray< localIndex > & allNodes = nodeSets[ "all" ];
  real64 xMax[3] = { std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max(), std::numeric_limits< real64 >::max(), std::numeric_limits< real64 >::max() };

  vtkIdTypeArray const & globalPointIdDataArray = getDataArray< vtkIdTypeArray >( mesh->GetPointData(), "GlobalPointIds" );

  for( vtkIdType v = 0; v < mesh->GetNumberOfPoints(); v++ )
  {
    double * point = mesh->GetPoint( v );
    nodeLocalToGlobal[v] = globalPointIdDataArray.GetValue( v );
    for( integer i = 0; i < 3; i++ )
    {
      X( v, i ) = point[i];
      xMax[i] = std::max( xMax[i], point[i] );
      xMin[i] = std::min( xMin[i], point[i] );
    }
    allNodes.insert( v );
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

/**
 * @brief Compute the potential rank neighbor list
 * @param[in] cuts the bounding boxes used by the VTK partitioner for all ranks
 * @details Compute the rank neighbors. The asssumption is that 2 ranks are neighbors if
 * the corresponding bounding boxes intersect.
 */
std::set< int > computeMPINeighborRanks( const std::vector< vtkBoundingBox > & cuts )
{
  unsigned int nbRanks = MpiWrapper::commSize();
  GEOSX_ERROR_IF( cuts.size() != nbRanks, "The number of VTK cuts does not match the number of ranks" );

  std::set< int > rankNeighbor;
  unsigned int currentRank = MpiWrapper::commRank();
  vtkBoundingBox currentBox = cuts[currentRank];

  for( unsigned int i = 0; i < cuts.size(); i++ )
  {
    if( i != currentRank && currentBox.Intersects( cuts[i] ) )
    {
      rankNeighbor.insert( i );
    }
  }

  return rankNeighbor;
}

/**
 * @brief Gathers all the data from all ranks, merge them, sort them, and remove duplicates.
 * @tparam T Type of the exchanged data.
 * @param data The data to be exchanged.
 * @return The merged data.
 * @note This function makes MPI calls.
 */
template< typename T >
std::vector< T > gatherData( std::vector< T > const & data )
{
  // Exchange the sizes of the data across all ranks.
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( data.size() ), dataSizes, MPI_COMM_GEOSX );
  // `totalDataSize` contains the total data size across all the MPI ranks.
  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );

  // Once the MPI exchange is done, `allData` will contain all the data of all the MPI ranks.
  // We want all ranks to get all the data. But each rank may have a different size of information.
  // Therefore, we use `allgatherv` that does not impose the same size across ranks like `allgather` does.
  std::vector< T > allData( totalDataSize );
  // `displacements` is the offset (relative to the receive buffer) to store the data for each rank.
  std::vector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( data.data(), data.size(), allData.data(), dataSizes.data(), displacements.data(), MPI_COMM_GEOSX );

  // Finalizing by sorting, removing duplicates and trimming the result vector at the proper size.
  std::sort( allData.begin(), allData.end() );
  auto newEnd = std::unique( allData.begin(), allData.end() );
  allData.erase( newEnd, allData.end() );

  return allData;
}

/**
 * @brief This method is used to preprocess the VTK mesh and count the number of cells, facets, regions and
 * surfaces over the current MPI rank
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[out] regionsHex map from region index to the hexahedron indexes
 * @param[out] regionsTetra map from region index to the tetra indexes
 * @param[out] regionsWedges map from region index to the wedges indexes
 * @param[out] regionsPyramids map from region index to the pyramids indexes
 * @return A map from the surface indices to the associated cell ids for the current rank.
 * The mapping contains the surface indices for the whole simulation across the MPI ranks,
 * even if there is no cell for current rank (then the list will be empty).
 */
std::map< int, std::vector< vtkIdType > > buildRegionToCellsAndFaces( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                                                                      std::map< int, std::vector< vtkIdType > > & regionsHex,
                                                                      std::map< int, std::vector< vtkIdType > > & regionsTetra,
                                                                      std::map< int, std::vector< vtkIdType > > & regionsWedges,
                                                                      std::map< int, std::vector< vtkIdType > > & regionsPyramids )
{
  std::map< int, std::vector< vtkIdType > > surfacesIdsToCellsIds;

  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( mesh->GetCellData(), "attribute" );

  auto countCell = [&]( localIndex c,
                        std::map< int, std::vector< vtkIdType > > & regionsElem )
  {
    if( attributeDataArray != nullptr )
    {
      int const region = attributeDataArray->GetValue( c );
      GEOSX_LOG_RANK_0_IF( region < 0, "Attribute value " << region << " not supported, please use value >=0 to describe regions" );
      regionsElem[region].push_back( c );
    }
    else
    {
      regionsElem[-1].push_back( c );
    }
  };


  for( vtkIdType c = 0; c < mesh->GetNumberOfCells(); c++ )
  {
    if( VTK_HEXAHEDRON == mesh->GetCellType( c ) )
    {
      countCell( c, regionsHex );
    }
    else if( VTK_TETRA == mesh->GetCellType( c ) )
    {
      countCell( c, regionsTetra );
    }
    else if( VTK_WEDGE == mesh->GetCellType( c ) )
    {
      countCell( c, regionsWedges );
    }
    else if( VTK_PYRAMID == mesh->GetCellType( c ) )
    {
      countCell( c, regionsPyramids );
    }
    else if( VTK_TRIANGLE == mesh->GetCellType( c ) || VTK_QUAD == mesh->GetCellType( c ) )
    {
      if( attributeDataArray != nullptr )
      {
        surfacesIdsToCellsIds[attributeDataArray->GetValue( c )].push_back( c );
      }
    }
  }

  // TODO Output logger information on the loaded external mesh similarly to PAMELA

  // Communicate all the cells blocks
  std::size_t const numRegions = regionsHex.size() + regionsTetra.size() + regionsPyramids.size() + regionsWedges.size();
  std::vector< int > cellBlockRegionIndex( numRegions );
  std::vector< ElementType > cellBlockElementType( numRegions );

  std::size_t count = 0;
  auto fillCellBlockNames = [&]( std::map< int, std::vector< vtkIdType > > & regionsElem,
                                 ElementType type ) -> void
  {
    for( auto const & region: regionsElem )
    {
      cellBlockRegionIndex[count] = region.first;
      cellBlockElementType[count] = type;
      count++;
    }
  };

  fillCellBlockNames( regionsHex, ElementType::Hexahedron );
  fillCellBlockNames( regionsTetra, ElementType::Tetrahedron );
  fillCellBlockNames( regionsWedges, ElementType::Prism );
  fillCellBlockNames( regionsPyramids, ElementType::Pyramid );

  // Communicate all the region names
  // TODO if in a mesh, there is for instance a region with tetra, and a region without tetra,
  //      they will both contain a cell block tetrahedron
  std::vector< int > const allCellBlockRegionIndex = gatherData( cellBlockRegionIndex );
  std::vector< ElementType > const allCellBlockElementType = gatherData( cellBlockElementType );
  for( std::size_t i = 0; i < allCellBlockRegionIndex.size(); i++ )
  {
    for( std::size_t j = 0; j < allCellBlockElementType.size(); j++ )
    {
      switch( allCellBlockElementType[j] )
      {
        case ElementType::Tetrahedron:
          regionsTetra[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Hexahedron:
          regionsHex[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Prism:
          regionsWedges[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Pyramid:
          regionsPyramids[allCellBlockRegionIndex[i]];
          break;
        default:
          GEOSX_ERROR( allCellBlockElementType[j] << " is not a recognized cell type" );
      }
    }
  }

  // We want to know the surface attributes from all MPI ranks.
  // Then we'll be able to allocate an entry in the `surfacesIdsToCellsIds` table.
  // The `cellsIds` can possibly be empty (if the rank does not hold any cells for the surface).
  std::vector< int > surfaces;
  surfaces.reserve( surfacesIdsToCellsIds.size() );
  for( auto const & s2c: surfacesIdsToCellsIds )
  {
    surfaces.push_back( s2c.first );
  }

  for( auto const & s: gatherData( surfaces ) )
  {
    surfacesIdsToCellsIds[s];
  }

  return surfacesIdsToCellsIds;
}

/**
 * @brief Find the properties that can be imported.
 * @details All the float and double vtkArrays can be imported.
 * @return The list of the arrays.
 */
std::vector< vtkDataArray * > findImportableArrays( vtkSmartPointer< vtkUnstructuredGrid > mesh )
{
  std::vector< vtkDataArray * > importableArrays;
  vtkCellData * cellData = mesh->GetCellData();
  for( int i = 0; i < cellData->GetNumberOfArrays(); i++ )
  {
    vtkAbstractArray * curArray = cellData->GetAbstractArray( i );
    // We support only float and double type. Good assumption?
    int dataType = curArray->GetDataType();
    if( !curArray->IsNumeric() )
    {
      GEOSX_LOG_RANK_0( curArray->GetName() << " is not a numeric array" );
      continue;
    }
    if( dataType != VTK_FLOAT && dataType != VTK_DOUBLE )
    {
      // TODO Do not show a message for GlobalIds
      // TODO check on "attribute" + make a specific constant for this.
      GEOSX_LOG_RANK_0( "Underlying data Type " << curArray->GetDataTypeAsString() << " for array "
                                                << curArray->GetName() << " is not supported by GEOSX and will be ignored "
                                                << "(Only double and float are supported)" );
      continue;
    }
    importableArrays.push_back( vtkDataArray::SafeDownCast( curArray ) );
  }
  return importableArrays;
}


/**
 * @brief Get the GEOSX element type
 * @param[in] cellType The vtk cell type
 * @return The GEOSX element type
 */
ElementType convertVtkToGeosxElementType( VTKCellType cellType )
{
  switch( cellType )
  {
    case VTK_HEXAHEDRON:
      return ElementType::Hexahedron;
    case VTK_TETRA:
      return ElementType::Tetrahedron;
    case VTK_WEDGE:
      return ElementType::Prism;
    case VTK_PYRAMID:
      return ElementType::Pyramid;
    default:
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator" );
  }
  // This `return` is never reached but I faced a compilation error...
  return ElementType::Polyhedron;
}


/**
 * @brief Write the hexahedron vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in] cellIds the hex indexes within this region
 * @param[out] cellToVertex list of nodes organized per cells
 * @param[out] localToGlobal the local index to global index map
 */
void writeHexahedronVertices( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                              std::vector< vtkIdType > const & cellIds,
                              arrayView2d< localIndex, cells::NODE_MAP_USD > cellToVertex,
                              arrayView1d< globalIndex > localToGlobal )
{
  localIndex cellCount = 0;
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( mesh->GetCellData(), "GlobalCellIds" );
  for( vtkIdType c: cellIds )
  {
    vtkCell * currentCell = mesh->GetCell( c );
    cellToVertex[cellCount][0] = currentCell->GetPointId( 0 );
    cellToVertex[cellCount][1] = currentCell->GetPointId( 1 );
    cellToVertex[cellCount][2] = currentCell->GetPointId( 3 );
    cellToVertex[cellCount][3] = currentCell->GetPointId( 2 );
    cellToVertex[cellCount][4] = currentCell->GetPointId( 4 );
    cellToVertex[cellCount][5] = currentCell->GetPointId( 5 );
    cellToVertex[cellCount][6] = currentCell->GetPointId( 7 );
    cellToVertex[cellCount][7] = currentCell->GetPointId( 6 );

    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

/**
 * @brief Write the wedge vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in,out] cellToVertex list of nodes organized per cells
 * @param[in] cellIds the wedges indexes within this region
 * @param[in,out] localToGLobal the local index to global index map
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeWedgeVertices( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                         std::vector< vtkIdType > const & cellIds,
                         arrayView2d< localIndex, cells::NODE_MAP_USD > cellToVertex,
                         arrayView1d< globalIndex > localToGlobal )
{
  localIndex cellCount = 0;
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( mesh->GetCellData(), "GlobalCellIds" );
  for( vtkIdType c: cellIds )
  {
    vtkCell * currentCell = mesh->GetCell( c );
    cellToVertex[cellCount][0] = currentCell->GetPointId( 0 );
    cellToVertex[cellCount][1] = currentCell->GetPointId( 3 );
    cellToVertex[cellCount][2] = currentCell->GetPointId( 1 );
    cellToVertex[cellCount][3] = currentCell->GetPointId( 4 );
    cellToVertex[cellCount][4] = currentCell->GetPointId( 2 );
    cellToVertex[cellCount][5] = currentCell->GetPointId( 5 );
    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

/**
 * @brief Fill @p cellBlock with the appropriate nodes and local/global mappings.
 * @param[in] cellType the vtk cell type for cells of the CellBlock being written
 * @param[in] cellIds the cell indexes of cell type \p cellType within this region
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in,out] cellBlock The cell block to be written
 */
void fillCellBlock( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                    int cellType,
                    std::vector< vtkIdType > const & cellIds,
                    CellBlock & cellBlock )
{
  localIndex const numNodesPerElement = cellBlock.numNodesPerElement();
  arrayView2d< localIndex, cells::NODE_MAP_USD > const cellToVertex = cellBlock.getElemToNode();
  arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

  // Writing connectivity and Local to Global
  if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
  {
    writeHexahedronVertices( mesh, cellIds, cellToVertex, localToGlobal );
    return;
  }
  if( cellType == VTK_WEDGE ) // Special case for wedge because of the ordering
  {
    writeWedgeVertices( mesh, cellIds, cellToVertex, localToGlobal );
    return;
  }
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( mesh->GetCellData(), "GlobalCellIds" );
  localIndex cellCount = 0;
  for( vtkIdType c: cellIds )
  {
    vtkCell * currentCell = mesh->GetCell( c );
    for( localIndex v = 0; v < numNodesPerElement; v++ )
    {
      cellToVertex[cellCount][v] = currentCell->GetPointId( v );
    }
    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

/**
 * @brief Builds the cell block name.
 * @param[in] elementName The element name.
 * @param[in] regionId The region Id.
 * @return The name.
 * @warning This name will be visible in the input file... Consider refactoring with great care.
 */
string buildCellBlockName( VTKCellType vtkCellType,
                           int regionId )
{
  GEOSX_ASSERT( regionId >= -1 );
  string const cellTypeName = getCellName( vtkCellType );
  return regionId != -1 ? std::to_string( regionId ) + "_" + cellTypeName : cellTypeName;
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

/**
 * @brief Imports 2d and 3d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importMaterialField( std::vector< vtkIdType > const & cellIds,
                          vtkDataArray * vtkArray,
                          WrapperBase & wrapper )
{
  // Scalar material fields are stored as 2D arrays, vector/tensor are 3D
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 2, 3 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex cellCount = 0;
    for( vtkIdType c: cellIds )
    {
      int componentIdx = 0;
      for( localIndex q = 0; q < view.size( 1 ); ++q )
      {
        // The same value is copied for all quadrature points.
        LvArray::forValuesInSlice( view[cellCount][q], [&]( auto & val )
        {
          val = vtkArray->GetComponent( c, componentIdx );
        } );
      }
      cellCount++;
      componentIdx++;
    }
  } );
}

/**
 * @brief Imports 1d and 2d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importRegularField( std::vector< vtkIdType > const & cellIds,
                         vtkDataArray * vtkArray,
                         WrapperBase & wrapper )
{
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 1, 2 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex cellCount = 0;
    for( vtkIdType c: cellIds )
    {
      LvArray::forValuesInSlice( view[cellCount], [&]( auto & val )
      {
        val = vtkArray->GetComponent( c, 0 );
      } );
      cellCount++;
    }
  } );
}

void importFieldOnCellElementSubRegion( int regionId,
                                        VTKCellType vtkCellType,
                                        std::vector< vtkIdType > const & cellIds,
                                        ElementRegionManager & elemManager,
                                        std::vector< vtkDataArray * > const & importableArrays )
{
  string const cellBlockName = buildCellBlockName( vtkCellType, regionId );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
  {
    if( subRegion.getName() == cellBlockName )
    {
      // Writing properties
      std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( subRegion );

      for( vtkDataArray * vtkArray: importableArrays )
      {
        if( materialWrapperNames.count( vtkArray->GetName() ) == 0 )
        {
          continue;
        }
        WrapperBase & wrapper = subRegion.getWrapperBase( vtkArray->GetName() );
        if( materialWrapperNames.count( vtkArray->GetName() ) > 0 && wrapper.numArrayDims() > 1 )
        {
          importMaterialField( cellIds, vtkArray, wrapper );
        }
        else
        {
          importRegularField( cellIds, vtkArray, wrapper );
        }
      }
    }
  } );
}

void VTKMeshGenerator::importFields( DomainPartition & domain ) const
{
  // TODO Having CellElementSubRegion and ConstitutiveBase... here in a pure geometric module is problematic.
  ElementRegionManager & elemManager = domain.getMeshBody( this->getName() ).getMeshLevel( 0 ).getElemManager();

  auto importFieldsForElementType = [&]( VTKCellType vtkCellType,
                                         std::map< int, std::vector< vtkIdType > > const & cellBlocks )
  {
    for( auto & cellBlock: cellBlocks )
    {
      importFieldOnCellElementSubRegion( cellBlock.first, vtkCellType, cellBlock.second, elemManager, m_importableArrays );
    }
  };

  importFieldsForElementType( VTK_HEXAHEDRON, m_regionsHex );
  importFieldsForElementType( VTK_TETRA, m_regionsTetra );
  importFieldsForElementType( VTK_WEDGE, m_regionsWedges );
  importFieldsForElementType( VTK_PYRAMID, m_regionsPyramids );

  std::map< string, string_array > fieldNames;
  for( vtkDataArray * vtkArray: m_importableArrays )
  {
    fieldNames["elems"].emplace_back( vtkArray->GetName() );
  }

  CommunicationTools::getInstance().synchronizeFields( fieldNames, domain.getMeshBody( this->getName() ).getMeshLevel( 0 ), domain.getNeighbors(), false );
}

/**
 * @brief Build all the cell blocks.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in] regionsHex map from region index to the hexahedron indexes in this region
 * @param[in] regionsTetra map from region index to the tetra indexes in this region
 * @param[in] regionsWedges map from region index to the wedges indexes in this region
 * @param[in] regionsPyramids map from region index to the pyramids indexes in this region
 * @param[out] cellBlockManager The instance that stores the cell blocks.
 */
void buildCellBlocks( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                      std::map< int, std::vector< vtkIdType > > const & regionsHex,
                      std::map< int, std::vector< vtkIdType > > const & regionsTetra,
                      std::map< int, std::vector< vtkIdType > > const & regionsWedges,
                      std::map< int, std::vector< vtkIdType > > const & regionsPyramids,
                      CellBlockManagerBase & cellBlockManager )
{
  // Creates a new cell block for each region and for each type of cell.
  auto fct = [&]( VTKCellType vtkType, std::map< int, std::vector< vtkIdType > > const & regionIdToCellIds ) -> void
  {
    for( auto const & r2c: regionIdToCellIds )
    {
      int const regionId = r2c.first;
      std::vector< vtkIdType > const & cellIds = r2c.second;

      string const cellBlockName = buildCellBlockName( vtkType, regionId );

      // Create and resize the cell block.
      CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );
      cellBlock.setElementType( convertVtkToGeosxElementType( vtkType ) );
      cellBlock.resize( cellIds.size() );

      fillCellBlock( mesh, vtkType, cellIds, cellBlock );
    }
  };

  fct( VTK_HEXAHEDRON, regionsHex );
  fct( VTK_TETRA, regionsTetra );
  fct( VTK_WEDGE, regionsWedges );
  fct( VTK_PYRAMID, regionsPyramids );
}

/**
 * @brief Build the "surface" node sets from the surface information.
 * @param[in] mesh The vtkUnstructuredGrid that is loaded
 * @param[in] surfacesIdsToCellsIds Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[out] cellBlockManager The instance that stores the node sets.
 * @note @p surfacesIdsToCellsIds will contain all the surface ids across all the MPI ranks, but only its cell ids.
 * If the current MPI rank has no cell id for a given surface, then an empty set will be created.
 */
void buildSurfaces( vtkSmartPointer< vtkUnstructuredGrid > mesh,
                    std::map< int, std::vector< vtkIdType > > const & surfacesIdsToCellsIds,
                    CellBlockManagerBase & cellBlockManager )
{
  std::map< string, SortedArray< localIndex > > & nodeSets = cellBlockManager.getNodeSets();

  for( auto const & s2c: surfacesIdsToCellsIds )
  {
    int const & surfaceId = s2c.first;
    std::vector< vtkIdType > const & cellIds = s2c.second;

    // Get or create all surfaces (even those which are empty in this rank)
    SortedArray< localIndex > & curNodeSet = nodeSets[ std::to_string( surfaceId ) ];

    for( vtkIdType const & c: cellIds )
    {
      vtkCell * currentCell = mesh->GetCell( c );
      for( localIndex v = 0; v < currentCell->GetNumberOfPoints(); v++ )
      {
        curNodeSet.insert( currentCell->GetPointId( v ) );
      }
    }
  }
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  // TODO refactor void MeshGeneratorBase::generateMesh( DomainPartition & domain )
  GEOSX_MARK_FUNCTION;

  int mpiSize = MpiWrapper::commSize();
  GEOSX_ERROR_IF( ( mpiSize & ( mpiSize - 1 ) ) != 0, "MPI size is not a power of 2. Can't be used with the VTKMeshGenerator" );
  vtkSmartPointer< vtkMultiProcessController > controller = getVTKController();
  vtkMultiProcessController::SetGlobalController( controller );

  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = loadVTKMesh( m_filePath );

  std::vector< vtkBoundingBox > cuts;
  m_vtkMesh = redistributeMesh( *loadedMesh, cuts );

  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
  meshBody.getMeshLevels().registerGroup< MeshLevel >( "Level0" );

  //CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

  // Ongoing test 
  // TODO: Modify CellBlockManagerABC so that it has the basic functions used by the MeshGenerator
  // The class CellBlockManager should not be a type for functions here
  HexCellBlockManager & cellBlockManager = meshBody.registerGroup< HexCellBlockManager >( keys::cellManager );
  
  
  double const globalLength = writeMeshNodes( cellBlockManager, m_vtkMesh );
  meshBody.setGlobalLengthScale( globalLength );

  // TODO Check that the neighbor information set is bulletproof
  domain.getMetisNeighborList() = computeMPINeighborRanks( cuts );

  std::map< int, std::vector< vtkIdType > > const 
  surfacesIdsToCellsIds = buildRegionToCellsAndFaces( m_vtkMesh, m_regionsHex, m_regionsTetra, m_regionsWedges, m_regionsPyramids );

  m_importableArrays = findImportableArrays( m_vtkMesh );

  buildCellBlocks( m_vtkMesh, m_regionsHex, m_regionsTetra, m_regionsWedges, m_regionsPyramids, cellBlockManager );

  buildSurfaces( m_vtkMesh, surfacesIdsToCellsIds, cellBlockManager );

  // TODO Check the memory usage that seems prohibitive - Do we need to build all connections?
  cellBlockManager.buildMaps();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geosx
