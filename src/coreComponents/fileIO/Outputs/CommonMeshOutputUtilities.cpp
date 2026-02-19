/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CommonMeshOutputUtilities.cpp
 */

#include "CommonMeshOutputUtilities.hpp"

#include "common/MpiWrapper.hpp"
#include "common/Path.hpp"
#include "common/TimingMacros.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "dataRepository/wrapperHelpers.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ObjectManagerBase.hpp"

#include <conduit.hpp>
#include <conduit_relay.hpp>

namespace geos
{
namespace outputUtilities
{
namespace internal
{

string toBlueprintShape( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Tetrahedron: return "tet";
    case ElementType::Hexahedron: return "hex";
    default:
    {
      GEOS_ERROR( "No Blueprint type for element type: " << elementType );
      return {};
    }
  }
}

static stdVector< int > getBlueprintNodeOrdering( ElementType const elementType )
{
  switch( elementType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 };
    case ElementType::Polygon:       return { };
    case ElementType::Tetrahedron:   return { 1, 0, 2, 3 };
    case ElementType::Pyramid:       return { 0, 3, 2, 1, 4, 0, 0, 0 };
    case ElementType::Wedge:         return { 0, 4, 2, 1, 5, 3, 0, 0 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Prism5:        return { };
    case ElementType::Prism6:        return { };
    case ElementType::Prism7:        return { };
    case ElementType::Prism8:        return { };
    case ElementType::Prism9:        return { };
    case ElementType::Prism10:       return { };
    case ElementType::Prism11:       return { };
    case ElementType::Polyhedron:    return { };
  }
  return {};
}

bool isGhostField( string const & fieldName )
{
  return fieldName == ObjectManagerBase::viewKeyStruct::ghostRankString();
}

void collectElementIndicesToWrite( CellElementSubRegion const & subRegion,
                                   int const writeGhostObjects,
                                   stdVector< localIndex > & elementIndices,
                                   bool & contiguousPrefix )
{
  localIndex const numElems = subRegion.size();
  elementIndices.clear();
  elementIndices.reserve( numElems );

  if( writeGhostObjects != 0 )
  {
    for( localIndex ei = 0; ei < numElems; ++ei )
    {
      elementIndices.emplace_back( ei );
    }
    contiguousPrefix = true;
    return;
  }

  arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
  for( localIndex ei = 0; ei < numElems; ++ei )
  {
    if( ghostRank[ei] < 0 )
    {
      elementIndices.emplace_back( ei );
    }
  }

  contiguousPrefix = true;
  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( elementIndices.size() ); ++i )
  {
    if( elementIndices[i] != i )
    {
      contiguousPrefix = false;
      break;
    }
  }
}

void reorderElementToNodeMap( CellElementSubRegion const & subRegion,
                              stdVector< localIndex > const & elementIndices,
                              conduit::Node & connectivity )
{
  GEOS_MARK_FUNCTION;

  arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();
  localIndex const numElems = LvArray::integerConversion< localIndex >( elementIndices.size() );
  localIndex const numNodesPerElem = elemToNodeMap.size( 1 );

  stdVector< int > const vtkOrdering = getBlueprintNodeOrdering( subRegion.getElementType() );
  GEOS_ERROR_IF_NE( localIndex( vtkOrdering.size() ), numNodesPerElem );

  constexpr int conduitTypeID = dataRepository::conduitTypeInfo< localIndex >::id;
  conduit::DataType const dtype( conduitTypeID, numElems * numNodesPerElem );
  connectivity.set( dtype );

  localIndex * const reorderedConnectivity = connectivity.value();
  forAll< serialPolicy >( numElems, [reorderedConnectivity, numNodesPerElem, elemToNodeMap, &vtkOrdering, &elementIndices] ( localIndex const i )
  {
    localIndex const elemIndex = elementIndices[i];
    for( localIndex j = 0; j < numNodesPerElem; ++j )
    {
      reorderedConnectivity[ i * numNodesPerElem + j ] = elemToNodeMap( elemIndex, vtkOrdering[ j ] );
    }
  } );
}

void trimFieldValuesNode( conduit::Node & valuesNode, localIndex const count )
{
  conduit::DataType dtype = valuesNode.dtype();
  conduit::index_t const originalCount = dtype.number_of_elements();
  GEOS_ERROR_IF_LT( originalCount, conduit::index_t( count ) );
  if( originalCount == conduit::index_t( count ) )
  {
    return;
  }

  void * dataPtr = const_cast< void * >( valuesNode.data_ptr() );
  dtype.set_number_of_elements( count );
  valuesNode.set_external( dtype, dataPtr );
}

void trimNewFieldsToCount( conduit::Node & fields,
                           stdVector< string > const & newFieldNames,
                           localIndex const count,
                           bool const contiguousPrefix )
{
  if( count < 0 )
  {
    return;
  }

  if( !contiguousPrefix )
  {
    for( string const & fieldName : newFieldNames )
    {
      fields.remove( fieldName );
    }

    GEOS_LOG_RANK_0( "Warning: dropped field output for a non-contiguous non-ghost element layout "
                     "(ghost filtering requested)." );
    return;
  }

  for( string const & fieldName : newFieldNames )
  {
    conduit::Node & field = fields.fetch_existing( fieldName );
    conduit::Node & values = field.fetch_existing( "values" );
    trimFieldValuesNode( values, count );
  }
}

void writeOutWrappersAsFields( dataRepository::Group const & group,
                               conduit::Node & fields,
                               string const & topology,
                               CommonMeshOptions const & options,
                               localIndex const countForFieldTrim,
                               bool const contiguousPrefix,
                               string const & prefix = "" )
{
  GEOS_MARK_FUNCTION;

  if( options.writeFieldData == 0 )
  {
    return;
  }

  group.forWrappers( [&] ( dataRepository::WrapperBase const & wrapper )
  {
    if( wrapper.getPlotLevel() > options.plotLevel || !wrapper.sizedFromParent() )
    {
      return;
    }

    if( options.writeGhostFieldData == 0 && isGhostField( wrapper.getName() ) )
    {
      return;
    }

    string const name = prefix.empty() ? wrapper.getName() : prefix + "-" + wrapper.getName();

    conduit::index_t const numChildrenBefore = fields.number_of_children();
    wrapper.addBlueprintField( fields, name, topology );

    stdVector< string > newFieldNames;
    newFieldNames.reserve( fields.number_of_children() - numChildrenBefore );
    for( conduit::index_t i = numChildrenBefore; i < fields.number_of_children(); ++i )
    {
      newFieldNames.emplace_back( fields.child( i ).name() );
    }
    trimNewFieldsToCount( fields, newFieldNames, countForFieldTrim, contiguousPrefix );
  } );
}

void writeOutConstitutiveData( dataRepository::Group const & constitutiveModel,
                               conduit::Node & fields,
                               string const & topology,
                               CommonMeshOptions const & options,
                               localIndex const countForFieldTrim,
                               bool const contiguousPrefix,
                               dataRepository::Group & averagedElementData )
{
  GEOS_MARK_FUNCTION;

  if( options.writeFieldData == 0 )
  {
    return;
  }

  dataRepository::Group & averagedConstitutiveData = averagedElementData.registerGroup( constitutiveModel.getName() );

  constitutiveModel.forWrappers( [&] ( dataRepository::WrapperBase const & wrapper )
  {
    if( wrapper.getPlotLevel() > options.plotLevel || !wrapper.sizedFromParent() )
    {
      return;
    }

    if( options.writeGhostFieldData == 0 && isGhostField( wrapper.getName() ) )
    {
      return;
    }

    string const fieldName = constitutiveModel.getName() + "-quadrature-averaged-" + wrapper.getName();
    conduit::index_t const numChildrenBefore = fields.number_of_children();
    averagedConstitutiveData.registerWrapper( wrapper.averageOverSecondDim( fieldName, averagedConstitutiveData ) ).
      addBlueprintField( fields, fieldName, topology );

    stdVector< string > newFieldNames;
    newFieldNames.reserve( fields.number_of_children() - numChildrenBefore );
    for( conduit::index_t i = numChildrenBefore; i < fields.number_of_children(); ++i )
    {
      newFieldNames.emplace_back( fields.child( i ).name() );
    }
    trimNewFieldsToCount( fields, newFieldNames, countForFieldTrim, contiguousPrefix );
  } );
}

void addNodalData( NodeManager const & nodeManager,
                   conduit::Node & coordset,
                   conduit::Node & topologies,
                   conduit::Node & fields,
                   CommonMeshOptions const & options )
{
  GEOS_MARK_FUNCTION;

  coordset[ "type" ] = "explicit";
  dataRepository::wrapperHelpers::populateMCArray( nodeManager.referencePosition(),
                                                   coordset[ "values" ],
                                                   { "x", "y", "z" } );

  string const coordsetName = coordset.name();
  conduit::Node & nodeTopology = topologies[ coordsetName ];
  nodeTopology[ "coordset" ] = coordsetName;
  nodeTopology[ "type" ] = "unstructured";
  nodeTopology[ "elements/shape" ] = "point";
  conduit::Node & connectivity = nodeTopology[ "elements/connectivity" ];

  localIndex const numNodes = nodeManager.size();
  constexpr int conduitTypeID = dataRepository::conduitTypeInfo< localIndex >::id;
  conduit::DataType const dtype( conduitTypeID, numNodes );
  connectivity.set( dtype );

  localIndex * const nodeIDs = connectivity.value();
  forAll< serialPolicy >( numNodes, [nodeIDs] ( localIndex const i )
  {
    nodeIDs[i] = i;
  } );

  writeOutWrappersAsFields( nodeManager, fields, coordsetName, options, -1, true );
}

void addElementData( ElementRegionManager const & elemRegionManager,
                     conduit::Node & coordset,
                     conduit::Node & topologies,
                     conduit::Node & fields,
                     CommonMeshOptions const & options,
                     dataRepository::Group & averagedElementData )
{
  GEOS_MARK_FUNCTION;

  elemRegionManager.forElementSubRegionsComplete< CellElementSubRegion >(
    [&] ( localIndex, localIndex, ElementRegionBase const & region, CellElementSubRegion const & subRegion )
  {
    stdVector< localIndex > elementIndices;
    bool contiguousPrefix = true;
    collectElementIndicesToWrite( subRegion, options.writeGhostObjects, elementIndices, contiguousPrefix );

    if( elementIndices.empty() )
    {
      return;
    }

    string const topologyName = region.getName() + "-" + subRegion.getName();

    conduit::Node & topology = topologies[ topologyName ];
    topology[ "coordset" ] = coordset.name();
    topology[ "type" ] = "unstructured";
    topology[ "elements/shape" ] = toBlueprintShape( subRegion.getElementType() );
    reorderElementToNodeMap( subRegion, elementIndices, topology[ "elements/connectivity" ] );

    localIndex const numElemsToWrite = LvArray::integerConversion< localIndex >( elementIndices.size() );
    writeOutWrappersAsFields( subRegion, fields, topologyName, options, numElemsToWrite, contiguousPrefix );

    dataRepository::Group & averagedSubRegionData = averagedElementData.registerGroup( topologyName );
    subRegion.getConstitutiveModels().forSubGroups( [&]( dataRepository::Group const & constitutiveModel )
    {
      writeOutConstitutiveData( constitutiveModel,
                                fields,
                                topologyName,
                                options,
                                numElemsToWrite,
                                contiguousPrefix,
                                averagedSubRegionData );

      if( options.outputFullQuadratureData )
      {
        writeOutWrappersAsFields( constitutiveModel,
                                  fields,
                                  topologyName,
                                  options,
                                  numElemsToWrite,
                                  contiguousPrefix,
                                  constitutiveModel.getName() );
      }
    } );
  } );
}

} // namespace internal

void buildCommonMeshTree( MeshLevel const & meshLevel,
                          real64 const time_n,
                          integer const cycleNumber,
                          CommonMeshOptions const & options,
                          dataRepository::Group & parentGroup,
                          conduit::Node & meshRoot )
{
  GEOS_MARK_FUNCTION;

  conduit::Node & mesh = meshRoot[ "mesh" ];
  conduit::Node & coordset = mesh[ "coordsets/nodes" ];
  conduit::Node & topologies = mesh[ "topologies" ];

  mesh[ "state/time" ] = time_n;
  mesh[ "state/cycle" ] = cycleNumber;

  internal::addNodalData( meshLevel.getNodeManager(), coordset, topologies, mesh[ "fields" ], options );

  dataRepository::Group averagedElementData( "averagedElementData", &parentGroup );
  internal::addElementData( meshLevel.getElemManager(),
                            coordset,
                            topologies,
                            mesh[ "fields" ],
                            options,
                            averagedElementData );

  if( mesh[ "fields" ].number_of_children() == 0 )
  {
    mesh.remove( "fields" );
  }
}

CommonHDFWriteResult writeCommonMeshHDF5( conduit::Node & meshRoot,
                                          string const & outputDirectory,
                                          integer const cycleNumber,
                                          CommonHDFWriteOptions const & options,
                                          conduit::Node * rootNode )
{
  GEOS_MARK_FUNCTION;

  GEOS_ERROR_IF( options.outputSubdirectory.empty(), "Output subdirectory must not be empty." );

  CommonHDFWriteResult result;
  result.numberOfFiles = options.writeParallelFiles != 0 ? MpiWrapper::commSize() : 1;
  result.fileRank = options.writeParallelFiles != 0 ? MpiWrapper::commRank() : 0;
  result.wroteFileOnThisRank = options.writeParallelFiles != 0 || MpiWrapper::commRank() == 0;
  result.cyclePath = GEOS_FMT( "{}/{}/cycle_{:07}", outputDirectory, options.outputSubdirectory, cycleNumber );
  result.filePathForRank = GEOS_FMT( "{}/rank_{:07}.hdf5", result.cyclePath, result.fileRank );
  result.hdfFileName = splitPath( result.filePathForRank ).second;

  conduit::Node fallbackRoot;
  conduit::Node & root = rootNode != nullptr ? *rootNode : fallbackRoot;

  if( options.writeParallelFiles != 0 )
  {
    result.filePathForRank = dataRepository::writeRootFile( root, result.cyclePath );
    result.hdfFileName = splitPath( result.filePathForRank ).second;
    conduit::relay::io::save( meshRoot, result.filePathForRank, "hdf5" );
    return result;
  }

  if( MpiWrapper::commRank() == 0 )
  {
    makeDirsForPath( result.cyclePath );
    string const rootFileName = splitPath( result.cyclePath ).second;

    root[ "protocol/name" ] = "hdf5";
    root[ "protocol/version" ] = CONDUIT_VERSION;
    root[ "number_of_files" ] = 1;
    root[ "file_pattern" ] = rootFileName + "/rank_%07d.hdf5";
    root[ "number_of_trees" ] = 1;
    root[ "tree_pattern" ] = "/";

    conduit::relay::io::save( root, result.cyclePath + ".root", "hdf5" );
    conduit::relay::io::save( meshRoot, result.filePathForRank, "hdf5" );
  }

  MpiWrapper::barrier( MPI_COMM_GEOS );
  return result;
}

} /* namespace outputUtilities */
} /* namespace geos */
