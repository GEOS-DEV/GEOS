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
 * @file XDMFOutput.cpp
 */

#include "XDMFOutput.hpp"

#include "CommonMeshOutputUtilities.hpp"
#include "common/Path.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"

#include <conduit.hpp>

#include <fstream>
#include <sstream>
#include <unordered_map>

namespace geos
{
namespace internal
{

struct TopologyInfo
{
  string xdmfType;
  localIndex nodesPerElement;
};

TopologyInfo toXdmfTopology( string const & blueprintShape )
{
  if( blueprintShape == "hex" )
  {
    return { "Hexahedron", 8 };
  }
  if( blueprintShape == "tet" )
  {
    return { "Tetrahedron", 4 };
  }

  GEOS_ERROR( "Unsupported Blueprint shape for XDMF: " << blueprintShape );
  return {};
}

} // namespace internal

XDMFOutput::XDMFOutput( string const & name,
                        dataRepository::Group * const parent ):
  OutputBase( name, parent )
{
  registerWrapper( "plotLevel", &m_plotLevel ).
    setApplyDefaultValue( dataRepository::PlotLevel::LEVEL_1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setRTTypeName( rtTypes::CustomTypes::plotLevel ).
    setDescription( "Determines which fields to write." );

  registerWrapper( "outputFullQuadratureData", &m_outputFullQuadratureData ).
    setApplyDefaultValue( false ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If true writes out data associated with every quadrature point." );

  registerWrapper( "writeParallelFiles", &m_writeParallelFiles ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If non-zero, write one heavy-data file per rank." );

  registerWrapper( "writeGhostObjects", &m_writeGhostObjects ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If zero, skip ghost-owned cell objects." );

  registerWrapper( "writeFieldData", &m_writeFieldData ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If zero, skip all field data." );

  registerWrapper( "writeGhostFieldData", &m_writeGhostFieldData ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "If zero, skip fields that contain ghost metadata." );
}

bool XDMFOutput::execute( real64 const time_n,
                          real64 const GEOS_UNUSED_PARAM( dt ),
                          integer const cycleNumber,
                          integer const GEOS_UNUSED_PARAM( eventCounter ),
                          real64 const GEOS_UNUSED_PARAM( eventProgress ),
                          DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  {
    Timer timer( m_outputTimer );

    MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();

    conduit::Node meshRoot;
    outputUtilities::buildCommonMeshTree( meshLevel,
                                          time_n,
                                          cycleNumber,
                                          outputUtilities::CommonMeshOptions{
                                            m_plotLevel,
                                            m_outputFullQuadratureData,
                                            m_writeGhostObjects,
                                            m_writeFieldData,
                                            m_writeGhostFieldData
                                          },
                                          *this,
                                          meshRoot );

    outputUtilities::CommonHDFWriteResult const writeResult =
      outputUtilities::writeCommonMeshHDF5( meshRoot,
                                            getOutputDirectory(),
                                            cycleNumber,
                                            outputUtilities::CommonHDFWriteOptions{
                                              "xdmfFiles",
                                              m_writeParallelFiles
                                            } );

    if( writeResult.wroteFileOnThisRank == 0 )
    {
      return false;
    }

    string const xdmfPath = GEOS_FMT( "{}/rank_{:07}.xmf", writeResult.cyclePath, writeResult.fileRank );
    std::ofstream xdmfFile( xdmfPath );
    GEOS_ERROR_IF( !xdmfFile.is_open(), GEOS_FMT( "Failed to open XDMF file for writing: {}", xdmfPath ) );
    xdmfFile << buildXdmfDocument( meshRoot.fetch_existing( "mesh" ), writeResult.hdfFileName, time_n );
  }

  return false;
}

string XDMFOutput::buildXdmfDocument( conduit::Node const & mesh,
                                      string const & hdfFileName,
                                      real64 const time )
{
  conduit::Node const & coordset = mesh.fetch_existing( "coordsets/nodes" );
  conduit::Node const & topologies = mesh.fetch_existing( "topologies" );
  conduit::Node const * fields = mesh.has_path( "fields" ) ? &mesh.fetch_existing( "fields" ) : nullptr;

  string const coordsetName = "nodes";
  string const xPath = hdfFileName + ":/mesh/coordsets/nodes/values/x";
  string const yPath = hdfFileName + ":/mesh/coordsets/nodes/values/y";
  string const zPath = hdfFileName + ":/mesh/coordsets/nodes/values/z";
  localIndex const numberOfNodes = coordset.fetch_existing( "values/x" ).dtype().number_of_elements();

  std::unordered_map< string, stdVector< conduit::Node const * > > elementFieldsByTopology;
  stdVector< conduit::Node const * > nodalFields;
  if( fields != nullptr )
  {
    for( conduit::index_t i = 0; i < fields->number_of_children(); ++i )
    {
      conduit::Node const & field = fields->child( i );
      string const topology = field.fetch_existing( "topology" ).as_string();
      if( topology == coordsetName )
      {
        nodalFields.emplace_back( &field );
      }
      else
      {
        elementFieldsByTopology[topology].emplace_back( &field );
      }
    }
  }

  std::ostringstream os;
  os << "<?xml version=\"1.0\" ?>\n";
  os << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  os << "<Xdmf Version=\"3.0\">\n";
  os << "  <Domain>\n";

  for( conduit::index_t i = 0; i < topologies.number_of_children(); ++i )
  {
    conduit::Node const & topology = topologies.child( i );
    string const topologyName = topology.name();
    string const shape = topology.fetch_existing( "elements/shape" ).as_string();
    if( shape == "point" )
    {
      continue;
    }

    internal::TopologyInfo const topoInfo = internal::toXdmfTopology( shape );
    conduit::Node const & connectivity = topology.fetch_existing( "elements/connectivity" );
    localIndex const totalConnectivitySize = connectivity.dtype().number_of_elements();
    localIndex const numberOfElements = totalConnectivitySize / topoInfo.nodesPerElement;
    string const connectivityPath = hdfFileName + ":/mesh/topologies/" + topologyName + "/elements/connectivity";

    os << "    <Grid Name=\"" << topologyName << "\" GridType=\"Uniform\">\n";
    os << "      <Time Value=\"" << time << "\" />\n";
    os << "      <Topology TopologyType=\"" << topoInfo.xdmfType << "\" NumberOfElements=\"" << numberOfElements << "\">\n";
    os << "        <DataItem Format=\"HDF\" Dimensions=\"" << numberOfElements << " " << topoInfo.nodesPerElement << "\">"
       << connectivityPath << "</DataItem>\n";
    os << "      </Topology>\n";
    os << "      <Geometry GeometryType=\"X_Y_Z\">\n";
    os << "        <DataItem Format=\"HDF\" Dimensions=\"" << numberOfNodes << "\">" << xPath << "</DataItem>\n";
    os << "        <DataItem Format=\"HDF\" Dimensions=\"" << numberOfNodes << "\">" << yPath << "</DataItem>\n";
    os << "        <DataItem Format=\"HDF\" Dimensions=\"" << numberOfNodes << "\">" << zPath << "</DataItem>\n";
    os << "      </Geometry>\n";

    for( conduit::Node const * field : nodalFields )
    {
      string const fieldName = field->name();
      localIndex const fieldSize = field->fetch_existing( "values" ).dtype().number_of_elements();
      string const fieldPath = hdfFileName + ":/mesh/fields/" + fieldName + "/values";
      os << "      <Attribute Name=\"" << fieldName << "\" Center=\"Node\" AttributeType=\"Scalar\">\n";
      os << "        <DataItem Format=\"HDF\" Dimensions=\"" << fieldSize << "\">" << fieldPath << "</DataItem>\n";
      os << "      </Attribute>\n";
    }

    auto const topologyFieldsIter = elementFieldsByTopology.find( topologyName );
    if( topologyFieldsIter != elementFieldsByTopology.end() )
    {
      for( conduit::Node const * field : topologyFieldsIter->second )
      {
        string const fieldName = field->name();
        localIndex const fieldSize = field->fetch_existing( "values" ).dtype().number_of_elements();
        string const fieldPath = hdfFileName + ":/mesh/fields/" + fieldName + "/values";
        os << "      <Attribute Name=\"" << fieldName << "\" Center=\"Cell\" AttributeType=\"Scalar\">\n";
        os << "        <DataItem Format=\"HDF\" Dimensions=\"" << fieldSize << "\">" << fieldPath << "</DataItem>\n";
        os << "      </Attribute>\n";
      }
    }

    os << "    </Grid>\n";
  }

  os << "  </Domain>\n";
  os << "</Xdmf>\n";
  return os.str();
}

REGISTER_CATALOG_ENTRY( OutputBase, XDMFOutput, string const &, dataRepository::Group * const )

} /* namespace geos */
