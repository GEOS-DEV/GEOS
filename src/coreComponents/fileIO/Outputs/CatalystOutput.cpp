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
 * @file CatalystOutput.cpp
 */

//GEOS
#include "CatalystOutput.hpp"

#include "fileIO/Catalyst/CatalystActions.hpp"
#include "fileIO/Catalyst/GenericConduitCapsule.tpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

//TPL
#include <conduit.hpp>
#include <conduit_cpp_to_c.hpp>

//std
#include <cstdlib>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

// Change all the node associated fields to have vertex association and delete the "nodes" element entry
void SanitizeNode(conduit::Node& blueprintNode)
{
  if (blueprintNode["topologies"].number_of_children() > 2)
  {
    GEOS_ERROR("CatalystOutput: Not yet support for multiple regions");
    return;
  }
  blueprintNode.remove("topologies/nodes");
  std::string topoName = blueprintNode["topologies"].child(0).name();
  auto& fields = blueprintNode["fields"];
  for ( conduit::index_t iField = 0; iField < fields.number_of_children(); ++iField )
  {
    auto& field = fields.child(iField);
    if (field["topology"].as_char8_str() == std::string("nodes"))
    {
      field["topology"] = topoName;
      field["association"] = "vertex";
    }
  }
}

}

namespace geos
{

struct CatalystOutput::CatalystInternals
{
  bool initialized = false;
  std::string scripts = "";
  std::string implementation = "";
  std::string implementationPath = "";
  std::string adiosConfig = "";
  std::string channelName = "fullfield";

  void initializeCatalyst()
  {

    std::string scriptsCopy = this->scripts;
    std::vector<std::string> scriptList;
    std::string delimiter = ":";
    std::size_t pos = scriptsCopy.find(delimiter);
    scriptList.emplace_back(scriptsCopy.substr(0, pos));
    scriptsCopy.erase(0, pos + delimiter.length());
    while ((pos = scriptsCopy.find(delimiter)) != std::string::npos) {
      scriptList.emplace_back(scriptsCopy.substr(0, pos));
      scriptsCopy.erase(0, pos + delimiter.length());
    }

    if ( scriptList.empty() || scriptList[0] == "" )
    {
      GEOS_ERROR("CatalystOutput: Constructor: no catalyst scripts found.");
    }

    conduit::Node initializer;
    for (std::size_t iScr = 0; iScr < scriptList.size(); ++iScr)
    {
      initializer["catalyst/scripts/script" + std::to_string(iScr)] = scriptList[iScr];
    }

    if ( !this->implementation.empty() )
    {
      initializer["catalyst_load/implementation"] = this->implementation;
      if ( !this->implementationPath.empty() )
      {
        initializer["catalyst_load/search_paths/" + this->implementation] = this->implementationPath;
      }
    }

    if ( this->implementation == "adios" || 
        (this->implementation.empty() && 
         std::string(std::getenv("CATALYST_IMPLEMENTATION_NAME")) == "adios"))
    {
      if ( this->adiosConfig == "" )
      {
        GEOS_ERROR("CatalystOutput: Constructor: Cannot use the catalyst/adios2 implementation without providing an adios2 configuration file.");
      }
      initializer["adios/config"] = this->adiosConfig;
    }

    auto capsule = GenericConduitCapsule<conduit::Node>(&initializer);
    if( !CatalystInitialize(&capsule) )
    {
      GEOS_ERROR("CatalystOutput: Constructor: catalyst failed to initialize.");
    }

    this->initialized = true;
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
CatalystOutput::CatalystOutput( std::string const & name,
                                Group * const parent )
               : BlueprintOutput( name, 
                                  parent )
               , internal( std::unique_ptr<CatalystInternals>( new CatalystInternals() ) )
{

  this->registerWrapper( "scripts", &this->internal->scripts ).
    setApplyDefaultValue( "" ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Column separated paths to the catalyst scripts." );

  this->registerWrapper( "implementation", &this->internal->implementation ).
    setApplyDefaultValue( "" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name of the catalyst implementation to use." );

  this->registerWrapper( "implementationPath", &this->internal->implementationPath ).
    setApplyDefaultValue( "" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Path to the catalyst the implementation to use." );

  this->registerWrapper( "adiosConfig", &this->internal->adiosConfig ).
    setApplyDefaultValue( "" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Path to the adios configuration file when using the catalyst-adios implementation." );

  this->registerWrapper( "fullFieldChannelName", &this->internal->channelName ).
    setApplyDefaultValue( "" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Name to give to the channel passing the full field data." );

}

///////////////////////////////////////////////////////////////////////////////////////////////////
CatalystOutput::~CatalystOutput() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
bool CatalystOutput::execute( real64 const time_n,
                              real64 const /*dt*/,
                              integer const cycleNumber,
                              integer const /*eventCounter*/,
                              real64 const /*eventProgress*/,
                              DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  if ( !this->internal->initialized )
  {
    this->internal->initializeCatalyst();
  }

  conduit::Node executeRoot;
  auto& catalystState = executeRoot["catalyst/state"];
  catalystState["timestep"].set(cycleNumber);
  catalystState["time"].set(time_n);

  auto& channel = executeRoot["catalyst/channels/" + this->internal->channelName];
  channel["type"] = "mesh";

  auto& meshGEOSRoot = channel["data"];
  this->mapMesh(time_n, cycleNumber, domain, meshGEOSRoot);
  ::SanitizeNode(meshGEOSRoot);

  // Add additional info if the implementation is adios based
  if ( this->internal->implementation == "adios" || 
        (this->internal->implementation.empty() && 
         std::string(std::getenv("CATALYST_IMPLEMENTATION_NAME")) == "adios") )
  {
    int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
    int const commSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

    MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();

    localIndex_array numberOfNodesPerRank( commSize );

    MpiWrapper::allGather( 
        meshLevel.getNodeManager().size(), 
        numberOfNodesPerRank );

    globalIndex offsetPoints = std::accumulate(
        numberOfNodesPerRank.begin(),
        numberOfNodesPerRank.begin() + thisRank,
        0 );

    globalIndex totalNumberOfNodes = std::accumulate(
        numberOfNodesPerRank.begin() + thisRank, 
        numberOfNodesPerRank.end(), 
        offsetPoints );

    auto& mesh = channel[ "data" ];
    mesh[ "coordsets/coords/points_shape" ].set_uint64(totalNumberOfNodes);
    mesh[ "coordsets/coords/points_start" ].set_uint64(offsetPoints);
    mesh[ "coordsets/coords/points_count" ].set_uint64(numberOfNodesPerRank[ thisRank ]);

    // WARNING
    // what follows on the elements only works if all the MPI partitions have the same number
    // of corresponding subRegions. I do not know if this is the case. If it is not, we expose
    // ourselves to hanging MPI communication here.
    auto& topologies = mesh[ "topologies" ];
    meshLevel.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >(
        [&] ( localIndex, localIndex, ElementRegionBase const & region, CellElementSubRegion const & subRegion )
        {
        std::string const topologyName = region.getName() + "-" + subRegion.getName();
        auto& topology = topologies[ topologyName ];

        localIndex_array numberOfElementsPerRank( commSize );

        MpiWrapper::allGather(
            subRegion.size(), 
            numberOfElementsPerRank );

        globalIndex offsetElements = std::accumulate(
            numberOfElementsPerRank.begin(),
            numberOfElementsPerRank.begin() + thisRank,
            0 );

        globalIndex totalNumberOfElements = std::accumulate(
            numberOfElementsPerRank.begin() + thisRank,
            numberOfElementsPerRank.end(),
            offsetElements );

        topology["cells_shape"].set_uint64(totalNumberOfElements);
        topology["cells_start"].set_uint64(offsetElements);
        topology["cells_count"].set_uint64(numberOfElementsPerRank[ thisRank ]);
        }); 
  }

  auto capsule = GenericConduitCapsule<conduit::Node>(&executeRoot);
  if ( !CatalystExecute(&capsule) )
  {
    GEOS_ERROR("CatalystOutput: execute: catalyst failed to execute.");
    return true;
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void CatalystOutput::cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const eventCounter,
                        real64 const eventProgress,
                        DomainPartition & domain )
{ 
  GEOS_MARK_FUNCTION;

  if ( this->execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain ) )
  {
    GEOS_ERROR("CatalystOutput: cleanup: last execute failed.");
  }

  conduit::Node emptyNode;
  auto capsule = GenericConduitCapsule<conduit::Node>(&emptyNode);
  if ( !CatalystFinalize(&capsule) )
  {
    GEOS_ERROR("CatalystOutput: cleanup: finalize catalyst failed");
  }

  this->internal->initialized = false;
}

REGISTER_CATALOG_ENTRY( OutputBase, CatalystOutput, string const &, dataRepository::Group * const )

}
