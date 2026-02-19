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
 * @file BlueprintOutput.cpp
 */

#include "BlueprintOutput.hpp"

#include "CommonMeshOutputUtilities.hpp"
#include "common/TimingMacros.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"

#include <conduit.hpp>
#include <conduit_blueprint.hpp>

namespace geos
{

BlueprintOutput::BlueprintOutput( string const & name,
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

bool BlueprintOutput::execute( real64 const time_n,
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

    conduit::Node const & mesh = meshRoot.fetch_existing( "mesh" );

    conduit::Node info;
    GEOS_ASSERT_MSG( conduit::blueprint::verify( "mesh", meshRoot, info ), info.to_json() );

    conduit::Node fileRoot;
    conduit::Node & index = fileRoot[ "blueprint_index/mesh" ];
    conduit::index_t const numberOfDomains = m_writeParallelFiles != 0 ? MpiWrapper::commSize() : 1;
    conduit::blueprint::mesh::generate_index( mesh, "mesh", numberOfDomains, index );

    info.reset();
    GEOS_ASSERT_MSG( conduit::blueprint::mesh::index::verify( index, info ), info.to_json() );

    outputUtilities::writeCommonMeshHDF5( meshRoot,
                                          getOutputDirectory(),
                                          cycleNumber,
                                          outputUtilities::CommonHDFWriteOptions{
                                            "blueprintFiles",
                                            m_writeParallelFiles
                                          },
                                          &fileRoot );
  }

  return false;
}

REGISTER_CATALOG_ENTRY( OutputBase, BlueprintOutput, string const &, dataRepository::Group * const )

} /* namespace geos */
