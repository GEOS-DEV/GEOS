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
 * @file MeshPartitioner.cpp
 */



#include "MeshPartitioner.hpp"
#include "GraphPartitionEngine.hpp"

#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include <vtkBoundingBox.h>

#ifdef GEOS_USE_TRILINOS
#include "mesh/graphs/ZoltanGraphColoring.hpp"
#else
#include "mesh/graphs/RLFGraphColoringMPI.hpp"
#endif

#include "mesh/mpiCommunications/NoOpEngine.hpp"
#ifdef GEOS_USE_PARMETIS
#include "mesh/mpiCommunications/ParMetisEngine.hpp"
#endif
#ifdef GEOS_USE_PTSCOTCH
#include "mesh/mpiCommunications/PTScotchEngine.hpp"
#endif


namespace geos
{
MeshPartitioner::MeshPartitioner( string const & name, Group * const parent )
  : DomainPartitioner( name, parent ),
  m_engine( nullptr ),
  m_engineType( "parmetis" ),
  m_numRefinements( 0 )
{
  // Register engine type
  registerWrapper( MeshPartitioner::viewKeyStruct::engineTypeString(), &m_engineType ).
    setApplyDefaultValue( "parmetis" ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Graph partitioning engine (parmetis, ptscotch, or noop)" );

  // Register number of refinements
  registerWrapper( MeshPartitioner::viewKeyStruct::numRefinementsString(), &m_numRefinements ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Number of refinement iterations for graph partitioning (ParMETIS only, ignored by other engines)" );
}


void MeshPartitioner::setEngine( GraphPartitionEngine & engine )
{
  m_engine = &engine;
}

GraphPartitionEngine & MeshPartitioner::getEngine() const
{
  GEOS_ERROR_IF( m_engine == nullptr,
                 "Graph partition engine has not been set. "
                 "Did you forget to call postInputInitialization()?" );
  return *m_engine;
}

void MeshPartitioner::processCommandLineOverrides( unsigned int const xparCL,
                                                   unsigned int const yparCL,
                                                   unsigned int const zparCL )
{
#if 0
  // If user provided command-line overrides for partition counts...
  if( xparCL != 0 || yparCL != 0 || zparCL != 0 )
  {
    // For mesh partitioners, we can support this by changing target partition count
    setPartitionCounts( xparCL, yparCL, zparCL );

    GEOS_LOG_RANK_0( GEOS_FMT( "MeshPartitioner: Command-line override - partitioning into {} parts ({} x {} x {})",
                               m_numPartitions, xparCL, yparCL, zparCL ) );
  }
  else
  {
    // Default: partition into commSize parts (1 partition per rank)
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );
    setPartitionCounts( 1, 1, mpiSize );

    GEOS_LOG_RANK_0( GEOS_FMT( "MeshPartitioner: Partitioning into {} parts (1 per rank)", mpiSize ) );
  }
#else
  // If user provided command-line overrides for partition counts...
  if( xparCL != 0 || yparCL != 0 || zparCL != 0 )
  {
    // ...warn that they are ignored for non-geometric partitioners
    GEOS_WARNING( "Partition counts from the command line (-x, -y, -z) are ignored when using a mesh-based partitioner." );
  }

  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  setPartitionCounts( 1, 1, mpiSize );

#endif

}



vtk::AllMeshes MeshPartitioner::partitionMeshes( vtk::AllMeshes & mesh, MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  // ====================================================================
  // STEP 1: Compute partitioning, specific to subclass
  // ====================================================================
  array1d< int64_t > const partitioning = computePartitioning( mesh, comm );

  // ====================================================================
  // STEP 2: Apply partitioning to redistribute VTK meshes
  // ====================================================================
  vtk::AllMeshes redistributedMesh =
    vtk::applyPartitioning( mesh, partitioning.toViewConst(), comm );

  // ====================================================================
  // STEP 3: Extract neighbor information from redistributed mesh
  // ====================================================================
  extractNeighborsFromMesh( redistributedMesh, comm );

  // ====================================================================
  // STEP 4: Compute graph coloring for communication scheduling
  // ====================================================================
  color();

  return redistributedMesh;
}



void MeshPartitioner::extractNeighborsFromMesh( vtk::AllMeshes const & mesh, MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  // Exchange bounding boxes to determine neighbors
  stdVector< vtkBoundingBox > const boxes =
    vtk::exchangeBoundingBoxes( *mesh.getMainMesh(), comm );

  // Find neighbor ranks based on bounding box overlaps
  stdVector< int > const neighbors = vtk::findNeighborRanks( std::move( boxes ) );

  // Store neighbors (if needed)
  setNeighborsRank( std::move( neighbors ) );
}



void MeshPartitioner::color()
{
  GEOS_MARK_FUNCTION;

  // Build adjacency list from neighbor ranks
  std::vector< camp::idx_t > adjncy;
  adjncy.reserve( m_neighborsRank.size() );
  std::copy( m_neighborsRank.begin(), m_neighborsRank.end(), std::back_inserter( adjncy ) );

#ifdef GEOS_USE_TRILINOS
  geos::graph::ZoltanGraphColoring coloring;
#else
  geos::graph::RLFGraphColoringMPI coloring;
#endif

  m_color = coloring.colorGraph( adjncy );

  // Verify that the coloring is valid (no two neighbors have the same color)
  if( !coloring.isColoringValid( adjncy, m_color ) )
  {
    GEOS_ERROR( "Graph coloring failed: an adjacent partition has the same color." );
  }

  m_numColors = coloring.getNumberOfColors( m_color );

  GEOS_LOG_RANK_0( GEOS_FMT( "MeshPartitioner: Coloring complete. Using {} colors", m_numColors ) );
}


void MeshPartitioner::postInputInitialization()
{
  DomainPartitioner::postInputInitialization();

  // Check if an engine was already created as a child group
  m_engine = this->getGroupPointer< GraphPartitionEngine >( viewKeyStruct::engineTypeString() );

  if( m_engine == nullptr )
  {
    // Engine not created as child group - create from engineType attribute
    GEOS_LOG_RANK_0( "Creating '" << m_engineType << "' engine for " << getName() );

    string const engineTypeLower = stringutilities::toLower( m_engineType );

    if( engineTypeLower == "noop" )
    {
      m_engine = &this->registerGroup< NoOpEngine >( viewKeyStruct::engineTypeString() );
    }
#ifdef GEOS_USE_PARMETIS
    else if( engineTypeLower == "parmetis" )
    {
      m_engine = &this->registerGroup< ParMetisEngine >( viewKeyStruct::engineTypeString() );
    }
#endif
#ifdef GEOS_USE_PTSCOTCH
    else if( engineTypeLower == "ptscotch" )
    {
      m_engine = &this->registerGroup< PTScotchEngine >( viewKeyStruct::engineTypeString() );
    }
#endif
    else
    {
      GEOS_THROW( GEOS_FMT( "Unknown graph partitioner engine type: {}'", m_engineType ), InputError );
    }
  }

  GEOS_ERROR_IF( m_engine == nullptr,
                 "Failed to create graph partition engine for " << getCatalogName() );

  // Transfer numRefinements to the engine
  m_engine->setNumRefinements( m_numRefinements );
}

GraphPartitionEngine * MeshPartitioner::createDefaultEngine()
{
  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );

  if( mpiSize == 1 )
  {
    GEOS_LOG_RANK_0( "Single-rank run: using NoOp engine (no mesh redistribution)" );
    return &this->registerGroup< NoOpEngine >( viewKeyStruct::engineTypeString() );
  }

#ifdef GEOS_USE_PARMETIS
  GEOS_LOG_RANK_0( "Multi-rank run: using ParMetis engine" );
  return &this->registerGroup< ParMetisEngine >( viewKeyStruct::engineTypeString() );
#elif defined(GEOS_USE_PTSCOTCH)
  GEOS_LOG_RANK_0( "Multi-rank run: using PT-Scotch engine" );
  return &this->registerGroup< PTScotchEngine >( viewKeyStruct::engineTypeString() );
#else
  GEOS_ERROR( getCatalogName() << " requires ParMetis or PT-Scotch for multi-rank runs.\n"
                               << "Please rebuild with PARMETIS or PTSCOTCH enabled." );
  return nullptr;
#endif
}

} // namespace geos
