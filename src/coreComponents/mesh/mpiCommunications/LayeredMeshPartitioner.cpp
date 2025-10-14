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

 #include "LayeredMeshPartitioner.hpp"
#include "mesh/generators/VTKUtilities.hpp"

namespace geos
{
using namespace dataRepository;

LayeredMeshPartitioner::LayeredMeshPartitioner( string const & name, Group * const parent )
  : MeshPartitioner( name, parent ),
  m_indexArrayName( "" ),
  m_numPartZ( 1 )
{
  // Only register parameters specific to LayeredMeshPartitioner
  registerWrapper( viewKeyStruct::indexArrayNameString(), &m_indexArrayName ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "Name of VTK cell data array containing [area, layer] indices" );

  registerWrapper( viewKeyStruct::numPartZString(), &m_numPartZ ).
    setApplyDefaultValue( 1 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Number of partitions in layer Z-direction" );
}

LayeredMeshPartitioner::~LayeredMeshPartitioner() = default;

void LayeredMeshPartitioner::postInputInitialization()
{
  // Call base class initialization (creates and sets engine)
  MeshPartitioner::postInputInitialization();

  // Validate LayeredMeshPartitioner-specific parameters
  GEOS_THROW_IF( m_indexArrayName.empty(),
                 "LayeredMeshPartitioner requires 'indexArrayName' to be specified", InputError );

  GEOS_THROW_IF_LE_MSG( m_numPartZ, 0,
                        GEOS_FMT( "LayeredMeshPartitioner requires 'numPartZ' > 0, got {}", m_numPartZ ), InputError );

  int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  GEOS_THROW_IF_NE_MSG( mpiSize % m_numPartZ, 0,
                        GEOS_FMT( "Total MPI ranks ( {} ) must be evenly divisible by numPartZ ({})", mpiSize, m_numPartZ ), InputError );

  GEOS_LOG_RANK_0( GEOS_FMT( "LayeredMeshPartitioner: {} ranks will be partitioned into {} area partitions x {} layer partitions",
                             mpiSize, mpiSize / m_numPartZ, m_numPartZ ) );
}

void LayeredMeshPartitioner::processCommandLineOverrides( unsigned int const xparCL,
                                                          unsigned int const yparCL,
                                                          unsigned int const zparCL )
{
  // For layered partitioner:
  // - zparCL -> numPartZ (layer partitions)
  // - xparCL * yparCL -> area partitions

  if( zparCL != 0 )
  {
    m_numPartZ = zparCL;
    GEOS_LOG_RANK_0( GEOS_FMT( "LayeredMeshPartitioner: Command-line override - numPartZ set to {}", m_numPartZ ) );
  }

  if( xparCL != 0 || yparCL != 0 )
  {
    unsigned int const xPart = ( xparCL != 0 ? xparCL : 1 );
    unsigned int const yPart = ( yparCL != 0 ? yparCL : 1 );
    unsigned int const areaPartitions = xPart * yPart;
    unsigned int const totalPartitions = areaPartitions * m_numPartZ;

    setPartitionCounts( xPart, yPart, m_numPartZ );

    GEOS_LOG_RANK_0( GEOS_FMT( "LayeredMeshPartitioner: Command-line override - "
                               "partitioning into {} total parts ({} area x {} layers)",
                               totalPartitions, areaPartitions, m_numPartZ ) );
  }
  else
  {
    // Default: partition into commSize parts
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOS );

    GEOS_THROW_IF_NE_MSG( mpiSize % m_numPartZ, 0,
                          GEOS_FMT( "Total MPI ranks ({}) must be evenly divisible by numPartZ ({})", mpiSize, m_numPartZ ), InputError );

    int const areaPartitions = mpiSize / m_numPartZ;
    setPartitionCounts( 1, areaPartitions, m_numPartZ );

    GEOS_LOG_RANK_0( GEOS_FMT( "LayeredMeshPartitioner: Partitioning into {} parts "
                               "({} area x {} layers)",
                               mpiSize, areaPartitions, m_numPartZ ) );
  }
}

array1d< int64_t > LayeredMeshPartitioner::computePartitioning( vtk::AllMeshes & mesh,
                                                                MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  GraphPartitionEngine & engine = getEngine();
  return vtk::partitionByAreaAndLayer( mesh,
                                       m_indexArrayName,
                                       engine,
                                       m_numPartZ,
                                       comm );
}

REGISTER_CATALOG_ENTRY( DomainPartitioner, LayeredMeshPartitioner, string const &, Group * const )

} // namespace geos
