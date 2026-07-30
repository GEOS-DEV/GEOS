/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolidMechanicsMPMRigidBody.cpp
 *
 * MPI and periodic-domain support for the color-partitioned RigidBodyMPM event.
 * The implementation is kept separate from SolidMechanicsMPM.cpp so the main
 * continuum explicit-step algorithm remains readable.
 */

#include "SolidMechanicsMPM.hpp"

#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace geos
{

using namespace dataRepository;

namespace
{

constexpr int rigidBodyRootRank = 0;
constexpr globalIndex rigidBodyFieldKeyEmpty = std::numeric_limits< globalIndex >::max();

struct LocalRigidBodyRegistryRecord
{
  globalIndex minimumParticleID = std::numeric_limits< globalIndex >::max();
  integer minimumContactGroup = std::numeric_limits< integer >::max();
  integer maximumContactGroup = std::numeric_limits< integer >::lowest();
};

using RigidBodyColorIndexMap = std::unordered_map< integer, localIndex >;

RigidBodyColorIndexMap makeRigidBodyColorIndexMap( array1d< integer > const & colors )
{
  RigidBodyColorIndexMap result;
  result.reserve( static_cast< std::size_t >( colors.size() ) );
  for( localIndex bodyIndex = 0; bodyIndex < colors.size(); ++bodyIndex )
  {
    result.emplace( colors[bodyIndex], bodyIndex );
  }
  return result;
}

localIndex findRigidBodyIndex( RigidBodyColorIndexMap const & colorToBodyIndex,
                               integer const color )
{
  auto const iterator = colorToBodyIndex.find( color );
  return iterator == colorToBodyIndex.end() ? -1 : iterator->second;
}

enum RigidBodyAccumulatorIndex : int
{
  rigidMass = 0,
  rigidFirstMomentX = 1,
  rigidFirstMomentY = 2,
  rigidFirstMomentZ = 3,
  rigidMomentumX = 4,
  rigidMomentumY = 5,
  rigidMomentumZ = 6,
  rigidSecondMomentXY = 7,
  rigidAngularMomentumOrigin = 8,
  rigidMappedContactMass = 9,
  rigidContactImpulseX = 10,
  rigidContactImpulseY = 11,
  rigidContactImpulseZ = 12,
  rigidContactTorqueImpulseOrigin = 13,
  rigidBodyAccumulatorStride = 14
};

enum RigidBodySolutionIndex : int
{
  rigidCenterOldX = 0,
  rigidCenterOldY = 1,
  rigidCenterOldZ = 2,
  rigidCenterNewX = 3,
  rigidCenterNewY = 4,
  rigidCenterNewZ = 5,
  rigidVelocityNewX = 6,
  rigidVelocityNewY = 7,
  rigidVelocityNewZ = 8,
  rigidOmegaNewZ = 9,
  rigidBodyMass = 10,
  rigidBodyInertiaZ = 11,
  rigidBodyMappedContactMass = 12,
  rigidBodyContactMassRatio = 13,
  rigidBodyForceX = 14,
  rigidBodyForceY = 15,
  rigidBodyForceZ = 16,
  rigidBodyTorqueZ = 17,
  rigidBodyCosDeltaTheta = 18,
  rigidBodySinDeltaTheta = 19,
  rigidBodySolutionStride = 20
};

/**
 * @brief Deterministically sums one fixed-size buffer on rank zero.
 *
 * MPI_Reduce is sufficient for mathematical correctness, but its reduction tree
 * is implementation dependent. RigidBodyMPM has few bodies, so gathering one
 * compact buffer per rank and summing in rank order is inexpensive and makes
 * rank zero the authoritative solution within a fixed decomposition.
 */
void gatherAndDeterministicallySumToRoot( std::vector< real64 > const & localValues,
                                         std::vector< real64 > & globalValues )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  int const valueCount = static_cast< int >( localValues.size() );

  std::vector< real64 > gatheredValues;
  if( rank == rigidBodyRootRank )
  {
    gatheredValues.resize( static_cast< std::size_t >( commSize ) * localValues.size(), 0.0 );
    globalValues.assign( localValues.size(), 0.0 );
  }

  MPI_Gather( valueCount > 0 ? localValues.data() : nullptr,
              valueCount,
              MPI_DOUBLE,
              rank == rigidBodyRootRank && valueCount > 0 ? gatheredValues.data() : nullptr,
              valueCount,
              MPI_DOUBLE,
              rigidBodyRootRank,
              MPI_COMM_GEOS );

  if( rank == rigidBodyRootRank )
  {
    for( int valueIndex = 0; valueIndex < valueCount; ++valueIndex )
    {
      long double sum = 0.0L;
      long double compensation = 0.0L;
      for( int sourceRank = 0; sourceRank < commSize; ++sourceRank )
      {
        long double const value = static_cast< long double >(
          gatheredValues[static_cast< std::size_t >( sourceRank ) * valueCount + valueIndex] );
        long double const correctedValue = value - compensation;
        long double const nextSum = sum + correctedValue;
        compensation = ( nextSum - sum ) - correctedValue;
        sum = nextSum;
      }
      globalValues[valueIndex] = static_cast< real64 >( sum );
    }
  }
}

GEOS_FORCE_INLINE
real64 minimumImageCoordinate( real64 const wrappedCoordinate,
                               real64 const referenceCoordinate,
                               real64 const period,
                               bool const periodic )
{
  if( !periodic || period <= 0.0 )
  {
    return wrappedCoordinate;
  }

  real64 displacement = wrappedCoordinate - referenceCoordinate;
  displacement -= period * std::round( displacement / period );
  return referenceCoordinate + displacement;
}

GEOS_FORCE_INLINE
real64 wrapCoordinate( real64 const coordinate,
                       real64 const minimum,
                       real64 const period )
{
  if( period <= 0.0 )
  {
    return coordinate;
  }

  real64 wrapped = std::fmod( coordinate - minimum, period );
  if( wrapped < 0.0 )
  {
    wrapped += period;
  }
  return minimum + wrapped;
}

GEOS_FORCE_INLINE
globalIndex makeRigidBodyFieldKey( localIndex const bodyIndex,
                                   integer const color )
{
  static_assert( sizeof( globalIndex ) >= sizeof( std::uint64_t ),
                 "RigidBodyMPM deterministic field keys require a 64-bit globalIndex." );

  std::uint64_t const priority = static_cast< std::uint64_t >( bodyIndex );
  std::uint64_t const unsignedColor = static_cast< std::uint32_t >( color );
  return static_cast< globalIndex >( ( priority << 32 ) | unsignedColor );
}

GEOS_FORCE_INLINE
localIndex rigidBodyIndexFromFieldKey( globalIndex const key )
{
  return static_cast< localIndex >( static_cast< std::uint64_t >( key ) >> 32 );
}

GEOS_FORCE_INLINE
integer rigidBodyColorFromFieldKey( globalIndex const key )
{
  return static_cast< integer >( static_cast< std::uint32_t >( static_cast< std::uint64_t >( key ) ) );
}

} // namespace

void SolidMechanicsMPM::resetRigidBodyMPIState()
{
  m_rigidBodyRegistryInitialized = 0;
  m_rigidBodyRegistrySynchronized = 0;
  m_rigidBodyColors.resize( 0 );
  m_rigidBodyContactGroups.resize( 0 );
  m_rigidBodyAnchorParticleIDs.resize( 0 );
  m_rigidBodyUnwrappedCenters.resize( 0 );
}

/**
 * @brief Returns the MPI-consistent timestep limit used only by RigidBodyMPM.
 *
 * Continuum MPM retains its existing local wavespeed/advection estimate. This
 * helper is called only while a rigid event is active and uses global particle
 * speed and cell-length reductions so every rank requests the same rigid-step
 * limit before the outer event machinery combines solver requests.
 */
real64 SolidMechanicsMPM::computeRigidBodyStableTimeStep(
  real64 const localMaximumParticleSpeed,
  real64 const localMinimumCellLength ) const
{
  real64 const globalMaximumParticleSpeed =
    MpiWrapper::max( localMaximumParticleSpeed );
  real64 const globalMinimumCellLength =
    MpiWrapper::min( localMinimumCellLength );

  real64 dtReturn = globalMaximumParticleSpeed > 1.0e-16
                    ? m_cflFactor * globalMinimumCellLength / globalMaximumParticleSpeed
                    : DBL_MAX;

  real64 maximumBoundarySpeed = 0.0;
  for( int i = 0; i < m_numDims; ++i )
  {
    maximumBoundarySpeed = LvArray::math::max(
      maximumBoundarySpeed,
      LvArray::math::abs( m_domainL[i] * m_xGlobalMin[i] ) );
    maximumBoundarySpeed = LvArray::math::max(
      maximumBoundarySpeed,
      LvArray::math::abs( m_domainL[i] * m_xGlobalMax[i] ) );
  }

  if( maximumBoundarySpeed > 1.0e-16 )
  {
    dtReturn = LvArray::math::min(
      dtReturn,
      m_cflFactor * globalMinimumCellLength / maximumBoundarySpeed );
  }

  if( m_rigidBodyContactCFL > 0.0 )
  {
    real64 const maximumClosureSpeed = LvArray::math::max(
      2.0 * globalMaximumParticleSpeed,
      globalMaximumParticleSpeed + maximumBoundarySpeed );
    if( maximumClosureSpeed > 1.0e-16 )
    {
      dtReturn = LvArray::math::min(
        dtReturn,
        m_rigidBodyContactCFL * globalMinimumCellLength / maximumClosureSpeed );
    }
  }

  if( m_rigidBodyMaxTimeStep > 0.0 )
  {
    dtReturn = LvArray::math::min( dtReturn, m_rigidBodyMaxTimeStep );
  }

  return dtReturn;
}

/**
 * @brief Builds a compact global body registry and periodic unwrapped centers.
 *
 * Body order is deterministic: bodies are sorted by the minimum global particle
 * id belonging to each color, with color as a tie break. The same order is used
 * for nodal field priority, rank-zero accumulation, history output, and the
 * broadcast rigid solution.
 */
void SolidMechanicsMPM::initializeRigidBodyMPIState( ParticleManager & particleManager,
                                                      SpatialPartition & partition )
{
  if( m_rigidBodyRegistrySynchronized == 1 )
  {
    return;
  }

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOS );
  int const registryInitializedMinimum =
    MpiWrapper::min( m_rigidBodyRegistryInitialized );
  int const registryInitializedMaximum =
    MpiWrapper::max( m_rigidBodyRegistryInitialized );

  GEOS_ERROR_IF( registryInitializedMinimum != registryInitializedMaximum,
                 "RigidBodyMPM restart state is inconsistent across MPI ranks." );

  if( registryInitializedMaximum == 1 )
  {
    m_rigidBodyColors.move( LvArray::MemorySpace::host, true );
    m_rigidBodyContactGroups.move( LvArray::MemorySpace::host, true );
    m_rigidBodyAnchorParticleIDs.move( LvArray::MemorySpace::host, true );
    m_rigidBodyUnwrappedCenters.move( LvArray::MemorySpace::host, true );

    int numberOfBodies = rank == rigidBodyRootRank
                         ? static_cast< int >( m_rigidBodyColors.size() )
                         : 0;
    MPI_Bcast( &numberOfBodies,
               1,
               MPI_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );
    GEOS_ERROR_IF( numberOfBodies <= 0,
                   "Restarted RigidBodyMPM registry has no bodies." );

    if( rank != rigidBodyRootRank )
    {
      m_rigidBodyColors.resize( numberOfBodies );
      m_rigidBodyContactGroups.resize( numberOfBodies );
      m_rigidBodyAnchorParticleIDs.resize( numberOfBodies );
      m_rigidBodyUnwrappedCenters.resize( 3 * numberOfBodies );
    }

    std::vector< long long > anchorParticleIDs( numberOfBodies, 0 );
    int invalidRestartRegistry = 0;
    if( rank == rigidBodyRootRank )
    {
      invalidRestartRegistry =
        m_rigidBodyContactGroups.size() != m_rigidBodyColors.size() ||
        m_rigidBodyAnchorParticleIDs.size() != m_rigidBodyColors.size() ||
        m_rigidBodyUnwrappedCenters.size() != 3 * m_rigidBodyColors.size();
      if( invalidRestartRegistry == 0 )
      {
        for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
        {
          anchorParticleIDs[bodyIndex] =
            static_cast< long long >( m_rigidBodyAnchorParticleIDs[bodyIndex] );
        }
      }
    }
    MPI_Bcast( &invalidRestartRegistry,
               1,
               MPI_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );
    GEOS_ERROR_IF( invalidRestartRegistry != 0,
                   "Restarted RigidBodyMPM registry arrays are inconsistent." );

    MPI_Bcast( m_rigidBodyColors.data(),
               numberOfBodies,
               MPI_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );
    MPI_Bcast( m_rigidBodyContactGroups.data(),
               numberOfBodies,
               MPI_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );
    MPI_Bcast( anchorParticleIDs.data(),
               numberOfBodies,
               MPI_LONG_LONG_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );
    MPI_Bcast( m_rigidBodyUnwrappedCenters.data(),
               3 * numberOfBodies,
               MPI_DOUBLE,
               rigidBodyRootRank,
               MPI_COMM_GEOS );

    if( rank != rigidBodyRootRank )
    {
      for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
      {
        m_rigidBodyAnchorParticleIDs[bodyIndex] =
          static_cast< globalIndex >( anchorParticleIDs[bodyIndex] );
      }
    }
    m_rigidBodyRegistrySynchronized = 1;
    return;
  }

  std::map< integer, LocalRigidBodyRegistryRecord > localRegistry;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.getWrapperBase( fields::mpm::particleColor::key() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleIDString() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleGroupString() ).move( LvArray::MemorySpace::host, true );

    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    arrayView1d< localIndex const > const particleGroup = subRegion.getParticleGroup();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      integer const color = particleColor[p];
      GEOS_ERROR_IF( color < 0, "RigidBodyMPM particleColor values must be non-negative." );

      LocalRigidBodyRegistryRecord & record = localRegistry[color];
      record.minimumParticleID = std::min( record.minimumParticleID, particleID[p] );
      record.minimumContactGroup = std::min( record.minimumContactGroup,
                                             static_cast< integer >( particleGroup[p] ) );
      record.maximumContactGroup = std::max( record.maximumContactGroup,
                                             static_cast< integer >( particleGroup[p] ) );
    }
  } );

  std::vector< long long > localRecords;
  localRecords.reserve( 4 * localRegistry.size() );
  for( auto const & entry : localRegistry )
  {
    localRecords.push_back( static_cast< long long >( entry.first ) );
    localRecords.push_back( static_cast< long long >( entry.second.minimumParticleID ) );
    localRecords.push_back( static_cast< long long >( entry.second.minimumContactGroup ) );
    localRecords.push_back( static_cast< long long >( entry.second.maximumContactGroup ) );
  }

  int const localRecordCount = static_cast< int >( localRegistry.size() );
  std::vector< int > recordCounts;
  if( rank == rigidBodyRootRank )
  {
    recordCounts.resize( commSize, 0 );
  }
  MPI_Gather( &localRecordCount,
              1,
              MPI_INT,
              rank == rigidBodyRootRank ? recordCounts.data() : nullptr,
              1,
              MPI_INT,
              rigidBodyRootRank,
              MPI_COMM_GEOS );

  std::vector< int > valueCounts;
  std::vector< int > valueDisplacements;
  std::vector< long long > gatheredRecords;
  if( rank == rigidBodyRootRank )
  {
    valueCounts.resize( commSize, 0 );
    valueDisplacements.resize( commSize, 0 );
    int totalRecords = 0;
    for( int sourceRank = 0; sourceRank < commSize; ++sourceRank )
    {
      valueDisplacements[sourceRank] = 4 * totalRecords;
      valueCounts[sourceRank] = 4 * recordCounts[sourceRank];
      totalRecords += recordCounts[sourceRank];
    }
    gatheredRecords.resize( 4 * totalRecords );
  }

  MPI_Gatherv( localRecords.empty() ? nullptr : localRecords.data(),
               static_cast< int >( localRecords.size() ),
               MPI_LONG_LONG_INT,
               rank == rigidBodyRootRank && !gatheredRecords.empty() ? gatheredRecords.data() : nullptr,
               rank == rigidBodyRootRank ? valueCounts.data() : nullptr,
               rank == rigidBodyRootRank ? valueDisplacements.data() : nullptr,
               MPI_LONG_LONG_INT,
               rigidBodyRootRank,
               MPI_COMM_GEOS );

  std::vector< integer > colors;
  std::vector< integer > contactGroups;
  std::vector< long long > anchorParticleIDs;
  integer invalidContactGroupColor = -1;

  if( rank == rigidBodyRootRank )
  {
    std::map< integer, LocalRigidBodyRegistryRecord > globalRegistry;
    for( std::size_t offset = 0; offset + 3 < gatheredRecords.size(); offset += 4 )
    {
      integer const color = static_cast< integer >( gatheredRecords[offset] );
      globalIndex const minimumParticleID = static_cast< globalIndex >( gatheredRecords[offset + 1] );
      integer const minimumContactGroup = static_cast< integer >( gatheredRecords[offset + 2] );
      integer const maximumContactGroup = static_cast< integer >( gatheredRecords[offset + 3] );

      LocalRigidBodyRegistryRecord & record = globalRegistry[color];
      record.minimumParticleID = std::min( record.minimumParticleID, minimumParticleID );
      record.minimumContactGroup = std::min( record.minimumContactGroup, minimumContactGroup );
      record.maximumContactGroup = std::max( record.maximumContactGroup, maximumContactGroup );
    }

    std::vector< std::tuple< globalIndex, integer, integer > > orderedBodies;
    orderedBodies.reserve( globalRegistry.size() );
    for( auto const & entry : globalRegistry )
    {
      if( entry.second.minimumContactGroup != entry.second.maximumContactGroup )
      {
        invalidContactGroupColor = entry.first;
      }
      orderedBodies.emplace_back( entry.second.minimumParticleID,
                                  entry.first,
                                  entry.second.minimumContactGroup );
    }
    std::sort( orderedBodies.begin(), orderedBodies.end() );

    colors.resize( orderedBodies.size() );
    contactGroups.resize( orderedBodies.size() );
    anchorParticleIDs.resize( orderedBodies.size() );
    for( std::size_t bodyIndex = 0; bodyIndex < orderedBodies.size(); ++bodyIndex )
    {
      anchorParticleIDs[bodyIndex] = static_cast< long long >( std::get< 0 >( orderedBodies[bodyIndex] ) );
      colors[bodyIndex] = std::get< 1 >( orderedBodies[bodyIndex] );
      contactGroups[bodyIndex] = std::get< 2 >( orderedBodies[bodyIndex] );
    }
  }

  MPI_Bcast( &invalidContactGroupColor,
             1,
             MPI_INT,
             rigidBodyRootRank,
             MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidContactGroupColor >= 0,
                 "RigidBodyMPM requires a uniform contact group within each particleColor. Color "
                 << invalidContactGroupColor << " contains multiple contact groups. Any number of "
                 "particleColor bodies may share the same catch-all contact group." );

  int numberOfBodies = rank == rigidBodyRootRank ? static_cast< int >( colors.size() ) : 0;
  MPI_Bcast( &numberOfBodies,
             1,
             MPI_INT,
             rigidBodyRootRank,
             MPI_COMM_GEOS );
  GEOS_ERROR_IF( numberOfBodies <= 0,
                 "RigidBodyMPM did not find any active particleColor bodies." );

  if( rank != rigidBodyRootRank )
  {
    colors.resize( numberOfBodies );
    contactGroups.resize( numberOfBodies );
    anchorParticleIDs.resize( numberOfBodies );
  }
  MPI_Bcast( colors.data(), numberOfBodies, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  MPI_Bcast( contactGroups.data(), numberOfBodies, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  MPI_Bcast( anchorParticleIDs.data(), numberOfBodies, MPI_LONG_LONG_INT, rigidBodyRootRank, MPI_COMM_GEOS );

  m_rigidBodyColors.resize( numberOfBodies );
  m_rigidBodyContactGroups.resize( numberOfBodies );
  m_rigidBodyAnchorParticleIDs.resize( numberOfBodies );
  for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
  {
    m_rigidBodyColors[bodyIndex] = colors[bodyIndex];
    m_rigidBodyContactGroups[bodyIndex] = contactGroups[bodyIndex];
    m_rigidBodyAnchorParticleIDs[bodyIndex] = static_cast< globalIndex >( anchorParticleIDs[bodyIndex] );
  }

  RigidBodyColorIndexMap const colorToBodyIndex =
    makeRigidBodyColorIndexMap( m_rigidBodyColors );

  std::vector< real64 > anchorLocal( 4 * numberOfBodies, 0.0 );
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleCenterString() ).move( LvArray::MemorySpace::host, true );
    arrayView1d< globalIndex const > const particleID = subRegion.getParticleID();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
      {
        if( particleID[p] == m_rigidBodyAnchorParticleIDs[bodyIndex] )
        {
          anchorLocal[4 * bodyIndex] += 1.0;
          for( int i = 0; i < 3; ++i )
          {
            anchorLocal[4 * bodyIndex + 1 + i] += particlePosition[p][i];
          }
        }
      }
    }
  } );

  std::vector< real64 > anchorGlobal;
  gatherAndDeterministicallySumToRoot( anchorLocal, anchorGlobal );
  if( rank != rigidBodyRootRank )
  {
    anchorGlobal.resize( anchorLocal.size(), 0.0 );
  }
  MPI_Bcast( anchorGlobal.data(),
             static_cast< int >( anchorGlobal.size() ),
             MPI_DOUBLE,
             rigidBodyRootRank,
             MPI_COMM_GEOS );

  integer invalidAnchorBody = -1;
  if( rank == rigidBodyRootRank )
  {
    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      if( LvArray::math::abs( anchorGlobal[4 * bodyIndex] - 1.0 ) > 0.5 )
      {
        invalidAnchorBody = bodyIndex;
      }
    }
  }
  MPI_Bcast( &invalidAnchorBody,
             1,
             MPI_INT,
             rigidBodyRootRank,
             MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidAnchorBody >= 0,
                 "RigidBodyMPM could not identify exactly one master anchor particle for color "
                 << m_rigidBodyColors[invalidAnchorBody] << "." );

  arrayView1d< int const > const periodic = partition.getPeriodic();
  real64 period[3] = { m_domainExtent[0], m_domainExtent[1], m_domainExtent[2] };
  std::vector< real64 > centerAccumulatorLocal( 4 * numberOfBodies, 0.0 );

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.getWrapperBase( fields::mpm::particleMass::key() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleSurfaceFlagString() ).move( LvArray::MemorySpace::host, true );

    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< integer const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      if( !isRigidBodySurfaceFlag( particleSurfaceFlag[p] ) )
      {
        continue;
      }

      localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, particleColor[p] );
      GEOS_ERROR_IF( bodyIndex < 0, "RigidBodyMPM body registry is missing an active particle color." );
      real64 const mass = particleMass[p];
      centerAccumulatorLocal[4 * bodyIndex] += mass;
      for( int i = 0; i < 3; ++i )
      {
        real64 const anchorCoordinate = anchorGlobal[4 * bodyIndex + 1 + i];
        real64 const unwrappedCoordinate = minimumImageCoordinate( particlePosition[p][i],
                                                                   anchorCoordinate,
                                                                   period[i],
                                                                   periodic[i] == 1 );
        centerAccumulatorLocal[4 * bodyIndex + 1 + i] += mass * unwrappedCoordinate;
      }
    }
  } );

  std::vector< real64 > centerAccumulatorGlobal;
  gatherAndDeterministicallySumToRoot( centerAccumulatorLocal, centerAccumulatorGlobal );
  std::vector< real64 > initialCenters( 3 * numberOfBodies, 0.0 );
  integer invalidSurfaceBody = -1;
  if( rank == rigidBodyRootRank )
  {
    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      real64 const mass = centerAccumulatorGlobal[4 * bodyIndex];
      if( mass <= 0.0 )
      {
        invalidSurfaceBody = bodyIndex;
        continue;
      }
      for( int i = 0; i < 3; ++i )
      {
        initialCenters[3 * bodyIndex + i] = centerAccumulatorGlobal[4 * bodyIndex + 1 + i] / mass;
      }
    }
  }
  MPI_Bcast( &invalidSurfaceBody,
             1,
             MPI_INT,
             rigidBodyRootRank,
             MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidSurfaceBody >= 0,
                 "RigidBodyMPM color " << m_rigidBodyColors[invalidSurfaceBody]
                 << " has no surface-shell mass." );
  MPI_Bcast( initialCenters.data(),
             static_cast< int >( initialCenters.size() ),
             MPI_DOUBLE,
             rigidBodyRootRank,
             MPI_COMM_GEOS );

  m_rigidBodyUnwrappedCenters.resize( initialCenters.size() );
  for( localIndex i = 0; i < m_rigidBodyUnwrappedCenters.size(); ++i )
  {
    m_rigidBodyUnwrappedCenters[i] = initialCenters[i];
  }

  std::vector< real64 > extentLocal( 3 * numberOfBodies, 0.0 );
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, particleColor[p] );
      GEOS_ERROR_IF( bodyIndex < 0,
                     "RigidBodyMPM body registry is missing an active particle color." );
      for( int i = 0; i < m_numDims; ++i )
      {
        real64 const coordinate = minimumImageCoordinate( particlePosition[p][i],
                                                          initialCenters[3 * bodyIndex + i],
                                                          period[i],
                                                          periodic[i] == 1 );
        extentLocal[3 * bodyIndex + i] = std::max( extentLocal[3 * bodyIndex + i],
                                                   LvArray::math::abs( coordinate - initialCenters[3 * bodyIndex + i] ) );
      }
    }
  } );

  std::vector< real64 > extentGlobal( extentLocal.size(), 0.0 );
  MPI_Allreduce( extentLocal.data(),
                 extentGlobal.data(),
                 static_cast< int >( extentLocal.size() ),
                 MPI_DOUBLE,
                 MPI_MAX,
                 MPI_COMM_GEOS );

  integer invalidPeriodicBody = -1;
  integer invalidPeriodicDirection = -1;
  if( rank == rigidBodyRootRank )
  {
    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      for( int i = 0; i < m_numDims; ++i )
      {
        if( periodic[i] == 1 && 2.0 * extentGlobal[3 * bodyIndex + i] >= 0.98 * period[i] )
        {
          invalidPeriodicBody = bodyIndex;
          invalidPeriodicDirection = i;
        }
      }
    }
  }
  MPI_Bcast( &invalidPeriodicBody, 1, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  MPI_Bcast( &invalidPeriodicDirection, 1, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidPeriodicBody >= 0,
                 "RigidBodyMPM color " << m_rigidBodyColors[invalidPeriodicBody]
                 << " reaches the half-period minimum-image limit in direction "
                 << invalidPeriodicDirection
                 << "; a unique unwrapped rigid body cannot be constructed." );

  m_rigidBodyRegistryInitialized = 1;
  m_rigidBodyRegistrySynchronized = 1;
}

/**
 * @brief Assigns globally deterministic color-defined nodal fields.
 *
 * The first pure fields are selected by repeated shared-node MPI_MIN reductions
 * of a key ordered by global body priority. This avoids rank-local particle
 * traversal order and is valid at ordinary and periodic partition interfaces.
 * The last field is the weighted overflow field.
 */
void SolidMechanicsMPM::computeRigidBodyColorFieldMappings( DomainPartition & domain,
                                                            ParticleManager & particleManager,
                                                            NodeManager & nodeManager,
                                                            MeshLevel & mesh )
{
  GEOS_MARK_FUNCTION;

  SpatialPartition & partition =
    dynamic_cast< SpatialPartition & >( domain.getGroup( domain.groupKeys.partitionManager ) );
  initializeRigidBodyMPIState( particleManager, partition );

  RigidBodyColorIndexMap const colorToBodyIndex =
    makeRigidBodyColorIndexMap( m_rigidBodyColors );

  int const maxGridFields = m_rigidBodyMaxGridFields > 0 ? m_rigidBodyMaxGridFields : m_numVelocityFields;
  GEOS_ERROR_IF( maxGridFields < 2,
                 "RigidBodyMPM requires maxGridFields >= 2 so at least one pure field and one overflow field exist." );
  GEOS_ERROR_IF( maxGridFields > m_numVelocityFields,
                 "RigidBodyMPM maxGridFields exceeds the number of allocated MPM velocity fields." );

  localIndex const numVelocityFields = m_numVelocityFields;

  // The deterministic field-selection pass is host based because shared-node
  // reductions are currently packed/unpacked on the host. Explicitly stage all
  // particle labels and map arrays needed by that pass; later P2G code moves
  // the effective maps back to the selected execution space.
  localIndex stagingSubRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.getWrapperBase( fields::mpm::particleColor::key() ).move(
      LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase(
      ParticleSubRegion::viewKeyStruct::particleSurfaceFlagString() ).move(
      LvArray::MemorySpace::host, true );

    m_mappedNodes[stagingSubRegionIndex].move( LvArray::MemorySpace::host, true );
    m_mappedFields[stagingSubRegionIndex].move( LvArray::MemorySpace::host, true );
    m_numEffectiveMappedNodes[stagingSubRegionIndex].move( LvArray::MemorySpace::host, true );
    m_effectiveMappedNodes[stagingSubRegionIndex].move( LvArray::MemorySpace::host, true );
    m_effectiveMappedFields[stagingSubRegionIndex].move( LvArray::MemorySpace::host, true );
    ++stagingSubRegionIndex;
  } );

  WrapperBase & fieldKeyWrapper = nodeManager.getWrapperBase( viewKeyStruct::gridRigidBodyFieldKeyString() );
  WrapperBase & fieldColorWrapper = nodeManager.getWrapperBase( viewKeyStruct::gridRigidBodyFieldColorString() );
  WrapperBase & fieldGroupWrapper = nodeManager.getWrapperBase( viewKeyStruct::gridRigidBodyFieldContactGroupString() );
  fieldKeyWrapper.move( LvArray::MemorySpace::host, true );
  fieldColorWrapper.move( LvArray::MemorySpace::host, true );
  fieldGroupWrapper.move( LvArray::MemorySpace::host, true );

  arrayView2d< globalIndex > const nodeFieldKey =
    nodeManager.getReference< array2d< globalIndex > >( viewKeyStruct::gridRigidBodyFieldKeyString() );
  arrayView2d< integer > const nodeFieldColor =
    nodeManager.getReference< array2d< integer > >( viewKeyStruct::gridRigidBodyFieldColorString() );
  arrayView2d< integer > const nodeFieldContactGroup =
    nodeManager.getReference< array2d< integer > >( viewKeyStruct::gridRigidBodyFieldContactGroupString() );

  for( localIndex nodeIndex = 0; nodeIndex < nodeManager.size(); ++nodeIndex )
  {
    for( localIndex fieldIndex = 0; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      nodeFieldKey[nodeIndex][fieldIndex] = rigidBodyFieldKeyEmpty;
      nodeFieldColor[nodeIndex][fieldIndex] = -1;
      nodeFieldContactGroup[nodeIndex][fieldIndex] = -1;
    }
  }

  int const numberOfPureFields = maxGridFields - 1;
  for( int targetField = 0; targetField < numberOfPureFields; ++targetField )
  {
    for( localIndex nodeIndex = 0; nodeIndex < nodeManager.size(); ++nodeIndex )
    {
      nodeFieldKey[nodeIndex][targetField] = rigidBodyFieldKeyEmpty;
    }

    localIndex subRegionIndex = 0;
    particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[subRegionIndex];
      arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
      arrayView1d< integer const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
      int const mappedNodeCount = 8 * subRegion.numberOfVerticesPerParticle();
      SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

      for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
      {
        localIndex const p = activeParticleIndices[pp];
        if( !isRigidBodySurfaceFlag( particleSurfaceFlag[p] ) )
        {
          continue;
        }

        integer const color = particleColor[p];
        localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, color );
        GEOS_ERROR_IF( bodyIndex < 0, "RigidBodyMPM body registry is missing an active color." );
        globalIndex const candidateKey = makeRigidBodyFieldKey( bodyIndex, color );

        for( int mappedIndex = 0; mappedIndex < mappedNodeCount; ++mappedIndex )
        {
          localIndex const nodeIndex = mappedNodes[pp][mappedIndex];
          bool alreadySelected = false;
          for( int fieldIndex = 0; fieldIndex < targetField; ++fieldIndex )
          {
            alreadySelected |= nodeFieldColor[nodeIndex][fieldIndex] == color;
          }
          if( !alreadySelected )
          {
            nodeFieldKey[nodeIndex][targetField] =
              std::min( nodeFieldKey[nodeIndex][targetField], candidateKey );
          }
        }
      }
      ++subRegionIndex;
    } );

    syncGridFields( { viewKeyStruct::gridRigidBodyFieldKeyString() },
                    domain,
                    nodeManager,
                    mesh,
                    MPI_MIN );

    for( localIndex nodeIndex = 0; nodeIndex < nodeManager.size(); ++nodeIndex )
    {
      globalIndex const key = nodeFieldKey[nodeIndex][targetField];
      nodeFieldColor[nodeIndex][targetField] =
        key == rigidBodyFieldKeyEmpty ? -1 : rigidBodyColorFromFieldKey( key );
    }
  }

  int const overflowField = maxGridFields - 1;
  for( localIndex nodeIndex = 0; nodeIndex < nodeManager.size(); ++nodeIndex )
  {
    nodeFieldKey[nodeIndex][overflowField] = 0;

    // Temporarily use the overflow entries of the two metadata arrays to
    // accumulate the minimum and maximum contact group among overflow colors.
    // This preserves the common one-catch-all-group case exactly while still
    // detecting the rare mixed-group overflow node.
    nodeFieldContactGroup[nodeIndex][overflowField] =
      std::numeric_limits< integer >::max();
    nodeFieldColor[nodeIndex][overflowField] =
      std::numeric_limits< integer >::lowest();
  }

  localIndex overflowSubRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[overflowSubRegionIndex];
    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView1d< integer const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    int const mappedNodeCount = 8 * subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      if( !isRigidBodySurfaceFlag( particleSurfaceFlag[p] ) )
      {
        continue;
      }
      integer const color = particleColor[p];
      localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, color );
      GEOS_ERROR_IF( bodyIndex < 0,
                     "RigidBodyMPM body registry is missing an overflow-field color." );
      integer const contactGroup = m_rigidBodyContactGroups[bodyIndex];

      for( int mappedIndex = 0; mappedIndex < mappedNodeCount; ++mappedIndex )
      {
        localIndex const nodeIndex = mappedNodes[pp][mappedIndex];
        bool selected = false;
        for( int fieldIndex = 0; fieldIndex < numberOfPureFields; ++fieldIndex )
        {
          selected |= nodeFieldColor[nodeIndex][fieldIndex] == color;
        }
        if( !selected )
        {
          nodeFieldKey[nodeIndex][overflowField] = 1;
          nodeFieldContactGroup[nodeIndex][overflowField] =
            std::min( nodeFieldContactGroup[nodeIndex][overflowField], contactGroup );
          nodeFieldColor[nodeIndex][overflowField] =
            std::max( nodeFieldColor[nodeIndex][overflowField], contactGroup );
        }
      }
    }
    ++overflowSubRegionIndex;
  } );

  syncGridFields( { viewKeyStruct::gridRigidBodyFieldKeyString() },
                  domain,
                  nodeManager,
                  mesh,
                  MPI_MAX );
  syncGridFields( { viewKeyStruct::gridRigidBodyFieldContactGroupString() },
                  domain,
                  nodeManager,
                  mesh,
                  MPI_MIN );
  syncGridFields( { viewKeyStruct::gridRigidBodyFieldColorString() },
                  domain,
                  nodeManager,
                  mesh,
                  MPI_MAX );

  for( localIndex nodeIndex = 0; nodeIndex < nodeManager.size(); ++nodeIndex )
  {
    for( int fieldIndex = 0; fieldIndex < numberOfPureFields; ++fieldIndex )
    {
      globalIndex const key = nodeFieldKey[nodeIndex][fieldIndex];
      if( key == rigidBodyFieldKeyEmpty )
      {
        nodeFieldColor[nodeIndex][fieldIndex] = -1;
        nodeFieldContactGroup[nodeIndex][fieldIndex] = -1;
      }
      else
      {
        localIndex const bodyIndex = rigidBodyIndexFromFieldKey( key );
        nodeFieldColor[nodeIndex][fieldIndex] = rigidBodyColorFromFieldKey( key );
        nodeFieldContactGroup[nodeIndex][fieldIndex] = m_rigidBodyContactGroups[bodyIndex];
      }
    }

    bool const overflowActive = nodeFieldKey[nodeIndex][overflowField] != 0;
    if( overflowActive )
    {
      integer const minimumOverflowGroup =
        nodeFieldContactGroup[nodeIndex][overflowField];
      integer const maximumOverflowGroup =
        nodeFieldColor[nodeIndex][overflowField];
      nodeFieldColor[nodeIndex][overflowField] = -2;
      nodeFieldContactGroup[nodeIndex][overflowField] =
        minimumOverflowGroup == maximumOverflowGroup
        ? minimumOverflowGroup
        : -2;
    }
    else
    {
      nodeFieldColor[nodeIndex][overflowField] = -1;
      nodeFieldContactGroup[nodeIndex][overflowField] = -1;
    }

    for( localIndex fieldIndex = maxGridFields; fieldIndex < numVelocityFields; ++fieldIndex )
    {
      nodeFieldColor[nodeIndex][fieldIndex] = -1;
      nodeFieldContactGroup[nodeIndex][fieldIndex] = -1;
    }
  }

  localIndex assignmentSubRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    arrayView2d< integer > const mappedFields = m_mappedFields[assignmentSubRegionIndex];
    arrayView2d< localIndex const > const mappedNodes = m_mappedNodes[assignmentSubRegionIndex];
    arrayView1d< localIndex const > const numEffectiveMappedNodes = m_numEffectiveMappedNodes[assignmentSubRegionIndex];
    arrayView2d< integer > const effectiveMappedFields = m_effectiveMappedFields[assignmentSubRegionIndex];
    arrayView2d< localIndex const > const effectiveMappedNodes = m_effectiveMappedNodes[assignmentSubRegionIndex];
    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView1d< integer const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    int const mappedNodeCount = 8 * subRegion.numberOfVerticesPerParticle();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      bool const surfaceParticle = isRigidBodySurfaceFlag( particleSurfaceFlag[p] );
      integer const color = particleColor[p];

      for( int mappedIndex = 0; mappedIndex < mappedNodeCount; ++mappedIndex )
      {
        if( !surfaceParticle )
        {
          mappedFields[pp][mappedIndex] = 0;
          continue;
        }

        localIndex const nodeIndex = mappedNodes[pp][mappedIndex];
        integer fieldIndex = overflowField;
        for( int candidateField = 0; candidateField < numberOfPureFields; ++candidateField )
        {
          if( nodeFieldColor[nodeIndex][candidateField] == color )
          {
            fieldIndex = candidateField;
            break;
          }
        }
        mappedFields[pp][mappedIndex] = fieldIndex;
      }

      for( localIndex mappedIndex = 0; mappedIndex < numEffectiveMappedNodes[pp]; ++mappedIndex )
      {
        if( !surfaceParticle )
        {
          effectiveMappedFields[pp][mappedIndex] = 0;
          continue;
        }

        localIndex const nodeIndex = effectiveMappedNodes[pp][mappedIndex];
        integer fieldIndex = overflowField;
        for( int candidateField = 0; candidateField < numberOfPureFields; ++candidateField )
        {
          if( nodeFieldColor[nodeIndex][candidateField] == color )
          {
            fieldIndex = candidateField;
            break;
          }
        }
        effectiveMappedFields[pp][mappedIndex] = fieldIndex;
      }
    }
    ++assignmentSubRegionIndex;
  } );

#ifdef GEOS_USE_DEVICE
  fieldColorWrapper.move( parallelDeviceMemorySpace, true );
  fieldGroupWrapper.move( parallelDeviceMemorySpace, true );

  localIndex deviceSubRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Registry construction and deterministic field selection are host-side.
    // Restore every particle/map array read by the immediately following
    // device P2G kernel, not only the field-number arrays.
    stdVector< string > const rigidMappingParticleFields =
    {
      fields::mpm::particleColor::key(),
      fields::mpm::particleMass::key(),
      ParticleSubRegion::viewKeyStruct::particleIDString(),
      ParticleSubRegion::viewKeyStruct::particleGroupString(),
      ParticleSubRegion::viewKeyStruct::particleCenterString(),
      ParticleSubRegion::viewKeyStruct::particleSurfaceFlagString()
    };
    for( string const & fieldName : rigidMappingParticleFields )
    {
      subRegion.getWrapperBase( fieldName ).move( parallelDeviceMemorySpace, true );
    }

    m_mappedNodes[deviceSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_mappedFields[deviceSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_numEffectiveMappedNodes[deviceSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_effectiveMappedNodes[deviceSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_effectiveMappedFields[deviceSubRegionIndex].move( parallelDeviceMemorySpace, true );
    ++deviceSubRegionIndex;
  } );
#endif
}

/**
 * @brief Solves the surface-shell rigid bodies on rank zero and broadcasts the result.
 *
 * Every rank accumulates raw moments and contact impulses from its master
 * particles. The compact per-body buffers are gathered and deterministically
 * summed on rank zero. Rank zero advances the body state once; all ranks then
 * apply the same broadcast center, velocity, rotation, and diagnostics.
 */
void SolidMechanicsMPM::rigidBodyParticleUpdate( real64 const time_n,
                                                 real64 const dt,
                                                 ParticleManager & particleManager,
                                                 NodeManager & nodeManager,
                                                 SpatialPartition & partition )
{
  GEOS_MARK_FUNCTION;
  GEOS_ERROR_IF( m_planeStrain != 1,
                 "RigidBodyMPM currently implements planar rigid-body updates for 2D plane-strain MPM tests." );
  GEOS_ERROR_IF( dt <= 0.0, "RigidBodyMPM requires a positive time step." );

  initializeRigidBodyMPIState( particleManager, partition );
  RigidBodyColorIndexMap const colorToBodyIndex =
    makeRigidBodyColorIndexMap( m_rigidBodyColors );
  localIndex const numberOfBodies = m_rigidBodyColors.size();
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  arrayView1d< int const > const periodic = partition.getPeriodic();

  nodeManager.getWrapperBase( viewKeyStruct::gridMassString() ).move( LvArray::MemorySpace::host, true );
  nodeManager.getWrapperBase( viewKeyStruct::gridUncontactedVelocityString() ).move( LvArray::MemorySpace::host, true );
  nodeManager.getWrapperBase( viewKeyStruct::gridVelocityString() ).move( LvArray::MemorySpace::host, true );

  arrayView2d< real64 const > const gridMass =
    nodeManager.getReference< array2d< real64 > >( viewKeyStruct::gridMassString() );
  arrayView3d< real64 const > const gridUncontactedVelocity =
    nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridUncontactedVelocityString() );
  arrayView3d< real64 const > const gridVelocity =
    nodeManager.getReference< array3d< real64 > >( viewKeyStruct::gridVelocityString() );

  for( std::size_t subRegionIndex = 0; subRegionIndex < m_numEffectiveMappedNodes.size(); ++subRegionIndex )
  {
    m_numEffectiveMappedNodes[subRegionIndex].move( LvArray::MemorySpace::host, true );
    m_effectiveMappedFields[subRegionIndex].move( LvArray::MemorySpace::host, true );
    m_effectiveMappedNodes[subRegionIndex].move( LvArray::MemorySpace::host, true );
    m_effectiveShapeFunctionValues[subRegionIndex].move( LvArray::MemorySpace::host, true );
  }

  std::vector< real64 > localAccumulator( rigidBodyAccumulatorStride * numberOfBodies, 0.0 );
  // Per body: directional extents x/y/z followed by the planar radial extent.
  std::vector< real64 > localPeriodicExtent( 4 * numberOfBodies, 0.0 );
  real64 const period[3] = { m_domainExtent[0], m_domainExtent[1], m_domainExtent[2] };

  localIndex subRegionIndex = 0;
  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    subRegion.getWrapperBase( fields::mpm::particleColor::key() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( fields::mpm::particleMass::key() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleSurfaceFlagString() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleCenterString() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleVelocityString() ).move( LvArray::MemorySpace::host, true );
    subRegion.getWrapperBase( ParticleSubRegion::viewKeyStruct::particleSurfacePositionString() ).move( LvArray::MemorySpace::host, true );

    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView1d< real64 const > const particleMass = subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< integer const > const particleSurfaceFlag = subRegion.getParticleSurfaceFlag();
    arrayView2d< real64 const > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 const > const particleVelocity = subRegion.getParticleVelocity();
    arrayView2d< real64 const > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView1d< localIndex const > const numEffectiveMappedNodes = m_numEffectiveMappedNodes[subRegionIndex];
    arrayView2d< integer const > const effectiveMappedFields = m_effectiveMappedFields[subRegionIndex];
    arrayView2d< localIndex const > const effectiveMappedNodes = m_effectiveMappedNodes[subRegionIndex];
    arrayView2d< real64 const > const effectiveShapeFunctionValues = m_effectiveShapeFunctionValues[subRegionIndex];
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      if( !isRigidBodySurfaceFlag( particleSurfaceFlag[p] ) )
      {
        continue;
      }

      localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, particleColor[p] );
      GEOS_ERROR_IF( bodyIndex < 0, "RigidBodyMPM body registry is missing an active color." );
      int const offset = rigidBodyAccumulatorStride * bodyIndex;
      real64 x[3] = {};
      for( int i = 0; i < 3; ++i )
      {
        x[i] = minimumImageCoordinate( particlePosition[p][i],
                                       m_rigidBodyUnwrappedCenters[3 * bodyIndex + i],
                                       period[i],
                                       periodic[i] == 1 );
      }

      real64 radialDistanceSquared = 0.0;
      for( int i = 0; i < m_numDims; ++i )
      {
        real64 const displacement =
          x[i] - m_rigidBodyUnwrappedCenters[3 * bodyIndex + i];
        radialDistanceSquared += displacement * displacement;
        if( periodic[i] == 1 )
        {
          localPeriodicExtent[4 * bodyIndex + i] = LvArray::math::max(
            localPeriodicExtent[4 * bodyIndex + i],
            LvArray::math::abs( displacement ) );
        }
      }
      localPeriodicExtent[4 * bodyIndex + 3] = LvArray::math::max(
        localPeriodicExtent[4 * bodyIndex + 3],
        LvArray::math::sqrt( radialDistanceSquared ) );

      real64 const mass = particleMass[p];
      localAccumulator[offset + rigidMass] += mass;
      localAccumulator[offset + rigidFirstMomentX] += mass * x[0];
      localAccumulator[offset + rigidFirstMomentY] += mass * x[1];
      localAccumulator[offset + rigidFirstMomentZ] += mass * x[2];
      localAccumulator[offset + rigidMomentumX] += mass * particleVelocity[p][0];
      localAccumulator[offset + rigidMomentumY] += mass * particleVelocity[p][1];
      localAccumulator[offset + rigidMomentumZ] += mass * particleVelocity[p][2];
      localAccumulator[offset + rigidSecondMomentXY] += mass * ( x[0] * x[0] + x[1] * x[1] );
      localAccumulator[offset + rigidAngularMomentumOrigin] +=
        mass * ( x[0] * particleVelocity[p][1] - x[1] * particleVelocity[p][0] );

      real64 const surfaceX = x[0] + particleSurfacePosition[p][0];
      real64 const surfaceY = x[1] + particleSurfacePosition[p][1];
      for( localIndex mappedIndex = 0; mappedIndex < numEffectiveMappedNodes[pp]; ++mappedIndex )
      {
        localIndex const nodeIndex = effectiveMappedNodes[pp][mappedIndex];
        integer const fieldIndex = effectiveMappedFields[pp][mappedIndex];
        real64 const nodalMass = gridMass[nodeIndex][fieldIndex];
        if( nodalMass <= 0.0 )
        {
          continue;
        }

        real64 deltaVelocity[3] = {};
        real64 deltaVelocityNormSquared = 0.0;
        for( int i = 0; i < 3; ++i )
        {
          deltaVelocity[i] = gridVelocity[nodeIndex][fieldIndex][i]
                             - gridUncontactedVelocity[nodeIndex][fieldIndex][i];
          deltaVelocityNormSquared += deltaVelocity[i] * deltaVelocity[i];
        }
        if( deltaVelocityNormSquared <= 1.0e-30 )
        {
          continue;
        }

        real64 const mappedMass = mass * effectiveShapeFunctionValues[pp][mappedIndex];
        real64 const impulseX = mappedMass * deltaVelocity[0];
        real64 const impulseY = mappedMass * deltaVelocity[1];
        real64 const impulseZ = mappedMass * deltaVelocity[2];
        localAccumulator[offset + rigidMappedContactMass] += mappedMass;
        localAccumulator[offset + rigidContactImpulseX] += impulseX;
        localAccumulator[offset + rigidContactImpulseY] += impulseY;
        localAccumulator[offset + rigidContactImpulseZ] += impulseZ;
        localAccumulator[offset + rigidContactTorqueImpulseOrigin] +=
          surfaceX * impulseY - surfaceY * impulseX;
      }
    }
    ++subRegionIndex;
  } );

#ifdef GEOS_USE_DEVICE
  // The compact body accumulation reads the coalesced map on the host. Restore
  // the map arrays so the next mapping/P2G cycle cannot inherit host-only
  // storage when the active-particle count is unchanged.
  for( std::size_t mapSubRegionIndex = 0;
       mapSubRegionIndex < m_numEffectiveMappedNodes.size();
       ++mapSubRegionIndex )
  {
    m_numEffectiveMappedNodes[mapSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_effectiveMappedFields[mapSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_effectiveMappedNodes[mapSubRegionIndex].move( parallelDeviceMemorySpace, true );
    m_effectiveShapeFunctionValues[mapSubRegionIndex].move( parallelDeviceMemorySpace, true );
  }
#endif

  std::vector< real64 > globalPeriodicExtent( localPeriodicExtent.size(), 0.0 );
  MPI_Allreduce( localPeriodicExtent.data(),
                 globalPeriodicExtent.data(),
                 static_cast< int >( localPeriodicExtent.size() ),
                 MPI_DOUBLE,
                 MPI_MAX,
                 MPI_COMM_GEOS );

  integer invalidPeriodicBody = -1;
  integer invalidPeriodicDirection = -1;
  if( rank == rigidBodyRootRank )
  {
    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      for( int i = 0; i < m_numDims; ++i )
      {
        if( periodic[i] != 1 )
        {
          continue;
        }

        real64 minimumPeriod = period[i];
        bool const affinePeriodicCell =
          m_prescribedFTable == 1 || m_stressControl[i] == 1;
        if( affinePeriodicCell )
        {
          real64 const endOfStepPeriod = period[i] * ( 1.0 + m_domainL[i] * dt );
          if( endOfStepPeriod <= 0.0 )
          {
            invalidPeriodicBody = bodyIndex;
            invalidPeriodicDirection = i;
            continue;
          }
          minimumPeriod = LvArray::math::min( minimumPeriod, endOfStepPeriod );
        }

        real64 const rotationSafeExtent = LvArray::math::max(
          globalPeriodicExtent[4 * bodyIndex + i],
          globalPeriodicExtent[4 * bodyIndex + 3] );
        if( 2.0 * rotationSafeExtent >= 0.98 * minimumPeriod )
        {
          invalidPeriodicBody = bodyIndex;
          invalidPeriodicDirection = i;
        }
      }
    }
  }
  MPI_Bcast( &invalidPeriodicBody, 1, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  MPI_Bcast( &invalidPeriodicDirection, 1, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidPeriodicBody >= 0,
                 "RigidBodyMPM color " << m_rigidBodyColors[invalidPeriodicBody]
                 << " reaches the half-period minimum-image limit during the current step in direction "
                 << invalidPeriodicDirection
                 << "; reduce periodic compaction or enlarge the periodic cell." );

  std::vector< real64 > globalAccumulator;
  gatherAndDeterministicallySumToRoot( localAccumulator, globalAccumulator );

  std::vector< real64 > solution( rigidBodySolutionStride * numberOfBodies, 0.0 );
  real64 rigidKineticEnergy = 0.0;
  real64 observedMaximumForce = 0.0;
  integer invalidMassBody = -1;

  if( rank == rigidBodyRootRank )
  {
    real64 const linearDampingScale =
      LvArray::math::max( 0.0, 1.0 - dt * m_rigidBodyLinearDamping );
    real64 const angularDampingScale =
      LvArray::math::max( 0.0, 1.0 - dt * m_rigidBodyAngularDamping );

    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      int const accumulatorOffset = rigidBodyAccumulatorStride * bodyIndex;
      int const solutionOffset = rigidBodySolutionStride * bodyIndex;
      real64 const mass = globalAccumulator[accumulatorOffset + rigidMass];
      if( mass <= 0.0 )
      {
        invalidMassBody = bodyIndex;
        continue;
      }

      real64 center[3] =
      {
        globalAccumulator[accumulatorOffset + rigidFirstMomentX] / mass,
        globalAccumulator[accumulatorOffset + rigidFirstMomentY] / mass,
        globalAccumulator[accumulatorOffset + rigidFirstMomentZ] / mass
      };
      real64 velocity[3] =
      {
        globalAccumulator[accumulatorOffset + rigidMomentumX] / mass,
        globalAccumulator[accumulatorOffset + rigidMomentumY] / mass,
        globalAccumulator[accumulatorOffset + rigidMomentumZ] / mass
      };

      real64 inertia = globalAccumulator[accumulatorOffset + rigidSecondMomentXY]
                       - mass * ( center[0] * center[0] + center[1] * center[1] );
      inertia = LvArray::math::max( inertia, 0.0 );
      real64 const angularMomentumCenter =
        globalAccumulator[accumulatorOffset + rigidAngularMomentumOrigin]
        - ( center[0] * globalAccumulator[accumulatorOffset + rigidMomentumY]
            - center[1] * globalAccumulator[accumulatorOffset + rigidMomentumX] );
      real64 const omega = inertia > 0.0 ? angularMomentumCenter / inertia : 0.0;

      real64 const mappedContactMass = globalAccumulator[accumulatorOffset + rigidMappedContactMass];
      real64 const contactMassRatio = mappedContactMass > 0.0 ? mass / mappedContactMass : 1.0;
      real64 contactImpulse[3] =
      {
        contactMassRatio * globalAccumulator[accumulatorOffset + rigidContactImpulseX],
        contactMassRatio * globalAccumulator[accumulatorOffset + rigidContactImpulseY],
        contactMassRatio * globalAccumulator[accumulatorOffset + rigidContactImpulseZ]
      };
      real64 const contactTorqueImpulse = contactMassRatio *
        ( globalAccumulator[accumulatorOffset + rigidContactTorqueImpulseOrigin]
          - ( center[0] * globalAccumulator[accumulatorOffset + rigidContactImpulseY]
              - center[1] * globalAccumulator[accumulatorOffset + rigidContactImpulseX] ) );

      real64 centerNew[3] = {};
      real64 velocityNew[3] = {};
      for( int i = 0; i < 3; ++i )
      {
        bool const affinePeriodicDirection =
          i < m_numDims && periodic[i] == 1 &&
          ( m_prescribedFTable == 1 || m_stressControl[i] == 1 );
        velocityNew[i] = linearDampingScale * ( velocity[i] + contactImpulse[i] / mass );
        if( affinePeriodicDirection )
        {
          // Periodic prescribedFTable kinematics store particle velocity as the
          // peculiar velocity. The affine box motion is applied to positions,
          // matching applySuperimposedVelocityGradient() in continuum mode.
          real64 const incrementalStretch = 1.0 + m_domainL[i] * dt;
          centerNew[i] = incrementalStretch * center[i] + dt * velocityNew[i];
        }
        else
        {
          centerNew[i] = center[i] + dt * velocityNew[i];
        }
      }

      real64 const omegaNew = inertia > 0.0
                              ? angularDampingScale * ( omega + contactTorqueImpulse / inertia )
                              : 0.0;
      real64 const force[3] =
      {
        contactImpulse[0] / dt,
        contactImpulse[1] / dt,
        contactImpulse[2] / dt
      };
      real64 const torque = contactTorqueImpulse / dt;
      real64 const forceNorm = LvArray::math::sqrt( force[0] * force[0]
                                                    + force[1] * force[1]
                                                    + force[2] * force[2] );
      observedMaximumForce = LvArray::math::max( observedMaximumForce, forceNorm );
      rigidKineticEnergy += 0.5 * mass * ( velocityNew[0] * velocityNew[0]
                                          + velocityNew[1] * velocityNew[1]
                                          + velocityNew[2] * velocityNew[2] )
                            + 0.5 * inertia * omegaNew * omegaNew;

      for( int i = 0; i < 3; ++i )
      {
        solution[solutionOffset + rigidCenterOldX + i] = center[i];
        solution[solutionOffset + rigidCenterNewX + i] = centerNew[i];
        solution[solutionOffset + rigidVelocityNewX + i] = velocityNew[i];
        solution[solutionOffset + rigidBodyForceX + i] = force[i];
      }
      solution[solutionOffset + rigidOmegaNewZ] = omegaNew;
      solution[solutionOffset + rigidBodyMass] = mass;
      solution[solutionOffset + rigidBodyInertiaZ] = inertia;
      solution[solutionOffset + rigidBodyMappedContactMass] = mappedContactMass;
      solution[solutionOffset + rigidBodyContactMassRatio] = contactMassRatio;
      solution[solutionOffset + rigidBodyTorqueZ] = torque;
      solution[solutionOffset + rigidBodyCosDeltaTheta] = std::cos( dt * omegaNew );
      solution[solutionOffset + rigidBodySinDeltaTheta] = std::sin( dt * omegaNew );
    }
  }

  MPI_Bcast( &invalidMassBody, 1, MPI_INT, rigidBodyRootRank, MPI_COMM_GEOS );
  GEOS_ERROR_IF( invalidMassBody >= 0,
                 "RigidBodyMPM color " << m_rigidBodyColors[invalidMassBody]
                 << " lost all surface-shell mass." );

  MPI_Bcast( solution.data(),
             static_cast< int >( solution.size() ),
             MPI_DOUBLE,
             rigidBodyRootRank,
             MPI_COMM_GEOS );
  MPI_Bcast( &rigidKineticEnergy, 1, MPI_DOUBLE, rigidBodyRootRank, MPI_COMM_GEOS );
  MPI_Bcast( &observedMaximumForce, 1, MPI_DOUBLE, rigidBodyRootRank, MPI_COMM_GEOS );

  m_rigidBodyKineticEnergy = rigidKineticEnergy;
  m_rigidBodyObservedMaxForce = observedMaximumForce;
  for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
  {
    int const solutionOffset = rigidBodySolutionStride * bodyIndex;
    for( int i = 0; i < 3; ++i )
    {
      m_rigidBodyUnwrappedCenters[3 * bodyIndex + i] =
        solution[solutionOffset + rigidCenterNewX + i];
    }
  }

  particleManager.forParticleSubRegions( [&]( ParticleSubRegion & subRegion )
  {
    // Move only the kinematic/surface fields updated by the rigid transform.
    // Avoid staging the full constitutive state for large particle counts.
    stdVector< string > const rigidUpdateFields =
    {
      fields::mpm::particleColor::key(),
      fields::mpm::particleMass::key(),
      fields::mpm::particleKineticEnergy::key(),
      fields::mpm::particleReferenceSurfaceNormal::key(),
      fields::mpm::particleReferenceSurfacePosition::key(),
      fields::mpm::particleReferenceRVectors::key(),
      fields::mpm::particleDeformationGradient::key(),
      fields::mpm::particleVelocityGradient::key(),
      ParticleSubRegion::viewKeyStruct::particleCenterString(),
      ParticleSubRegion::viewKeyStruct::particleVelocityString(),
      ParticleSubRegion::viewKeyStruct::particleSurfaceFlagString(),
      ParticleSubRegion::viewKeyStruct::particleSurfaceNormalString(),
      ParticleSubRegion::viewKeyStruct::particleSurfacePositionString(),
      ParticleSubRegion::viewKeyStruct::particleRVectorsString()
    };
    for( string const & fieldName : rigidUpdateFields )
    {
      subRegion.getWrapperBase( fieldName ).move( LvArray::MemorySpace::host, true );
    }

    arrayView1d< integer const > const particleColor = subRegion.getField< fields::mpm::particleColor >();
    arrayView2d< real64 > const particlePosition = subRegion.getParticleCenter();
    arrayView2d< real64 > const particleVelocity = subRegion.getParticleVelocity();
    arrayView2d< real64 > const particleSurfaceNormal = subRegion.getParticleSurfaceNormal();
    arrayView2d< real64 > const particleSurfacePosition = subRegion.getParticleSurfacePosition();
    arrayView3d< real64 > const particleRVectors = subRegion.getParticleRVectors();
    ParticleType const particleType = subRegion.getParticleType();
    arrayView2d< real64 > const particleReferenceSurfaceNormal =
      subRegion.getField< fields::mpm::particleReferenceSurfaceNormal >();
    arrayView2d< real64 > const particleReferenceSurfacePosition =
      subRegion.getField< fields::mpm::particleReferenceSurfacePosition >();
    arrayView3d< real64 > const particleReferenceRVectors =
      subRegion.getField< fields::mpm::particleReferenceRVectors >();
    arrayView3d< real64 > const particleDeformationGradient =
      subRegion.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 > const particleVelocityGradient =
      subRegion.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 const > const particleMass =
      subRegion.getField< fields::mpm::particleMass >();
    arrayView1d< real64 > const particleKineticEnergy =
      subRegion.getField< fields::mpm::particleKineticEnergy >();
    SortedArrayView< localIndex const > const activeParticleIndices = subRegion.activeParticleIndices();

    for( localIndex pp = 0; pp < activeParticleIndices.size(); ++pp )
    {
      localIndex const p = activeParticleIndices[pp];
      localIndex const bodyIndex = findRigidBodyIndex( colorToBodyIndex, particleColor[p] );
      GEOS_ERROR_IF( bodyIndex < 0, "RigidBodyMPM body registry is missing an active color." );
      int const solutionOffset = rigidBodySolutionStride * bodyIndex;
      real64 const cosTheta =
        solution[solutionOffset + rigidBodyCosDeltaTheta];
      real64 const sinTheta =
        solution[solutionOffset + rigidBodySinDeltaTheta];

      real64 const centerOld[3] =
      {
        solution[solutionOffset + rigidCenterOldX],
        solution[solutionOffset + rigidCenterOldY],
        solution[solutionOffset + rigidCenterOldZ]
      };
      real64 const centerNew[3] =
      {
        solution[solutionOffset + rigidCenterNewX],
        solution[solutionOffset + rigidCenterNewY],
        solution[solutionOffset + rigidCenterNewZ]
      };

      real64 positionUnwrapped[3] = {};
      for( int i = 0; i < 3; ++i )
      {
        positionUnwrapped[i] = minimumImageCoordinate( particlePosition[p][i],
                                                       centerOld[i],
                                                       period[i],
                                                       periodic[i] == 1 );
      }
      real64 const rxOld = positionUnwrapped[0] - centerOld[0];
      real64 const ryOld = positionUnwrapped[1] - centerOld[1];
      real64 const rzOld = positionUnwrapped[2] - centerOld[2];
      real64 const rxNew = cosTheta * rxOld - sinTheta * ryOld;
      real64 const ryNew = sinTheta * rxOld + cosTheta * ryOld;
      real64 const rzNew = rzOld;
      real64 positionNew[3] =
      {
        centerNew[0] + rxNew,
        centerNew[1] + ryNew,
        centerNew[2] + rzNew
      };

      for( int i = 0; i < m_numDims; ++i )
      {
        if( periodic[i] == 1 )
        {
          real64 targetMinimum = m_xGlobalMin[i];
          real64 targetPeriod = m_domainExtent[i];
          bool const affinePeriodicCell =
            m_prescribedFTable == 1 || m_stressControl[i] == 1;
          if( affinePeriodicCell )
          {
            real64 const incrementalStretch = 1.0 + m_domainL[i] * dt;
            targetMinimum *= incrementalStretch;
            targetPeriod *= incrementalStretch;
          }
          positionNew[i] = wrapCoordinate( positionNew[i], targetMinimum, targetPeriod );
        }
      }

      real64 const omegaNew = solution[solutionOffset + rigidOmegaNewZ];
      particleVelocity[p][0] = solution[solutionOffset + rigidVelocityNewX] - omegaNew * ryNew;
      particleVelocity[p][1] = solution[solutionOffset + rigidVelocityNewY] + omegaNew * rxNew;
      particleVelocity[p][2] = solution[solutionOffset + rigidVelocityNewZ];
      for( int i = 0; i < 3; ++i )
      {
        particlePosition[p][i] = positionNew[i];
      }

      auto rotateVector2D = [=]( arraySlice1d< real64 > const vector )
      {
        real64 const x = vector[0];
        real64 const y = vector[1];
        vector[0] = cosTheta * x - sinTheta * y;
        vector[1] = sinTheta * x + cosTheta * y;
      };
      rotateVector2D( particleSurfaceNormal[p] );
      rotateVector2D( particleSurfacePosition[p] );
      if( particleType != ParticleType::SinglePoint )
      {
        for( int a = 0; a < 3; ++a )
        {
          rotateVector2D( particleRVectors[p][a] );
        }
      }

      for( int i = 0; i < 3; ++i )
      {
        particleReferenceSurfaceNormal[p][i] = particleSurfaceNormal[p][i];
        particleReferenceSurfacePosition[p][i] = particleSurfacePosition[p][i];
      }
      for( int a = 0; a < 3; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          particleReferenceRVectors[p][a][i] = particleRVectors[p][a][i];
        }
      }
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          particleDeformationGradient[p][i][j] = i == j ? 1.0 : 0.0;
          particleVelocityGradient[p][i][j] = 0.0;
        }
      }

      real64 const speedSquared = particleVelocity[p][0] * particleVelocity[p][0]
                                  + particleVelocity[p][1] * particleVelocity[p][1]
                                  + particleVelocity[p][2] * particleVelocity[p][2];
      particleKineticEnergy[p] = 0.5 * particleMass[p] * speedSquared;
    }

#ifdef GEOS_USE_DEVICE
    // The rigid transform is intentionally host-side for now. Move the updated
    // particle state back to the device before the next explicit step.
    for( string const & fieldName : rigidUpdateFields )
    {
      subRegion.getWrapperBase( fieldName ).move( parallelDeviceMemorySpace, true );
    }
#endif
  } );

  if( m_rigidBodyHistory == 1 && rank == rigidBodyRootRank
      && ( m_rigidBodyHistoryWriteInterval <= 0.0
           || time_n + dt >= m_nextRigidBodyHistoryWriteTime ) )
  {
    std::ofstream history( "rigidBodyHistory.csv", std::ios::app );
    for( localIndex bodyIndex = 0; bodyIndex < numberOfBodies; ++bodyIndex )
    {
      int const offset = rigidBodySolutionStride * bodyIndex;
      real64 const bodyKineticEnergy =
        0.5 * solution[offset + rigidBodyMass]
        * ( solution[offset + rigidVelocityNewX] * solution[offset + rigidVelocityNewX]
            + solution[offset + rigidVelocityNewY] * solution[offset + rigidVelocityNewY]
            + solution[offset + rigidVelocityNewZ] * solution[offset + rigidVelocityNewZ] )
        + 0.5 * solution[offset + rigidBodyInertiaZ]
          * solution[offset + rigidOmegaNewZ] * solution[offset + rigidOmegaNewZ];

      real64 const maxContactPenetrationRatio =
        m_rigidBodyContactLengthScale > 0.0
        ? m_rigidBodyObservedMaxPenetration / m_rigidBodyContactLengthScale
        : 0.0;

      history << time_n + dt << ',' << m_rigidBodyColors[bodyIndex] << ','
              << solution[offset + rigidBodyMass] << ','
              << solution[offset + rigidBodyMappedContactMass] << ','
              << solution[offset + rigidBodyContactMassRatio] << ','
              << m_rigidBodyObservedJammingStress << ','
              << m_rigidBodyContactLengthScale << ','
              << m_rigidBodyObservedMaxPenetration << ','
              << maxContactPenetrationRatio << ','
              << solution[offset + rigidCenterNewX] << ','
              << solution[offset + rigidCenterNewY] << ','
              << solution[offset + rigidCenterNewZ] << ','
              << solution[offset + rigidVelocityNewX] << ','
              << solution[offset + rigidVelocityNewY] << ','
              << solution[offset + rigidVelocityNewZ] << ','
              << solution[offset + rigidOmegaNewZ] << ','
              << solution[offset + rigidBodyForceX] << ','
              << solution[offset + rigidBodyForceY] << ','
              << solution[offset + rigidBodyForceZ] << ','
              << solution[offset + rigidBodyTorqueZ] << ','
              << bodyKineticEnergy << '\n';
    }
    m_nextRigidBodyHistoryWriteTime =
      time_n + dt + LvArray::math::max( 0.0, m_rigidBodyHistoryWriteInterval );
  }
}

} // namespace geos
