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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_

#include "dataRepository/Group.hpp"
#include "PartitionBase.hpp"
#include "mesh/DomainPartition.hpp"


#include <map>

constexpr int nsdof = 3;
namespace geos
{

// inline bool isEqual( const real64& val1, const real64& val2, const real64& tolfac=0.0 )
// {
//   realT tol = 0.0;
//   if( tolfac > 1.0e-15 )
//     tol = fabs(tolfac) * (fabs(val1)+fabs(val2))*0.5;
//   return val1<=(val2+tol) && val1>=(val2-tol);
// }

// CC: Taken from old geos
// Planar Sorter
// Sorts pairs of local and global indexes by the positions of their corresponding node points in a plane.
class PlanarSorter {

public:
	PlanarSorter(const arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >& refPos, int dim) :
               dimension(dim), 
               refPositions(refPos) {};

	// sort operator for pairs containing local indexes (sort based on 1st element in pair)
	bool operator()(const std::pair<localIndex, localIndex>& lhs,
	                const std::pair<localIndex, localIndex>& rhs) 
  {
		bool rv = false;
		int a = 0;
		int b = 2;
		if (dimension == 0)
			a = 1;
		if (dimension == 2)
			b = 1;

		const arraySlice1d<real64 const>& lhsVect = refPositions[lhs.first];
		const arraySlice1d<real64 const>& rhsVect = refPositions[rhs.first];

		if (lhsVect[a] < rhsVect[a]) {
			rv = true;
		} else if (isEqual(lhsVect[a], rhsVect[a])
				&& (lhsVect[b] < rhsVect[b])) {
			rv = true;
		};

		return rv;
	};

	// sort operator for local indexes
	bool operator()(const localIndex& lhs, 
                  const localIndex& rhs)
  {
		bool rv = false;
		int a = 0;
		int b = 2;
		if (dimension == 0)
			a = 1;
		if (dimension == 2)
			b = 1;

		const arraySlice1d<real64 const>& lhsVect = refPositions[lhs];
		const arraySlice1d<real64 const>& rhsVect = refPositions[rhs];

		if (lhsVect[a] < rhsVect[a]) {
			rv = true;
		} else if (isEqual(lhsVect[a], rhsVect[a])
				&& (lhsVect[b] < rhsVect[b])) {
			rv = true;
		};

		return rv;
	};

private:
	int dimension;
	const arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >& refPositions;
};

/**
 * @brief Concrete (cartesian?) partitioning.
 */
class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition(string const & name,
                   Group * const parent );

  ~SpatialPartition() override;

  bool isCoordInPartition( const real64 & coord, const int dir ) const override;

  bool isCoordInPartitionBoundingBox( const R1Tensor & elemCenter,
                                      const real64 & boundaryRadius ) const;

  void updateSizes( arrayView1d< real64 > const domainL,
                    real64 const dt );

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] ) override;

  // real64 * getLocalMin()
  arrayView1d< real64 > getLocalMin()
  {
    return m_min;
  }

  // real64 * getLocalMax()
  arrayView1d< real64 > getLocalMax()
  {
    return m_max;
  }

  // real64 * getGlobalMin()
  arrayView1d< real64 > getGlobalMin()
  {
    return m_gridMin;
  }

  // real64 * getGlobalMax()
  arrayView1d< real64 > getGlobalMax()
  {
    return m_gridMax;
  }

  arrayView1d< int const > const getPeriodic() const {
    return m_Periodic;
  }

  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  int getColor() override;

  void repartitionMasterParticles( ParticleSubRegion & subRegion,
                                   MPI_iCommData & commData );

  void getGhostParticlesFromNeighboringPartitions( DomainPartition & domain,
                                                   MPI_iCommData & commData,
                                                   const real64 & boundaryRadius );

  /**
   * @brief Send coordinates to neighbors as part of repartition.
   * @param[in] particleCoordinatesSendingToNeighbors Single list of coordinates sent to all neighbors
   * @param[in] commData Solver's MPI communicator
   * @param[in] particleCoordinatesReceivedFromNeighbors List of lists of coordinates received from each neighbor
   */
  void sendCoordinateListToNeighbors( arrayView1d< R1Tensor > const & particleCoordinatesSendingToNeighbors,
                                      MPI_iCommData & commData,
                                      std::vector< array1d< R1Tensor > > & particleCoordinatesReceivedFromNeighbors
                                      );

  template< typename indexType >
  void sendListOfIndicesToNeighbors( std::vector< array1d< indexType > > & listSendingToEachNeighbor,
                                     MPI_iCommData & commData,
                                     std::vector< array1d< indexType > > & listReceivedFromEachNeighbor );

  void sendParticlesToNeighbor( ParticleSubRegionBase & subRegion,
                                std::vector< int > const & newParticleStartingIndices,
                                std::vector< int > const & numberOfIncomingParticles,
                                MPI_iCommData & commData,
                                std::vector< array1d< localIndex > > const & particleLocalIndicesToSendToEachNeighbor );

  void setPeriodicDomainBoundaryObjects(  MeshBody & grid,
                                          NodeManager & nodeManager,
                                          EdgeManager & edgeManager,
                                          FaceManager & faceManager );

  struct viewKeyStruct
  {
    static constexpr char const * periodicString() { return "periodic"; }
    static constexpr char const * minString() { return "min"; }
    static constexpr char const * maxString() { return "max"; }
    static constexpr char const * partitionLocationsString() { return "partitionLocations"; }
    static constexpr char const * blockSizeString() { return "blockSize"; }
    static constexpr char const * gridSizeString() { return "gridSize"; }
    static constexpr char const * gridMinString() { return "gridMin"; }
    static constexpr char const * gridMaxString() { return "gridMax"; }
    static constexpr char const * contactGhostMinString() { return "contactGhostMin"; }
    static constexpr char const * contactGhostMaxString() { return "contactGhostMax"; }
  } partitionViewKeys;

  static string catalogName() { return "SpatialPartition"; }

  void postProcessInput() override; 

 /// number of partitions
  array1d< int > m_Partitions;
  
  /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;

  /// ijk partition indexes
  array1d< int > m_coords;

private:

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );

  /// Minimum extent of partition dimensions (excluding ghost objects)
  // real64 m_min[3];
  array1d< real64 > m_min;

  /// Maximum extent of partition dimensions (excluding ghost objects)
  // real64 m_max[3];
  array1d< real64 > m_max;

  /// Locations of partition boundaries
  array1d< array1d< real64 > > m_partitionLocations;

  /// Length of partition dimensions (excluding ghost objects).
  // real64 m_blockSize[3];
  array1d< real64 > m_blockSize;

  /// Total length of problem dimensions (excluding ghost objects).
  // real64 m_gridSize[3];
  array1d< real64 > m_gridSize;
  
  /// Minimum extent of problem dimensions (excluding ghost objects).
  // real64 m_gridMin[3];
  array1d< real64 > m_gridMin;

  /// Maximum extent of problem dimensions (excluding ghost objects).
  // real64 m_gridMax[3];
  array1d< real64 > m_gridMax;

  /**
   * @brief Ghost position (min).
   */
  // real64 m_contactGhostMin[3];
  array1d< real64 > m_contactGhostMin;

  /**
   * @brief Ghost position (max).
   */
  // real64 m_contactGhostMax[3];
  array1d< real64 > m_contactGhostMax;
};

}
#endif /* GEOS_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
