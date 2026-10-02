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
 * @file VTKRefinementCommunication.hpp
 */

#ifndef GEOS_VTK_REFINEMENT_COMMUNICATION_HPP
#define GEOS_VTK_REFINEMENT_COMMUNICATION_HPP

#include "VTKRefinementTopology.hpp"
#include "VTKRefinementAssociations.hpp"
#include "common/MpiWrapper.hpp"

#include <functional>
#include <map>
#include <string>

namespace geos::vtk::refinement
{

using Bytes = std::vector< unsigned char >;

struct PointCreation
{
  EntityKey key;
  Participants participants;
  Coordinates position;
  Bytes fields;            ///< Canonically encoded point data, interpreted by GEOS'
                           ///< transfer policies.
  double supportScale = 1; ///< Support extent; zero is allowed for isolated existing vertices (roundoff tolerance only).
};

struct PointRecord
{
  vtkIdType globalId;
  Coordinates position;
  Bytes fields;
};

struct IdRange
{
  vtkIdType first; ///< Only meaningful for a nonzero local count.
  std::uint64_t total;
};

/** Identity is generation + this full key, independently of allocated IDs. */
struct ChildCellKey
{
  std::uint64_t meshNamespace;
  vtkIdType parentId;
  std::uint64_t templateCode;
  std::uint64_t ordinal;
  bool operator<( ChildCellKey const & other ) const;
};

struct CellCreation
{
  ChildCellKey key;
  Participants participants;
  Bytes fields;
};

struct CellRecord
{
  vtkIdType globalId;
  Bytes fields;
};

struct CoarseSurface
{
  std::vector< Connectivity > cornerBuckets; ///< Clean main GIDs, in surface-corner order.
};

struct SupportLookup
{
  EntityKey key;
  int faceOwner; ///< A rank holding this actual incident volume face.
};

struct CommunicationStatistics
{
  std::uint64_t directoryExchanges{};
  std::uint64_t neighborExchanges{};
  std::uint64_t payloadChunksSent{};
  std::uint64_t payloadBytesSent{};
  std::uint64_t countBytesSent{};
};

/** GEOS mesh protocol, separate from LvArray storage and pure cell templates.
 * The communicator is duplicated to isolate tags. Discovery runs on coarse
 * entities; subsequent point resolution only exchanges with cached neighbors.
 * Every public protocol call is collective on the supplied communicator.
 */
class Communication
{
public:
  explicit Communication( MPI_Comm comm, std::uint64_t chunkBytes = UINT64_C( 1 ) << 30 );
  ~Communication();
  Communication( Communication const & ) = delete;
  Communication & operator=( Communication const & ) = delete;

  void checked( std::string const & phase, std::function< void() > const & work ) const;
  Sharing discoverSharing( std::vector< EntityKey > const & entities );
  /** One coarse full-face pass rejects global nonmanifold and crossed-cycle
   * incidence, including a locally internal face also used by another rank.
   * Three-corner quad probes also reject quad/triangle nonmatching interfaces.
   * Neighboring volume cells need not use consistently oriented face cycles:
   * refinement keys canonicalize face orientation independently of the input.
   * This validation cost is separate from boundary-only sharing discovery.
   */
  void validateVolumeFaces( std::vector< MainFace > const & faces, std::uint64_t mainNamespace = 0 );
  IdRange allocateRange( std::uint64_t localCount, vtkIdType base ) const;
  /** Validate/synchronize existing shared vertices while preserving their IDs. */
  std::map< EntityKey, PointRecord > reconcileExistingPoints( std::vector< PointCreation > const & points );
  std::map< EntityKey, PointRecord > resolvePoints( std::uint64_t generation, std::vector< PointCreation > const & points,
                                                    vtkIdType localExistingMaximum );
  /** Allocate/synchronize replicated marker or auxiliary surface children.
   * Optional allocatedRange receives the range total without another collective allocation.
   */
  std::map< ChildCellKey, CellRecord > resolveCells( std::uint64_t generation, std::vector< CellCreation > const & cells,
                                                    vtkIdType base, IdRange * allocatedRange = nullptr );
  std::map< ChildCellKey, CellRecord > reconcileExistingCells( std::vector< CellCreation > const & cells );
  /** Discover actual local/remote coarse side faces once, using vertex-ID routing.
   * Query holders are not treated as owners of main points or faces. The resulting
   * contact neighbors are cached for subsequent support-ID lookup.
   */
  std::vector< std::vector< SurfaceSide > > discoverSurfaceSides( std::vector< MainFace > const & localFaces,
                                                                  std::vector< CoarseSurface > const & surfaces,
                                                                  std::uint64_t mainNamespace = 0 );
  /** Resolve surface recipes from an actual face owner, without a fine directory. */
  std::vector< vtkIdType > resolveSupportIds( std::uint64_t generation, std::vector< SupportLookup > const & requests,
                                              std::unordered_map< EntityKey, vtkIdType, EntityKeyHash > const & localIds );

  Participants const & neighbors() const { return m_neighbors; }
  CommunicationStatistics const & statistics() const { return m_statistics; }
  int rank() const { return m_rank; }
  int size() const { return m_size; }

private:
  using Mail = std::map< int, Bytes >;
  Mail exchangeDirectory( Mail const & outgoing );
  Mail exchangeNeighbors( Mail const & outgoing );
  Mail exchangePayloads( Mail const & outgoing, std::map< int, std::uint64_t > const & incomingLengths );
  void initializeNeighbors( Participants neighbors );
  void includeContactNeighbors( Participants neighbors );
  std::map< EntityKey, PointRecord > resolvePointRecords( std::uint64_t generation, std::vector< PointCreation > const & points,
                                                          vtkIdType localExistingMaximum, bool existing );
  std::map< ChildCellKey, CellRecord > resolveCellRecords( std::uint64_t generation, std::vector< CellCreation > const & cells,
                                                           vtkIdType base, bool existing, IdRange * allocatedRange = nullptr );

  MPI_Comm m_comm;
  MPI_Comm m_neighborComm = MPI_COMM_NULL;
  int m_rank;
  int m_size;
  std::uint64_t m_chunkBytes;
  Participants m_neighbors, m_sources, m_destinations;
  bool m_discovered = false;
  CommunicationStatistics m_statistics;
};

} // namespace geos::vtk::refinement
#endif
