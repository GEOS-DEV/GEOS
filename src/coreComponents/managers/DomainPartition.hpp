/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DomainPartition.hpp
 */

#ifndef GEOSX_MANAGERS_DOMAINPARTITION_HPP_
#define GEOSX_MANAGERS_DOMAINPARTITION_HPP_

#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/Group.hpp"
#include "managers/ProblemManager.hpp"
#include "mesh/MeshBody.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/SpatialPartition.hpp"

namespace geosx
{

class SiloFile;

namespace dataRepository {
namespace keys {

string const partitionManager( "partitionManager" );

} /* end of namespace keys */
} /* end of namespace dataRepository */

class PartitionBase;

class DomainPartition : private dataRepository::GroupDownCastHelper<DomainPartition>
{
public:
  DomainPartition( std::string const & name,
                   Group * const parent );

  ~DomainPartition() override;

  DomainPartition() = delete;
  DomainPartition( DomainPartition const & ) = delete;
  DomainPartition( DomainPartition && ) = delete;
  DomainPartition & operator=( DomainPartition const & ) = delete;
  DomainPartition & operator=( DomainPartition && ) = delete;

  virtual void RegisterDataOnMeshRecursive( Group * const MeshBodies ) override final;

  void InitializationOrder( string_array & order ) override final;

  void GenerateSets();

  /**
   * @name MPI functionality
   */
  ///@{

  void SetupCommunications( bool use_nonblocking );

  void AddNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );
  ///@}

  void ReadSilo( const SiloFile & siloFile,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart );

  void WriteFiniteElementMesh( SiloFile & siloFile,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart );

  void ReadFiniteElementMesh( const SiloFile & siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart );

  constitutive::ConstitutiveManager const * getConstitutiveManager() const;
  constitutive::ConstitutiveManager * getConstitutiveManager();

  Group const * getMeshBodies() const;
  Group * getMeshBodies();

  MeshBody const * getMeshBody(string const & meshName ) const;
  MeshBody * getMeshBody(string const & meshName );

  MeshBody const * getMeshBody(localIndex const index ) const;
  MeshBody * getMeshBody(localIndex const index );

  ProblemManager const * GetProblemManager() const;
  ProblemManager * GetProblemManager();

  PartitionBase const & GetPartitionBase() const;
  PartitionBase & GetPartitionBase();

  CellBlockManager const * GetCellManager() const;
  CellBlockManager * GetCellManager();

  SpatialPartition const & GetSpatialPartition() const;
  SpatialPartition & GetSpatialPartition();

  void SetMetisNeighborList(std::set<int> const & others);

  std::vector< NeighborCommunicator > & getNeighbors()
  { return m_neighbors; }

  std::vector< NeighborCommunicator > const & getNeighbors() const
  { return m_neighbors; };

private:

  struct viewKeysStruct
  {
    dataRepository::ViewKey neighbors = { "Neighbors" };
  } viewKeys;

  struct groupKeysStruct
  {
    static constexpr auto meshBodiesString = "MeshBodies";
    static constexpr auto constitutiveManagerString = "Constitutive";

    dataRepository::GroupKey meshBodies           = { meshBodiesString };
    dataRepository::GroupKey constitutiveManager  = { constitutiveManagerString };
    dataRepository::GroupKey communicationManager    = { "communicationManager" };
  } groupKeys;

  std::set<int> m_metisNeighborList;
  std::vector< NeighborCommunicator > m_neighbors;

  };

} /* end of namespace geosx */

#endif /* GEOSX_MANAGERS_DOMAINPARTITION_HPP_ */
