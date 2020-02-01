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

string const partitionManager("partitionManager");

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
  DomainPartition( DomainPartition const &) = delete;
  DomainPartition( DomainPartition &&) = delete;
  DomainPartition& operator=( DomainPartition const & ) = delete;
  DomainPartition& operator=( DomainPartition && ) = delete;

  virtual void RegisterDataOnMeshRecursive( Group * const MeshBodies ) override final;

  void InitializationOrder( string_array & order ) override final;

  void GenerateSets();

  /**
   * @name MPI functionality
   */
  ///@{

//  void FindMatchedPartitionBoundaryObjects( ObjectManagerBase * const group,
//                                            array1d< array1d<localIndex> > & matchedPartitionBoundaryObjects );

//  static std::set<int> & getFreeCommIDs();
//  static int reserveCommID();
//  static void releaseCommID( int & ID );

  void SetupCommunications();

  void AddNeighbors(const unsigned int idim,
                    MPI_Comm& cartcomm,
                    int* ncoords);
  ///@}

  void ReadSilo( const SiloFile& siloFile,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart );

  void WriteFiniteElementMesh( SiloFile& siloFile,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart );

  void ReadFiniteElementMesh( const SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart );

  constitutive::ConstitutiveManager const * GetConstitutiveManager() const ;
  constitutive::ConstitutiveManager * GetConstitutiveManager() ;

  Group const * getMeshBodies() const ;
  Group * getMeshBodies() ;

  MeshBody const * getMeshBody(string const & meshName ) const ;
  MeshBody * getMeshBody(string const & meshName ) ;

  MeshBody const * getMeshBody(localIndex const index ) const ;
  MeshBody * getMeshBody(localIndex const index ) ;

  ProblemManager const * GetProblemManager() const ;
  ProblemManager * GetProblemManager() ;

  array1d<NeighborCommunicator> const & GetNeighborCommunicators() const;
  array1d<NeighborCommunicator> & GetNeighborCommunicators() ;

  PartitionBase const & GetPartitionBase() const ;
  PartitionBase & GetPartitionBase() ;

  CellBlockManager const * GetCellManager() const ;
  CellBlockManager * GetCellManager() ;

  SpatialPartition const & GetSpatialPartition() const ;
  SpatialPartition & GetSpatialPartition();

  void SetMetisNeighborList(std::set<int> const & others) ;

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

  };

} /* end of namespace geosx */

#endif /* GEOSX_MANAGERS_DOMAINPARTITION_HPP_ */
