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

#include "dataRepository/Group.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

namespace geosx
{

class SiloFile;
namespace dataRepository
{
namespace keys
{
string const partitionManager( "partitionManager" );
}
}

class ObjectManagerBase;
class PartitionBase;

class DomainPartition : public dataRepository::Group
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

  struct groupKeysStruct
  {
    static constexpr auto meshBodiesString = "MeshBodies";
    static constexpr auto constitutiveManagerString = "Constitutive";

    dataRepository::GroupKey meshBodies           = { meshBodiesString };
    dataRepository::GroupKey constitutiveManager  = { constitutiveManagerString };
    dataRepository::GroupKey communicationManager    = { "communicationManager" };
  } groupKeys;


  constitutive::ConstitutiveManager const * getConstitutiveManager() const
  { return this->GetGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }

  constitutive::ConstitutiveManager * getConstitutiveManager()
  { return this->GetGroup< constitutive::ConstitutiveManager >( groupKeys.constitutiveManager ); }


  Group const * getMeshBodies() const
  { return this->GetGroup( groupKeys.meshBodies ); }

  Group * getMeshBodies()
  { return this->GetGroup( groupKeys.meshBodies ); }

  MeshBody const * getMeshBody( string const & meshName ) const
  { return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( meshName ); }

  MeshBody * getMeshBody( string const & meshName )
  { return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( meshName ); }

  MeshBody const * getMeshBody( localIndex const index ) const
  { return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( index ); }

  MeshBody * getMeshBody( localIndex const index )
  { return this->GetGroup( groupKeys.meshBodies )->GetGroup< MeshBody >( index ); }

  std::set< int > & getMetisNeighborList()
  { return m_metisNeighborList; }

  std::set< int > const & getMetisNeighborList() const
  { return m_metisNeighborList; }

  std::vector< NeighborCommunicator > & getNeighbors()
  { return m_neighbors; }

  std::vector< NeighborCommunicator > const & getNeighbors() const
  { return m_neighbors; };

private:

  std::set< int > m_metisNeighborList;
  std::vector< NeighborCommunicator > m_neighbors;

};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_DOMAINPARTITION_HPP_ */
