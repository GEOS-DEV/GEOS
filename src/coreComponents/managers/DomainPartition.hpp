/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file DomainPartition.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_

#include "dataRepository/ManagedGroup.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/ConstitutiveManager.hpp"
namespace geosx
{

class SiloFile;
namespace dataRepository
{
namespace keys
{
string const partitionManager("partitionManager");
}
}

class ObjectManagerBase;
class PartitionBase;

class DomainPartition : public dataRepository::ManagedGroup
{
public:
  DomainPartition( std::string const & name,
                   ManagedGroup * const parent );

  ~DomainPartition() override;

  DomainPartition() = delete;
  DomainPartition( DomainPartition const &) = delete;
  DomainPartition( DomainPartition &&) = delete;
  DomainPartition& operator=( DomainPartition const & ) = delete;
  DomainPartition& operator=( DomainPartition && ) = delete;

  virtual void RegisterDataOnMeshRecursive( ManagedGroup * const MeshBodies ) override final;


  void InitializationOrder( string_array & order ) override final;

  void SetMaps();
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

  struct viewKeysStruct
  {
    dataRepository::ViewKey neighbors = { "Neighbors" };
  } viewKeys;

  struct groupKeysStruct
  {
    static constexpr auto meshBodiesString = "MeshBodies";
    static constexpr auto constitutiveManagerString = "ConstitutiveManager";

    dataRepository::GroupKey meshBodies           = { meshBodiesString };
    dataRepository::GroupKey constitutiveManager  = { constitutiveManagerString };
    dataRepository::GroupKey communicationManager    = { "communicationManager" };
  } groupKeys;


  constitutive::ConstitutiveManager const * getConstitutiveManager() const
  { return this->GetGroup<constitutive::ConstitutiveManager>(groupKeys.constitutiveManager); }

  constitutive::ConstitutiveManager * getConstitutiveManager()
  { return this->GetGroup<constitutive::ConstitutiveManager>(groupKeys.constitutiveManager); }


  ManagedGroup const * getMeshBodies() const
  { return this->GetGroup(groupKeys.meshBodies); }
  ManagedGroup * getMeshBodies()
  { return this->GetGroup(groupKeys.meshBodies); }

  MeshBody const * getMeshBody( string const & meshName ) const
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(meshName); }
  MeshBody * getMeshBody( string const & meshName )
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(meshName); }

  MeshBody const * getMeshBody( integer const index ) const
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(index); }
  MeshBody * getMeshBody( integer const index )
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(index); }


private:


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_ */
