/*
 * DomainPartition.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
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

class PartitionBase;

class DomainPartition : public dataRepository::ManagedGroup
{
public:
  DomainPartition( std::string const & name,
                   ManagedGroup * const parent );

  ~DomainPartition();

  DomainPartition() = delete;
  DomainPartition( DomainPartition const &) = delete;
  DomainPartition( DomainPartition &&) = delete;
  DomainPartition& operator=( DomainPartition const & ) = delete;
  DomainPartition& operator=( DomainPartition && ) = delete;

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const ) override;
  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  void InitializationOrder( string_array & order ) override final;

  void SetMaps();
  void GenerateSets();


  // THIS STUFF NEEDS TO GO SOMEWHERE ELSE
  void WriteSilo( SiloFile& siloFile,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart );

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
    dataRepository::ViewKey partitionManager    = { "partitionManager" };
  }viewKeys;

  struct groupKeysStruct
  {
    static constexpr auto meshBodiesString = "MeshBodies";
    static constexpr auto constitutiveManagerString = "ConstitutiveManager";

    dataRepository::GroupKey meshBodies           = { meshBodiesString };
    dataRepository::GroupKey constitutiveManager  = { constitutiveManagerString };
  }groupKeys;


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

  MeshBody const * getMeshBody( int32 const index ) const
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(index); }
  MeshBody * getMeshBody( int32 const index )
  { return this->GetGroup(groupKeys.meshBodies)->GetGroup<MeshBody>(index); }


private:


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_ */
