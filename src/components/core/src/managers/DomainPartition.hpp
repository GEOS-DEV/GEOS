/*
 * DomainPartition.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_

#include "dataRepository/ManagedGroup.hpp"

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

  //  PartitionBase * GetPartition() {return m_partitio}
private:


};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_ */
