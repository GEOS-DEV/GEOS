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

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const );


  PartitionBase * GetPartition() {return m_partition;}
private:

  PartitionBase * m_partition;

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_ */
