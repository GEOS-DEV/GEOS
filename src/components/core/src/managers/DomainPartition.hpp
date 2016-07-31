/*
 * DomainPartition.hpp
 *
 *  Created on: Jul 28, 2016
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_

#include "dataRepository/WrapperCollection.hpp"

namespace geosx
{


class DomainPartition : public dataRepository::WrapperCollection
{
public:
  DomainPartition( std::string const & name,
                   WrapperCollection * const parent );

  ~DomainPartition();

  DomainPartition() = delete;
  DomainPartition( DomainPartition const &) = delete;
  DomainPartition( DomainPartition &&) = delete;
  DomainPartition& operator=( DomainPartition const & ) = delete;
  DomainPartition& operator=( DomainPartition && ) = delete;

  virtual void Registration( dataRepository::WrapperCollection * const );

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_MANAGERS_DOMAINPARTITION_HPP_ */
