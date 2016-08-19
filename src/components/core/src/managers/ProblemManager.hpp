/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#include "../dataRepository/SynchronizedGroup.hpp"

namespace geosx
{

class DomainPartition;

class ProblemManager : public dataRepository::SynchronizedGroup
{
public:
  explicit ProblemManager( const std::string& name,
                           SynchronizedGroup * const parent );
  ~ProblemManager();

  static std::string CatalogName() { return "ProblemManager"; }


  virtual void Registration( dataRepository::SynchronizedGroup * const );

  void ParseCommandLineInput( int const& argc, char* const argv[]);

  void ParseInputFile();

  void ApplySchedulerEvent();

  DomainPartition & getDomainPartition();
  DomainPartition const & getDomainPartition() const;

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
