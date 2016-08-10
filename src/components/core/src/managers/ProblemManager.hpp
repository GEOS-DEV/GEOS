/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#include "dataRepository/WrapperCollection.hpp"

namespace geosx
{

class DomainPartition;

class ProblemManager : public dataRepository::WrapperCollection
{
public:
  explicit ProblemManager( const std::string& name,
                           WrapperCollection * const parent );
  ~ProblemManager();

  static std::string CatalogName() { return "ProblemManager"; }


  virtual void Registration( dataRepository::WrapperCollection * const );

  void ParseCommandLineInput( int const& argc, char* const argv[]);

  void ParseInputFile();

  void ApplySchedulerEvent();

  std::unordered_map<string,string>& simulationParameterMap()
  {
    return getReference< std::unordered_map<string,string> >(keys::simulationParameterMap);
  }

  std::unordered_map<string,string> m_simulationParameterMap;

  DomainPartition & getDomainPartition();
  DomainPartition const & getDomainPartition() const;

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
