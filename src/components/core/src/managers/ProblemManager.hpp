/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#include "dataRepository/WrapperCollection.hpp"
#include <pugixml.hpp>
// #include <Python.h>
// #include <numpy/arrayobject.h>

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

  void InitializePythonInterpreter();

  void ClosePythonInterpreter();

  void ParseInputFile();

  void ApplySchedulerEvent();

  DomainPartition & getDomainPartition();
  DomainPartition const & getDomainPartition() const;

  pugi::xml_document xmlDocument;
  pugi::xml_parse_result xmlResult;
  pugi::xml_node xmlProblemNode;

};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
