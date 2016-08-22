/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

<<<<<<< HEAD
// Note: the python header must be included first to avoid conflicting definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include "dataRepository/WrapperCollection.hpp"
#include <pugixml.hpp>
=======
#include "../dataRepository/SynchronizedGroup.hpp"
>>>>>>> feature/settgast/addLegacyCode

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
