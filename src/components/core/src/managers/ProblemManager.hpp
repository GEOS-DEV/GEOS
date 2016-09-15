/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#ifndef USE_PYTHON
#define USE_PYTHON 0
#endif

#if USE_PYTHON==1
// Note: the python header must be included first to avoid conflicting definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif
#include "pugixml/src/pugixml.hpp"

#include "dataRepository/ManagedGroup.hpp"
#include "PhysicsSolvers/PhysicsSolverManager.hpp"
#include "InputDocumentation.hpp"
#include "schema/SchemaUtilities.hpp"

namespace geosx
{

class DomainPartition;

class ProblemManager : public dataRepository::ManagedGroup
{
public:
  explicit ProblemManager( const std::string& name,
                           ManagedGroup * const parent );
  ~ProblemManager();

  static std::string CatalogName() { return "ProblemManager"; }


  virtual void Registration( dataRepository::ManagedGroup * const );

  void ParseCommandLineInput( int const& argc, char* const argv[]);

  void InitializePythonInterpreter();

  void ClosePythonInterpreter();

  void ParseInputFile();

  void InitializeObjects();

  void RunSimulation();

  void ApplySchedulerEvent();

  DomainPartition & getDomainPartition();
  DomainPartition const & getDomainPartition() const;

  pugi::xml_document xmlDocument;
  pugi::xml_parse_result xmlResult;
  pugi::xml_node xmlProblemNode;

private:
  PhysicsSolverManager & m_physicsSolverManager;
  cxx_utilities::InputDocumentation m_inputDocumentationHead;
};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
