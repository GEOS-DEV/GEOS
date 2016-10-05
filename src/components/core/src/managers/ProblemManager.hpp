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
#include "optionparser/src/optionparser.h"

#include "ObjectManagerBase.hpp"
#include "PhysicsSolvers/PhysicsSolverManager.hpp"
#include "EventManager.hpp"
#include "schema/SchemaUtilities.hpp"
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"

namespace geosx
{

struct Arg : public option::Arg
{
  static option::ArgStatus Unknown(const option::Option& option, bool /*error*/)
  {
    std::cout << "Unknown option: " << option.name << std::endl;
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus NonEmpty(const option::Option& option, bool /*error*/)
  {
    if ((option.arg != 0) && (option.arg[0] != 0))
    {
      return option::ARG_OK;
    }

    std::cout << "Error: " << option.name << " requires a non-empty argument!" << std::endl;
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus Numeric(const option::Option& option, bool /*error*/)
  {
    char* endptr = 0;
    if ((option.arg != 0) && strtol(option.arg, &endptr, 10)){};
    if ((endptr != option.arg) && (*endptr == 0))
    {
      return option::ARG_OK;
    }

    std::cout << "Error: " << option.name << " requires a long-int argument!" << std::endl;
    return option::ARG_ILLEGAL;
  }
  
};


class DomainPartition;

class ProblemManager : public ObjectManagerBase
{
public:
  explicit ProblemManager( const std::string& name,
                           ObjectManagerBase * const parent );
  ~ProblemManager();

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static std::string CatalogName() { return "ProblemManager"; }
  string getName() const override final
  {
    return ProblemManager::CatalogName();
  }
  ///@}

  virtual void Registration( dataRepository::ManagedGroup * const ) override;

  void ParseCommandLineInput( int & argc, char* argv[]);

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
  PhysicsSolverManager * m_physicsSolverManager;
  EventManager * m_eventManager;
};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
