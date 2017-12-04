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
// Note: the python header must be included first to avoid conflicting
// definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif
#include "optionparser.h"

#include "ObjectManagerBase.hpp"
#include "PhysicsSolvers/PhysicsSolverManager.hpp"
#include "EventManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "schema/SchemaUtilities.hpp"
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{
string const eventManager="EventManager";
}
}

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
    if ((option.arg != 0) && strtol(option.arg, &endptr, 10)) {};
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
                           ManagedGroup * const parent );

//  explicit ProblemManager( const std::string& name,
//                           ManagedGroup * const parent,
//                           cxx_utilities::DocumentationNode * docNode );
  ~ProblemManager();

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static std::string CatalogName() { return "ProblemManager"; }
  string getCatalogName() const override final
  {
    return ProblemManager::CatalogName();
  }
  ///@}



  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  virtual void CreateChild( string const & childKey, string const & childName ) override;

  void ParseCommandLineInput( int & argc, char* argv[]);

  void InitializePythonInterpreter();

  void ClosePythonInterpreter();

  void ParseInputFile();

  void InitializationOrder( string_array & order ) override final;

  void InitializePreSubGroups( ManagedGroup * const group ) override final;

  void InitializePostSubGroups( ManagedGroup * const group ) override final;

  void RunSimulation();

  void ApplySchedulerEvent();

  void WriteSilo( integer const cycleNumber, real64 const problemTime );

  void ApplyInitialConditions();

  DomainPartition * getDomainPartition();
  DomainPartition const * getDomainPartition() const;

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult;
  xmlWrapper::xmlNode xmlProblemNode;

  struct viewKeysStruct
  {
    dataRepository::ViewKey verbosity                = { "verbosityFlag" };
    dataRepository::ViewKey inputFileName            = {"inputFileName"};
    dataRepository::ViewKey restartFileName          = {"restartFileName"};
    dataRepository::ViewKey beginFromRestart         = {"beginFromRestart"};
    dataRepository::ViewKey xPartitionsOverride      = {"xPartitionsOverride"};
    dataRepository::ViewKey yPartitionsOverride      = {"yPartitionsOverride"};
    dataRepository::ViewKey zPartitionsOverride      = {"zPartitionsOverride"};
    dataRepository::ViewKey overridePartitionNumbers = {"overridePartitionNumbers"};
    dataRepository::ViewKey schemaLevel              = {"schemaLevel"};
  } viewKeys;

  struct groupKeysStruct
  {
    dataRepository::GroupKey domain    = { "domain" };
    dataRepository::GroupKey commandLine    = { "commandLine" };
    dataRepository::GroupKey boundaryConditionManager = { "BoundaryConditions" };
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };
    dataRepository::GroupKey elementRegionManager = { "ElementRegions" };
    dataRepository::GroupKey eventManager = { "Events" };
    dataRepository::GroupKey finiteElementManager = { "FiniteElements" };
    dataRepository::GroupKey geometricObjectManager = { "Geometry" };
    dataRepository::GroupKey meshManager = { "Mesh" };
    dataRepository::GroupKey physicsSolverManager = { "Solvers" };
  } groupKeys;



private:
  PhysicsSolverManager * m_physicsSolverManager;
  EventManager * m_eventManager;
  NewFunctionManager * m_functionManager;
};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
