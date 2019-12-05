/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProblemManager.hpp
 */

#ifndef GEOSX_MANAGERS_PROBLEMMANAGER_HPP_
#define GEOSX_MANAGERS_PROBLEMMANAGER_HPP_

#ifdef GEOSX_USE_PYTHON
// Note: the python header must be included first to avoid conflicting
// definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

#include "ObjectManagerBase.hpp"
#include "EventManager.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "fileIO/schema/SchemaUtilities.hpp"

namespace geosx
{

class PhysicsSolverManager;
class DomainPartition;

class ProblemManager : public ObjectManagerBase
{
public:
  explicit ProblemManager( const std::string& name,
                           Group * const parent );

  ~ProblemManager() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  const static string CatalogName() 
  { return "Problem"; }
  virtual const string getCatalogName() const override final
  { return ProblemManager::CatalogName(); }
  ///@}

  /**
   * This function is used to inform the schema generator of any
   * deviations between the xml and GEOS data structures.
   */
  virtual void SetSchemaDeviations(xmlWrapper::xmlNode schemaRoot,
                                   xmlWrapper::xmlNode schemaParent,
                                   integer documentationType) override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  void ParseCommandLineInput( int argc, char* argv[]);

  static bool ParseRestart( int argc, char* argv[], std::string& restartFileName );

  void InitializePythonInterpreter();

  void ClosePythonInterpreter();

  void GenerateDocumentation();

  void ParseInputFile();

  void GenerateMesh();

  void ApplyNumericalMethods();

  void InitializationOrder( string_array & order ) override final;

  /**
   * Function to setup the problem once the input has been read in, or the values
   * of the objects in the hierarchy have been sufficently set to generate a
   * mesh, etc.
   */
  void ProblemSetup();

  /**
   * Run the events in the scheduler.
   */
  void RunSimulation();


  void ReadRestartOverwrite();

  void ApplyInitialConditions();

  DomainPartition * getDomainPartition();
  DomainPartition const * getDomainPartition() const;

  const string & getProblemName() const
  { return GetGroup<Group>(groupKeys.commandLine)->getReference<string>(viewKeys.problemName); }

  const string & getInputFileName() const
  { return GetGroup<Group>(groupKeys.commandLine)->getReference<string>(viewKeys.inputFileName); }

  const string & getRestartFileName() const
  { return GetGroup<Group>(groupKeys.commandLine)->getReference<string>(viewKeys.restartFileName); }

  const string & getSchemaFileName() const
  { return GetGroup<Group>(groupKeys.commandLine)->getReference<string>(viewKeys.schemaFileName); }

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult;
  xmlWrapper::xmlNode xmlProblemNode;

  struct viewKeysStruct
  {
    dataRepository::ViewKey inputFileName            = {"inputFileName"};
    dataRepository::ViewKey restartFileName          = {"restartFileName"};
    dataRepository::ViewKey beginFromRestart         = {"beginFromRestart"};
    dataRepository::ViewKey xPartitionsOverride      = {"xPartitionsOverride"};
    dataRepository::ViewKey yPartitionsOverride      = {"yPartitionsOverride"};
    dataRepository::ViewKey zPartitionsOverride      = {"zPartitionsOverride"};
    dataRepository::ViewKey overridePartitionNumbers = {"overridePartitionNumbers"};
    dataRepository::ViewKey schemaFileName           = {"schemaFileName"};
    dataRepository::ViewKey problemName              = {"problemName"};
    dataRepository::ViewKey outputDirectory          = {"outputDirectory"};
  } viewKeys;

  struct groupKeysStruct
  {
//    constexpr auto eventManager="EventManager";
    dataRepository::GroupKey commandLine    = { "commandLine" };
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };
    dataRepository::GroupKey domain    = { "domain" };
    dataRepository::GroupKey eventManager = { "Events" };
    dataRepository::GroupKey fieldSpecificationManager = { "FieldSpecifications" };
    dataRepository::GroupKey functionManager = { "Functions" };
    dataRepository::GroupKey geometricObjectManager = { "Geometry" };
    dataRepository::GroupKey meshManager = { "Mesh" };
    dataRepository::GroupKey numericalMethodsManager = { "NumericalMethods" };
    dataRepository::GroupKey outputManager = { "Outputs" };
    dataRepository::GroupKey physicsSolverManager = { "Solvers" };
  } groupKeys;

  PhysicsSolverManager & GetPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  PhysicsSolverManager const & GetPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

protected:
  virtual void PostProcessInput() override final;

  virtual void InitializePostSubGroups( Group * const group ) override final;

private:

  PhysicsSolverManager * m_physicsSolverManager;
  EventManager * m_eventManager;
  FunctionManager * m_functionManager;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_PROBLEMMANAGER_HPP_ */
