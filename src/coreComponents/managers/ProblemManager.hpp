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

#include "EventManager.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "fileIO/schema/schemaUtilities.hpp"

namespace geosx
{

class PhysicsSolverManager;
class DomainPartition;

class ProblemManager : public dataRepository::Group
{
public:
  
  /**
   * @param name the name of this object manager
   * @param parent the parent Group
   */
  explicit ProblemManager( const std::string & name,
                           Group * const parent );

  /**
   * @brief Destructor, deletes all Groups and Wrappers owned by this Group
   */
  ~ProblemManager() override;

  /**
   * @brief Handles deviations between the datastructure and schema
   * @param schemaRoot schema root node handle
   * @param schemaParent schema parent node handle
   * @param documentationType flag to indicate the type of schema (0=input, 1=other)
   * 
   * This function handles deviations between the xml and data structure
   * on the Problem level (Functions, Mesh, etc.).  This can also be used to
   * add entries to the schema, which are not used during normal code execution
   * (e.g.: Benchmark)
   */
  virtual void SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                    xmlWrapper::xmlNode schemaParent,
                                    integer documentationType ) override;

  /**
   * @brief Creates a new sub-Group using the ObjectCatalog functionality.
   * @param childKey The name of the new object type's key in the
   *                 ObjectCatalog.
   * @param childName The name of the new object in the collection of
   *                  sub-Groups.
   * @return A pointer to the new Group created by this function.
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Parses command line input
   */
  void ParseCommandLineInput();

  /**
   * @brief Parses a restart file
   * @param restartFileName the name of the restart file
   */
  static bool ParseRestart( std::string & restartFileName );

  /**
   * @brief Initializes a python interpreter within GEOSX
   *
   * Note: This is not regularly used or tested, and may be removed in future versions.
   * To use this feature, the code must be compiled with the GEOSX_USE_PYTHON flag
   */
  void InitializePythonInterpreter();

  /**
   * @brief Closes the internal python interpreter
   *
   * Note: This is not regularly used or tested, and may be removed in future versions.
   * To use this feature, the code must be compiled with the GEOSX_USE_PYTHON flag
   */
  void ClosePythonInterpreter();

  /**
   * @brief Generates the xml schema documentation
   *
   * This function is called when the code is called with the -s schema_name option.
   * Before generating the schema, the code builds up a comprehensive datastructure.
   * (Note: catalog objects throughout the code will typically be registered via the
   * ExpandObjectCatalogs method.)  Once ready, SchemaUtilities will recusively walk
   * through the database, generating the xml schema.
   */
  void GenerateDocumentation();

  /**
   * @brief Parses the input xml file
   *
   * The name of the input file is indicated via the -i option on the command line
   */
  void ParseInputFile();

  /**
   * @brief Generates numerical meshes used throughout the code
   */
  void GenerateMesh();

  /**
   * @brief Applies numerical methods to objects throughout the code
   */
  void ApplyNumericalMethods();

  /**
   * @brief Defines the order in which objects should be initialized
   */
  void InitializationOrder( string_array & order ) override final;

  /**
   * @brief Sets up the problem after the input has been read in
   */
  void ProblemSetup();

  /**
   * @brief Run the events in the scheduler.
   */
  void RunSimulation();

  /**
   * @brief After initialization, overwrites data using a restart file
   */
  void ReadRestartOverwrite();

  /**
   * @brief Applies initial conditions indicated within the input file FieldSpecifications block
   */
  void ApplyInitialConditions();

  /**
   * @brief Returns a pointer to the DomainPartition
   */
  DomainPartition * getDomainPartition();

  /**
   * @brief Returns a pointer to the DomainPartition
   */
  DomainPartition const * getDomainPartition() const;

  /**
   * @brief Returns the problem name
   */
  const string & getProblemName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.problemName ); }

  /**
   * @brief Returns the input file name
   */
  const string & getInputFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.inputFileName ); }

  /**
   * @brief Returns the restart file name
   */
  const string & getRestartFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.restartFileName ); }

  /**
   * @brief Returns the schema file name
   */
  const string & getSchemaFileName() const
  { return GetGroup< Group >( groupKeys.commandLine )->getReference< string >( viewKeys.schemaFileName ); }

  /// Input file xml document handle
  xmlWrapper::xmlDocument xmlDocument;

  /// Input file parsing results
  xmlWrapper::xmlResult xmlResult;

  /// Input file Problem node handle
  xmlWrapper::xmlNode xmlProblemNode;

  /// Command line input viewKeys
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
    dataRepository::ViewKey useNonblockingMPI        = {"useNonblockingMPI"};
  } viewKeys;

  /// Child group viewKeys
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

  /**
   * @brief Returns the PhysicsSolverManager
   */
  PhysicsSolverManager & GetPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the PhysicsSolverManager
   */
  PhysicsSolverManager const & GetPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

protected:
  /**
   * @brief Post process the command line input
   */
  virtual void PostProcessInput() override final;

private:

  /// The PhysicsSolverManager
  PhysicsSolverManager * m_physicsSolverManager;
  
  /// The EventManager
  EventManager * m_eventManager;
  
  /// The FunctionManager
  FunctionManager * m_functionManager;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_PROBLEMMANAGER_HPP_ */
