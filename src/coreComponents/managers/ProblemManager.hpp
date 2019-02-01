/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*
 * ProblemManager.hpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#ifndef COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_
#define COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_

#ifdef GEOSX_USE_PYTHON
// Note: the python header must be included first to avoid conflicting
// definitions of _posix_c_source
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

#include "ObjectManagerBase.hpp"
#include "EventManager.hpp"
#include "managers/Functions/NewFunctionManager.hpp"
#include "fileIO/schema/SchemaUtilities.hpp"

namespace geosx
{

class PhysicsSolverManager;
namespace dataRepository
{
namespace keys
{
string const eventManager="EventManager";
}
}


class DomainPartition;

class ProblemManager : public ObjectManagerBase
{
public:
  explicit ProblemManager( const std::string& name,
                           ManagedGroup * const parent );

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
                                   xmlWrapper::xmlNode schemaParent) override;

  virtual void RegisterDataOnMeshRecursive( ManagedGroup * const MeshBodies ) override final;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

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


  void ReadRestartOverwrite( const std::string& restartFileName );

  void ApplyInitialConditions();

  DomainPartition * getDomainPartition();
  DomainPartition const * getDomainPartition() const;

  const string & getProblemName() const
  { return GetGroup<ManagedGroup>(groupKeys.commandLine)->getReference<string>(viewKeys.problemName); }

  const string & getInputFileName() const
  { return GetGroup<ManagedGroup>(groupKeys.commandLine)->getReference<string>(viewKeys.inputFileName); }

  const string & getRestartFileName() const
  { return GetGroup<ManagedGroup>(groupKeys.commandLine)->getReference<string>(viewKeys.restartFileName); }

  const string & getSchemaFileName() const
  { return GetGroup<ManagedGroup>(groupKeys.commandLine)->getReference<string>(viewKeys.schemaFileName); }

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
    dataRepository::ViewKey schemaFileName           = {"schemaFileName"};
    dataRepository::ViewKey problemName              = {"problemName"};
    dataRepository::ViewKey outputDirectory          = {"outputDirectory"};
  } viewKeys;

  struct groupKeysStruct
  {
    dataRepository::GroupKey commandLine    = { "commandLine" };
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };
    dataRepository::GroupKey domain    = { "domain" };
    dataRepository::GroupKey elementRegionManager = { "ElementRegions" };
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

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override final;

private:

  PhysicsSolverManager * m_physicsSolverManager;
  EventManager * m_eventManager;
  NewFunctionManager * m_functionManager;
};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
