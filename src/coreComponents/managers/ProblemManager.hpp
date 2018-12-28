/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
#include "../../../cxx-utilities/src/src/DocumentationNode.hpp"

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

//  explicit ProblemManager( const std::string& name,
//                           ManagedGroup * const parent,
//                           cxx_utilities::DocumentationNode * docNode );
  ~ProblemManager() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  const static string CatalogName() 
  { return "ProblemManager"; }
  virtual const string getCatalogName() const override final
  { return ProblemManager::CatalogName(); }
  ///@}

  virtual void RegisterDataOnMeshRecursive( ManagedGroup * const MeshBodies ) override final;

  virtual ManagedGroup * CreateChild( string const & childKey, string const & childName ) override;

  void ParseCommandLineInput( int argc, char* argv[]);

  static bool ParseRestart( int argc, char* argv[], std::string& restartFileName );

  void InitializePythonInterpreter();

  void ClosePythonInterpreter();

  void ParseInputFile();

  void InitializationOrder( string_array & order ) override final;

  void InitializePreSubGroups( ManagedGroup * const group ) override final;

  void InitializePostSubGroups( ManagedGroup * const group ) override final;

  void RunSimulation();

  void ApplySchedulerEvent();

  void WriteSilo( integer const cycleNumber, real64 const problemTime );

  // function to create and dump the restart file
  void WriteRestart( integer const cycleNumber );

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
    dataRepository::ViewKey problemName              = {"problemName"};
    dataRepository::ViewKey outputDirectory          = {"outputDirectory"};
  } viewKeys;

  struct groupKeysStruct
  {
    dataRepository::GroupKey domain    = { "domain" };
    dataRepository::GroupKey commandLine    = { "commandLine" };
    dataRepository::GroupKey boundaryConditionManager = { "BoundaryConditions" };
    dataRepository::GroupKey constitutiveManager = { "Constitutive" };
    dataRepository::GroupKey elementRegionManager = { "ElementRegions" };
    dataRepository::GroupKey eventManager = { "Events" };
    dataRepository::GroupKey numericalMethodsManager = { "NumericalMethods" };
    dataRepository::GroupKey geometricObjectManager = { "Geometry" };
    dataRepository::GroupKey meshManager = { "Mesh" };
    dataRepository::GroupKey physicsSolverManager = { "Solvers" };
    dataRepository::GroupKey outputManager = { "Outputs" };
  } groupKeys;

  PhysicsSolverManager & GetPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  PhysicsSolverManager const & GetPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

private:
  PhysicsSolverManager * m_physicsSolverManager;
  //SolverBase * m_physicsSolverManager;
  EventManager * m_eventManager;
  NewFunctionManager * m_functionManager;
};

} /* namespace geosx */

#endif /* COMPONENTS_CORE_SRC_MANAGERS_PROBLEMMANAGER_HPP_ */
