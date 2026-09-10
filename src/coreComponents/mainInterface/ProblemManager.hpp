/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ProblemManager.hpp
 */


#ifndef GEOS_MAININTERFACE_PROBLEMMANAGER_HPP_
#define GEOS_MAININTERFACE_PROBLEMMANAGER_HPP_

#include "dataRepository/ProblemRepository.hpp"
#include "common/initializeEnvironment.hpp"

// for helper functions
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "functions/FunctionManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "events/EventManager.hpp"
#include "events/tasks/TasksManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/MeshManager.hpp"
#include "fileIO/Outputs/OutputManager.hpp"
#include "mesh/simpleGeometricObjects/GeometricObjectManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geos
{

/**
 * @brief A Group to contain the command line options within the data-repository
 */
class CommandLine : public dataRepository::Group
{
public:

  /// @cond DO_NOT_DOCUMENT

  struct viewKeysStruct
  {
    dataRepository::ViewKey inputFileName            = {"inputFileName"};            ///< Input file name key
    dataRepository::ViewKey restartFileName          = {"restartFileName"};          ///< Restart file name key
    dataRepository::ViewKey beginFromRestart         = {"beginFromRestart"};         ///< Flag to begin from restart key
    dataRepository::ViewKey xPartitionsOverride      = {"xPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in x key
    dataRepository::ViewKey yPartitionsOverride      = {"yPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in y key
    dataRepository::ViewKey zPartitionsOverride      = {"zPartitionsOverride"};      ///< Override of number of
                                                                                     ///< subdivisions in z key
    dataRepository::ViewKey overridePartitionNumbers = {"overridePartitionNumbers"}; ///< Flag to override partitioning
                                                                                     ///< key
    dataRepository::ViewKey schemaFileName           = {"schemaFileName"};           ///< Schema file name key
    dataRepository::ViewKey problemName              = {"problemName"};              ///< Problem name key
    dataRepository::ViewKey outputDirectory          = {"outputDirectory"};          ///< Output directory key
    dataRepository::ViewKey useNonblockingMPI        = {"useNonblockingMPI"};        ///< Flag to use non-block MPI key
    dataRepository::ViewKey suppressPinned           = {"suppressPinned"};           ///< Flag to suppress use of pinned
                                                                                     ///< memory key
  } m_vks; ///< Command line input viewKeys

  /// @endcond

  /**
   * @brief Construct a new CommandLine Group to contain the command line options within the data-repository.
   * @param name
   * @param parent
   */
  CommandLine( string const & name, Group * parent );

  /**
   * @brief Setup all the command line Group values from the provided inputs
   * @param options provided input options
   * @param inputFileName final composed main input filename (returned by `xmlWrapper::buildMultipleInputXML()`)
   */
  void setValues( CommandLineOptions const & options,
                  string_view inputFileName );

};

// CommandLine Group is available through the ProblemRepository as a mutable problem-unique manager.
template<> inline CommandLine & dataRepository::ProblemRepository::getManager()
{ return getRootGroup().getGroup< CommandLine >( m_gks.commandLine ); }

// CommandLine Group is available through the ProblemRepository as a const problem-unique manager.
template<> inline CommandLine const & dataRepository::ProblemRepository::getManager() const
{ return getRootGroup().getGroup< CommandLine >( m_gks.commandLine ); }

namespace constitutive
{ class ConstitutiveManager; }

/**
 * @class ProblemManager
 * @brief This is the class handling the operation flow of the problem being ran in GEOS
 */
class ProblemManager : public dataRepository::Group, public dataRepository::ProblemRepository
{
public:

  /**
   * @brief Create a new ProblemManager, it must be created from the root conduit node.
   * @param root The root conduit node.
   */
  explicit ProblemManager( conduit::Node & root );

  /**
   * @brief Destructor, deletes all Groups and Wrappers owned by this Group
   */
  ~ProblemManager() override;

  /**
   * @brief Handles deviations between the data structure and schema
   * @param schemaRoot schema root node handle
   * @param schemaParent schema parent node handle
   * @param documentationType flag to indicate the type of schema (0=input, 1=other)
   * @details This function handles deviations between the xml and data structure
   * on the Problem level (Functions, Mesh, etc.).  This can also be used to
   * add entries to the schema, which are not used during normal code execution
   * (e.g.: Benchmark)
   */
  virtual void setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
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
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Parses command line input
   */
  void parseCommandLineInput();

  /**
   * @brief Parses a restart file.
   * @param restartFileName The name of the restart file.
   * @param options The command line options.
   * @return Flag indicating beginFromRestart status
   */
  static bool parseRestart( string & restartFileName, CommandLineOptions const & options );

  /**
   * @brief Generates the xml schema documentation
   * This function is called when the code is called with the -s schema_name option.
   * @details Before generating the schema, the code builds up a comprehensive datastructure.
   * (Note: catalog objects throughout the code will typically be registered via the
   * ExpandObjectCatalogs method.)  Once ready, SchemaUtilities will recusively walk
   * through the database, generating the xml schema.
   */
  void generateDocumentation();

  /**
   * @brief Parses the input xml file
   * @details The name of the input file is indicated via the -i option on the command line
   */
  void parseInputFile();

  /**
   * @brief Parses the input xml string
   * @param xmlString the contents of the xml file as a string
   * @details This is used primarily for testing purposes
   */
  void parseInputString( string const & xmlString );

  /**
   * @brief Parses the input xml document. Also add the includes content to the xmlDocument when
   * `Include` nodes are encountered.
   * @param xmlDocument The parsed xml document handle
   */
  void parseXMLDocument( xmlWrapper::xmlDocument & xmlDocument );

  /**
   * @brief Generates numerical meshes used throughout the code
   */
  void generateMesh();

  /**
   * @brief Import field data from external sources (e.g. mesh generator).
   */
  void importFields();

  /**
   * @brief Allocates constitutive relations according to the discretizations
   *   on each subregion.
   */
  void applyNumericalMethods();

  /**
   * @brief Defines the order in which objects should be initialized
   * @param order list defining ordering sequence
   */
  void initializationOrder( string_array & order ) override final;

  /**
   * @brief Sets up the problem after the input has been read in
   */
  void problemSetup();

  /**
   * @brief Run the events in the scheduler.
   * @return True iff the simulation exited early, and needs to be run again to completion.
   */
  bool runSimulation();

  /**
   * @brief After initialization, overwrites data using a restart file
   */
  void readRestartOverwrite();

  /**
   * @brief Applies initial conditions indicated within the input file FieldSpecifications block
   */
  void applyInitialConditions();

  /**
   * @brief Returns the problem name
   * @return The problem name
   */
  string const & getProblemName() const;

  /**
   * @brief Returns the input file name
   * @return The input file name
   */
  string const & getInputFileName() const;

  /**
   * @brief Returns the restart file name
   * @return The restart file name
   */
  string const & getRestartFileName() const;

  /**
   * @brief Returns the schema file name
   * @return The schema file name
   */
  string const & getSchemaFileName() const;

  /**
   *  @name Managers access-helper functions for tests
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT

  CommandLine & getCommandLine()
  { return getManager< CommandLine >(); }

  CommandLine const & getCommandLine() const
  { return getManager< CommandLine >(); }

  DomainPartition & getDomainPartition()
  { return getManager< DomainPartition >(); }

  DomainPartition const & getDomainPartition() const
  { return getManager< DomainPartition >(); }

  PhysicsSolverManager & getPhysicsSolverManager()
  { return getManager< PhysicsSolverManager >(); }

  PhysicsSolverManager const & getPhysicsSolverManager() const
  { return getManager< PhysicsSolverManager >(); }

  FunctionManager & getFunctionManager()
  { return getManager< FunctionManager >(); }

  FunctionManager const & getFunctionManager() const
  { return getManager< FunctionManager >(); }

  FieldSpecificationManager & getFieldSpecificationManager()
  { return getManager< FieldSpecificationManager >(); }

  FieldSpecificationManager const & getFieldSpecificationManager() const
  { return getManager< FieldSpecificationManager >(); }

  EventManager & getEventManager()
  { return getManager< EventManager >(); }

  EventManager const & getEventManager() const
  { return getManager< EventManager >(); }

  TasksManager & getTasksManager()
  { return getManager< TasksManager >(); }

  TasksManager const & getTasksManager() const
  { return getManager< TasksManager >(); }

  NumericalMethodsManager & getNumericalMethodsManager()
  { return getManager< NumericalMethodsManager >(); }

  NumericalMethodsManager const & getNumericalMethodsManager() const
  { return getManager< NumericalMethodsManager >(); }

  MeshManager & getMeshManager()
  { return getManager< MeshManager >(); }

  MeshManager const & getMeshManager() const
  { return getManager< MeshManager >(); }

  OutputManager & getOutputManager()
  { return getManager< OutputManager >(); }

  OutputManager const & getOutputManager() const
  { return getManager< OutputManager >(); }

  GeometricObjectManager & getGeometricObjectManager()
  { return getManager< GeometricObjectManager >(); }

  GeometricObjectManager const & getGeometricObjectManager() const
  { return getManager< GeometricObjectManager >(); }

  constitutive::ConstitutiveManager & getConstitutiveManager()
  { return getManager< constitutive::ConstitutiveManager >(); }

  constitutive::ConstitutiveManager const & getConstitutiveManager() const
  { return getManager< constitutive::ConstitutiveManager >(); }

  /// @endcond
  ///@}

protected:
  /**
   * @brief Post process the command line input
   */
  virtual void postInputInitialization() override final;

private:

  /**
   * @brief Determine the number of quadrature points required for each
   *   MeshBody/Region/SubRegion.
   * @param meshBodies Reference to the mesh bodies object.
   * @return A tuple containing the number of quadrature points for every
   *   MeshBody/region/subregion combination.
   *
   * Checks all physics solvers for targetRegions and constitutive models to
   * determine the minimum number of quadrature points for each subregion.
   */
  map< std::tuple< string, string, string, string >, localIndex > calculateRegionQuadrature( Group & meshBodies );


  map< std::pair< string, Group const * const >, string_array const & > getDiscretizations();

  void generateMeshLevel( MeshLevel & meshLevel,
                          CellBlockManagerABC const & cellBlockManager,
                          Group const * const discretization,
                          string_array const & targetRegions );

  void generateMeshLevel( MeshLevel & meshLevel,
                          ParticleBlockManagerABC & particleBlockManager,
                          string_array const & );

  /**
   * @brief Allocate constitutive relations on each subregion with appropriate
   *   number of quadrature point.
   * @param meshBodies Reference to the mesh bodies object.
   * @param constitutiveManager The constitutive manager object.
   * @param regionQuadrature The map containing the number of quadrature points for every
   *  MeshBody/ElementRegion/ElementSubRegion.
   */
  void setRegionQuadrature( Group & meshBodies,
                            constitutive::ConstitutiveManager const & constitutiveManager,
                            map< std::tuple< string, string, string, string >, localIndex > const & regionQuadrature );

};

} /* namespace geos */

#endif /* GEOS_MAININTERFACE_PROBLEMMANAGER_HPP_ */
