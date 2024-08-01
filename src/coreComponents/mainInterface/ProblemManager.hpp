/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
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

#include "dataRepository/Group.hpp"
#include "dataRepository/InputProcessing.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshManager.hpp"

namespace geos
{

class PhysicsSolverManager;
class DomainPartition;
class GeometricObjectManager;
class FiniteElementDiscretization;
class MeshLevel;
namespace constitutive
{
class ConstitutiveManager;
}
class EventManager;
class TasksManager;
class FunctionManager;
class FieldSpecificationManager;
struct CommandLineOptions;
class CellBlockManagerABC;
class ParticleBlockManagerABC;

/**
 * @class ProblemManager
 * @brief This is the class handling the operation flow of the problem being ran in GEOS
 */
class ProblemManager : public dataRepository::Group
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

  void processSchemaDeviations( xmlWrapper::xmlDocument & document, std::set< string > & mergable );

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
   * @brief Returns a reference to the DomainPartition
   * @return reference to the DomainPartition
   */
  DomainPartition & getDomainPartition()
  { return getGroup< DomainPartition >( groupKeys.domain ); }

  /**
   * @brief Returns a reference to the DomainPartition
   * @return Const reference to the DomainPartition
   */
  DomainPartition const & getDomainPartition() const
  { return getGroup< DomainPartition >( groupKeys.domain ); }

  /**
   * @brief Returns the problem name
   * @return The problem name
   */
  string const & getProblemName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.problemName ); }

  /**
   * @brief Returns the input file name
   * @return The input file name
   */
  string const & getInputFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.inputFileName ); }

  /**
   * @brief Returns the restart file name
   * @return The restart file name
   */
  string const & getRestartFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.restartFileName ); }

  /**
   * @brief Returns the schema file name
   * @return The schema file name
   */
  string const & getSchemaFileName() const
  { return getGroup< Group >( groupKeys.commandLine ).getReference< string >( viewKeys.schemaFileName ); }

  /// Command line input viewKeys
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
    dataRepository::ViewKey preprocessOnly           = {"preprocessOnly"};           ///< preprocess on flag key
  } viewKeys; ///< Command line input viewKeys

  /// Child group viewKeys
  struct groupKeysStruct
  {
    /// @return Numerical methods string
    static constexpr char const * commandLineString() { return "commandLine"; }
    static constexpr char const * constitutiveManagerString() { return "Constitutive"; }
    static constexpr char const * domainString() { return "domain"; }
    static constexpr char const * eventManagerString() { return "Events"; }
    static constexpr char const * fieldSpecificationManagerString() { return "FieldSpecifications"; }
    static constexpr char const * functionManagerString() { return "Functions"; }
    static constexpr char const * geometricObjectManagerString() { return "Geometry"; }
    static constexpr char const * meshManagerString() { return "Mesh"; }
    static constexpr char const * numericalMethodsManagerString() { return "NumericalMethods"; }
    static constexpr char const * outputManagerString() { return "Outputs"; }
    static constexpr char const * physicsSolverManagerString() { return "Solvers"; }
    static constexpr char const * tasksManagerString() { return "Tasks"; }
    dataRepository::GroupKey commandLine    = { commandLineString() };                          ///< Command line key
    dataRepository::GroupKey constitutiveManager = { constitutiveManagerString() };                    ///< Constitutive key
    dataRepository::GroupKey domain    = { domainString() };                                    ///< Domain key
    dataRepository::GroupKey eventManager = { eventManagerString() };                                 ///< Events key
    dataRepository::GroupKey fieldSpecificationManager = { fieldSpecificationManagerString() };       ///< Field specification key
    dataRepository::GroupKey functionManager = { functionManagerString() };                           ///< Functions key
    dataRepository::GroupKey geometricObjectManager = { geometricObjectManagerString() };                     ///< Geometry key
    dataRepository::GroupKey meshManager = { meshManagerString() };                                    ///< Mesh key
    dataRepository::GroupKey numericalMethodsManager = { numericalMethodsManagerString()}; ///< Numerical methods key
    dataRepository::GroupKey outputManager = { outputManagerString() };                               ///< Outputs key
    dataRepository::GroupKey physicsSolverManager = { physicsSolverManagerString() };                        ///< Solvers key
    dataRepository::GroupKey tasksManager = { tasksManagerString() };                                  ///< Tasks key
  } groupKeys; ///< Child group viewKeys


  /**
   * @brief Returns the PhysicsSolverManager
   * @return Reference to the PhysicsSolverManager
   */
  PhysicsSolverManager & getPhysicsSolverManager()
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the PhysicsSolverManager
   * @return Const reference to the PhysicsSolverManager
   */
  PhysicsSolverManager const & getPhysicsSolverManager() const
  {
    return *m_physicsSolverManager;
  }

  /**
   * @brief Returns the FunctionManager.
   * @return The FunctionManager.
   */
  FunctionManager & getFunctionManager()
  {
    GEOS_ERROR_IF( m_functionManager == nullptr, "Not initialized." );
    return *m_functionManager;
  }

  /**
   * @brief Returns the const FunctionManager.
   * @return The const FunctionManager.
   */
  FunctionManager const & getFunctionManager() const
  {
    GEOS_ERROR_IF( m_functionManager == nullptr, "Not initialized." );
    return *m_functionManager;
  }

  /**
   * @brief Returns the FieldSpecificationManager.
   * @return The FieldSpecificationManager.
   */
  FieldSpecificationManager & getFieldSpecificationManager()
  {
    GEOS_ERROR_IF( m_fieldSpecificationManager == nullptr, "Not initialized." );
    return *m_fieldSpecificationManager;
  }

  /**
   * @brief Returns the const FunctionManager.
   * @return The const FunctionManager.
   */
  FieldSpecificationManager const & getFieldSpecificationManager() const
  {
    GEOS_ERROR_IF( m_fieldSpecificationManager == nullptr, "Not initialized." );
    return *m_fieldSpecificationManager;
  }

  /**
   * @brief Returns the EventManager.
   * @return The EventManager.
   */
  EventManager & getEventManager()
  {return *m_eventManager;}

  /**
   * @brief Returns the TasksManager.
   * @return The TasksManager.
   */
  TasksManager & getTasksManager()
  {return *m_tasksManager;}

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


  map< std::pair< string, Group const * const >, arrayView1d< string const > const >
  getDiscretizations() const;

  void generateMeshLevel( MeshLevel & meshLevel,
                          CellBlockManagerABC const & cellBlockManager,
                          Group const * const discretization,
                          arrayView1d< string const > const & targetRegions );

  void generateMeshLevel( MeshLevel & meshLevel,
                          ParticleBlockManagerABC & particleBlockManager,
                          arrayView1d< string const > const & );

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

  /// The PhysicsSolverManager
  PhysicsSolverManager * m_physicsSolverManager;

  /// The EventManager
  EventManager * m_eventManager;

  /// The TasksManager
  TasksManager * m_tasksManager;

  /// The FunctionManager
  FunctionManager * m_functionManager;

  /// The FieldSpecificationManager
  FieldSpecificationManager * m_fieldSpecificationManager;
};

} /* namespace geos */

#endif /* GEOS_MAININTERFACE_PROBLEMMANAGER_HPP_ */
