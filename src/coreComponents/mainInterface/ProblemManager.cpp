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

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

// Source includes
#include "ProblemManager.hpp"
#include "GeosxState.hpp"
#include "initialization.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/Path.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "dataRepository/KeyNames.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "events/tasks/TasksManager.hpp"
#include "events/EventManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "fileIO/Outputs/OutputManager.hpp"
#include "functions/FunctionManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/simpleGeometricObjects/GeometricObjectManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "schema/schemaUtilities.hpp"

// System includes
#include <vector>
#include <regex>

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ProblemManager::ProblemManager( conduit::Node & root ):
  dataRepository::Group( dataRepository::keys::ProblemManager, root ),
  m_physicsSolverManager( nullptr ),
  m_eventManager( nullptr ),
  m_functionManager( nullptr ),
  m_fieldSpecificationManager( nullptr )
{
  // Groups that do not read from the xml
  registerGroup< DomainPartition >( groupKeys.domain );
  Group & commandLine = registerGroup< Group >( groupKeys.commandLine );
  commandLine.setRestartFlags( RestartFlags::WRITE );

  setInputFlags( InputFlags::PROBLEM_ROOT );

  m_fieldSpecificationManager = &registerGroup< FieldSpecificationManager >( groupKeys.fieldSpecificationManager );

  m_eventManager = &registerGroup< EventManager >( groupKeys.eventManager );
  registerGroup< NumericalMethodsManager >( groupKeys.numericalMethodsManager );
  registerGroup< GeometricObjectManager >( groupKeys.geometricObjectManager );
  registerGroup< MeshManager >( groupKeys.meshManager );
  registerGroup< OutputManager >( groupKeys.outputManager );
  m_physicsSolverManager = &registerGroup< PhysicsSolverManager >( groupKeys.physicsSolverManager );
  registerGroup< TasksManager >( groupKeys.tasksManager );
  m_functionManager = &registerGroup< FunctionManager >( groupKeys.functionManager );

  // Command line entries
  commandLine.registerWrapper< string >( viewKeys.inputFileName.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Name of the input xml file." );

  commandLine.registerWrapper< string >( viewKeys.restartFileName.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Name of the restart file." );

  commandLine.registerWrapper< integer >( viewKeys.beginFromRestart.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Flag to indicate restart run." );

  commandLine.registerWrapper< string >( viewKeys.problemName.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Used in writing the output files, if not specified defaults to the name of the input file." );

  commandLine.registerWrapper< string >( viewKeys.outputDirectory.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Directory in which to put the output files, if not specified defaults to the current directory." );

  commandLine.registerWrapper< integer >( viewKeys.xPartitionsOverride.key() ).
    setApplyDefaultValue( 1 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Number of partitions in the x-direction" );

  commandLine.registerWrapper< integer >( viewKeys.yPartitionsOverride.key() ).
    setApplyDefaultValue( 1 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Number of partitions in the y-direction" );

  commandLine.registerWrapper< integer >( viewKeys.zPartitionsOverride.key() ).
    setApplyDefaultValue( 1 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Number of partitions in the z-direction" );

  commandLine.registerWrapper< integer >( viewKeys.overridePartitionNumbers.key() ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Flag to indicate partition number override" );

  commandLine.registerWrapper< string >( viewKeys.schemaFileName.key() ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Name of the output schema" );

  commandLine.registerWrapper< integer >( viewKeys.useNonblockingMPI.key() ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Whether to prefer using non-blocking MPI communication where implemented (results in non-deterministic DOF numbering)." );

  commandLine.registerWrapper< integer >( viewKeys.suppressPinned.key( ) ).
    setApplyDefaultValue( 0 ).
    setRestartFlags( RestartFlags::WRITE ).
    setDescription( "Whether to disallow using pinned memory allocations for MPI communication buffers." );

}

ProblemManager::~ProblemManager()
{}


Group * ProblemManager::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{ return nullptr; }


void ProblemManager::problemSetup()
{
  GEOSX_MARK_FUNCTION;
  postProcessInputRecursive();

  generateMesh();

//  initialize_postMeshGeneration();

  applyNumericalMethods();

  registerDataOnMeshRecursive( getDomainPartition().getMeshBodies() );

  initialize();

  importFields();
}


void ProblemManager::parseCommandLineInput()
{
  Group & commandLine = getGroup< Group >( groupKeys.commandLine );

  CommandLineOptions const & opts = getGlobalState().getCommandLineOptions();

  commandLine.getReference< string >( viewKeys.restartFileName ) = opts.restartFileName;
  commandLine.getReference< integer >( viewKeys.beginFromRestart ) = opts.beginFromRestart;
  commandLine.getReference< integer >( viewKeys.xPartitionsOverride ) = opts.xPartitionsOverride;
  commandLine.getReference< integer >( viewKeys.yPartitionsOverride ) = opts.yPartitionsOverride;
  commandLine.getReference< integer >( viewKeys.zPartitionsOverride ) = opts.zPartitionsOverride;
  commandLine.getReference< integer >( viewKeys.overridePartitionNumbers ) = opts.overridePartitionNumbers;
  commandLine.getReference< integer >( viewKeys.useNonblockingMPI ) = opts.useNonblockingMPI;
  commandLine.getReference< integer >( viewKeys.suppressPinned ) = opts.suppressPinned;

  string & outputDirectory = commandLine.getReference< string >( viewKeys.outputDirectory );
  outputDirectory = opts.outputDirectory;
  OutputBase::setOutputDirectory( outputDirectory );

  string & inputFileName = commandLine.getReference< string >( viewKeys.inputFileName );
  inputFileName = xmlWrapper::buildMultipleInputXML( opts.inputFileNames, outputDirectory );

  string & schemaName = commandLine.getReference< string >( viewKeys.schemaFileName );
  schemaName = opts.schemaName;

  string & problemName = commandLine.getReference< string >( viewKeys.problemName );
  problemName = opts.problemName;
  OutputBase::setFileNameRoot( problemName );

  if( schemaName.empty())
  {
    inputFileName = getAbsolutePath( inputFileName );
    Path::pathPrefix() = splitPath( inputFileName ).first;
  }

  if( opts.traceDataMigration )
  {
    chai::ArrayManager::getInstance()->enableCallbacks();
  }
  else
  {
    chai::ArrayManager::getInstance()->disableCallbacks();
  }
}


bool ProblemManager::parseRestart( string & restartFileName, CommandLineOptions const & options )
{
  bool const beginFromRestart = options.beginFromRestart;
  restartFileName = options.restartFileName;

  if( beginFromRestart == 1 )
  {
    string dirname, basename;
    std::tie( dirname, basename ) = splitPath( restartFileName );

    std::vector< string > dir_contents = readDirectory( dirname );

    GEOSX_THROW_IF( dir_contents.empty(),
                    "Directory gotten from " << restartFileName << " " << dirname << " is empty.",
                    InputError );

    std::regex basename_regex( basename );

    string min_str;
    string & max_match = min_str;
    bool match_found = false;
    for( string & s : dir_contents )
    {
      if( std::regex_match( s, basename_regex ))
      {
        match_found = true;
        max_match = (s > max_match)? s : max_match;
      }
    }

    GEOSX_THROW_IF( !match_found,
                    "No matches found for pattern " << basename << " in directory " << dirname << ".",
                    InputError );

    restartFileName = getAbsolutePath( dirname + "/" + max_match );
  }

  return beginFromRestart;
}


void ProblemManager::generateDocumentation()
{
  // Documentation output
  GEOSX_LOG_RANK_0( "Trying to generate schema..." );
  Group & commandLine = getGroup< Group >( groupKeys.commandLine );
  string const & schemaName = commandLine.getReference< string >( viewKeys.schemaFileName );

  if( !schemaName.empty() )
  {
    // Generate an extensive data structure
    generateDataStructureSkeleton( 0 );

    MeshManager & meshManager = this->getGroup< MeshManager >( groupKeys.meshManager );
    DomainPartition & domain = getDomainPartition();
    meshManager.generateMeshLevels( domain );

    registerDataOnMeshRecursive( domain.getMeshBodies() );

    // Generate schema
    schemaUtilities::ConvertDocumentationToSchema( schemaName.c_str(), this, 0 );

    // Generate non-schema documentation
    schemaUtilities::ConvertDocumentationToSchema((schemaName + ".other").c_str(), this, 1 );
  }
}


void ProblemManager::setSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
                                          xmlWrapper::xmlNode schemaParent,
                                          integer documentationType )
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child( "xsd:choice" );
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child( "xsd:choice" );
    targetChoiceNode.append_attribute( "minOccurs" ) = "0";
    targetChoiceNode.append_attribute( "maxOccurs" ) = "unbounded";
  }

  // These objects are handled differently during the xml read step,
  // so we need to explicitly add them into the schema structure
  DomainPartition & domain = getDomainPartition();

  m_functionManager->generateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( *m_functionManager, schemaRoot, targetChoiceNode, documentationType );

  m_fieldSpecificationManager->generateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( *m_fieldSpecificationManager, schemaRoot, targetChoiceNode, documentationType );

  ConstitutiveManager & constitutiveManager = domain.getGroup< ConstitutiveManager >( groupKeys.constitutiveManager );
  schemaUtilities::SchemaConstruction( constitutiveManager, schemaRoot, targetChoiceNode, documentationType );

  MeshManager & meshManager = this->getGroup< MeshManager >( groupKeys.meshManager );
  meshManager.generateMeshLevels( domain );
  ElementRegionManager & elementManager = domain.getMeshBody( 0 ).getBaseDiscretization().getElemManager();
  elementManager.generateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( elementManager, schemaRoot, targetChoiceNode, documentationType );


  // Add entries that are only used in the pre-processor
  Group & IncludedList = this->registerGroup< Group >( xmlWrapper::includedListTag );
  IncludedList.setInputFlags( InputFlags::OPTIONAL );

  Group & includedFile = IncludedList.registerGroup< Group >( xmlWrapper::includedFileTag );
  includedFile.setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  schemaUtilities::SchemaConstruction( IncludedList, schemaRoot, targetChoiceNode, documentationType );

  Group & parameterList = this->registerGroup< Group >( "Parameters" );
  parameterList.setInputFlags( InputFlags::OPTIONAL );

  Group & parameter = parameterList.registerGroup< Group >( "Parameter" );
  parameter.setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  parameter.registerWrapper< string >( "value" ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Input parameter definition for the preprocessor" );

  schemaUtilities::SchemaConstruction( parameterList, schemaRoot, targetChoiceNode, documentationType );

  Group & benchmarks = this->registerGroup< Group >( "Benchmarks" );
  benchmarks.setInputFlags( InputFlags::OPTIONAL );

  for( string const machineName : {"quartz", "lassen"} )
  {
    Group & machine = benchmarks.registerGroup< Group >( machineName );
    machine.setInputFlags( InputFlags::OPTIONAL );

    Group & run = machine.registerGroup< Group >( "Run" );
    run.setInputFlags( InputFlags::OPTIONAL );

    run.registerWrapper< string >( "name" ).setInputFlag( InputFlags::REQUIRED ).
      setDescription( "The name of this benchmark." );

    run.registerWrapper< int >( "nodes" ).setInputFlag( InputFlags::REQUIRED ).
      setDescription( "The number of nodes needed to run the benchmark." );

    run.registerWrapper< int >( "tasksPerNode" ).setInputFlag( InputFlags::REQUIRED ).
      setDescription( "The number of tasks per node to run the benchmark with." );

    run.registerWrapper< int >( "threadsPerTask" ).setInputFlag( InputFlags::OPTIONAL ).
      setDescription( "The number of threads per task to run the benchmark with." );

    run.registerWrapper< int >( "timeLimit" ).setInputFlag( InputFlags::OPTIONAL ).
      setDescription( "The time limit of the benchmark." );

    run.registerWrapper< string >( "args" ).setInputFlag( InputFlags::OPTIONAL ).
      setDescription( "Any extra command line arguments to pass to GEOSX." );

    run.registerWrapper< string >( "autoPartition" ).setInputFlag( InputFlags::OPTIONAL ).
      setDescription( "May be 'Off' or 'On', if 'On' partitioning arguments are created automatically. Default is Off." );

    run.registerWrapper< array1d< int > >( "strongScaling" ).setInputFlag( InputFlags::OPTIONAL ).
      setDescription( "Repeat the benchmark N times, scaling the number of nodes in the benchmark by these values." );
  }

  schemaUtilities::SchemaConstruction( benchmarks, schemaRoot, targetChoiceNode, documentationType );
}


void ProblemManager::parseInputFile()
{
  Group & commandLine = getGroup( groupKeys.commandLine );
  string const & inputFileName = commandLine.getReference< string >( viewKeys.inputFileName );

  // Load preprocessed xml file
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult const xmlResult = xmlDocument.load_file( inputFileName.c_str() );
  GEOSX_THROW_IF( !xmlResult, GEOSX_FMT( "Errors found while parsing XML file {}\nDescription: {}\nOffset: {}",
                                         inputFileName, xmlResult.description(), xmlResult.offset ), InputError );

  // Add path information to the file
  xmlDocument.append_child( xmlWrapper::filePathString ).append_attribute( xmlWrapper::filePathString ).set_value( inputFileName.c_str() );

  // Parse the results
  parseXMLDocument( xmlDocument );
}


void ProblemManager::parseInputString( string const & xmlString )
{
  // Load preprocessed xml file
  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( xmlString.c_str(), xmlString.length() );
  GEOSX_THROW_IF( !xmlResult, GEOSX_FMT( "Errors found while parsing XML string\nDescription: {}\nOffset: {}",
                                         xmlResult.description(), xmlResult.offset ), InputError );

  // Parse the results
  parseXMLDocument( xmlDocument );
}


void ProblemManager::parseXMLDocument( xmlWrapper::xmlDocument const & xmlDocument )
{
  // Extract the problem node and begin processing the user inputs
  xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( this->getName().c_str() );
  processInputFileRecursive( xmlProblemNode );

  // The objects in domain are handled separately for now
  {
    DomainPartition & domain = getDomainPartition();
    ConstitutiveManager & constitutiveManager = domain.getGroup< ConstitutiveManager >( groupKeys.constitutiveManager );
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager.getName().c_str());
    constitutiveManager.processInputFileRecursive( topLevelNode );

    // Open mesh levels
    MeshManager & meshManager = this->getGroup< MeshManager >( groupKeys.meshManager );
    meshManager.generateMeshLevels( domain );
    Group & meshBodies = domain.getMeshBodies();
    xmlWrapper::xmlNode elementRegionsNode = xmlProblemNode.child( MeshLevel::groupStructKeys::elemManagerString() );

    for( xmlWrapper::xmlNode regionNode : elementRegionsNode.children() )
    {
      string const regionName = regionNode.attribute( "name" ).value();
      string const
      regionMeshBodyName = ElementRegionBase::verifyMeshBodyName( meshBodies,
                                                                  regionNode.attribute( "meshBody" ).value() );

      MeshBody & meshBody = domain.getMeshBody( regionMeshBodyName );
      meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
      {
        ElementRegionManager & elementManager = meshLevel.getElemManager();
        Group * newRegion = elementManager.createChild( regionNode.name(), regionName );
        newRegion->processInputFileRecursive( regionNode );
      } );
    }
  }
}


void ProblemManager::postProcessInput()
{
  DomainPartition & domain = getDomainPartition();

  Group const & commandLine = getGroup< Group >( groupKeys.commandLine );
  integer const & xparCL = commandLine.getReference< integer >( viewKeys.xPartitionsOverride );
  integer const & yparCL = commandLine.getReference< integer >( viewKeys.yPartitionsOverride );
  integer const & zparCL = commandLine.getReference< integer >( viewKeys.zPartitionsOverride );

  integer const & suppressPinned = commandLine.getReference< integer >( viewKeys.suppressPinned );
  setPreferPinned((suppressPinned == 0));

  PartitionBase & partition = domain.getReference< PartitionBase >( keys::partitionManager );
  bool repartition = false;
  integer xpar = 1;
  integer ypar = 1;
  integer zpar = 1;
  if( xparCL != 0 )
  {
    repartition = true;
    xpar = xparCL;
  }
  if( yparCL != 0 )
  {
    repartition = true;
    ypar = yparCL;
  }
  if( zparCL != 0 )
  {
    repartition = true;
    zpar = zparCL;
  }
  if( repartition )
  {
    partition.setPartitions( xpar, ypar, zpar );
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
    // Case : Using MPI domain decomposition and partition are not defined (mainly for external mesh readers)
    if( mpiSize > 1 && xpar == 1 && ypar == 1 && zpar == 1 )
    {
      //TODO  confirm creates no issues with MPI_Cart_Create
      partition.setPartitions( 1, 1, mpiSize );
    }
  }
}


void ProblemManager::initializationOrder( string_array & order )
{
  SortedArray< string > usedNames;


  {
    order.emplace_back( groupKeys.numericalMethodsManager.key() );
    usedNames.insert( groupKeys.numericalMethodsManager.key() );
  }

  {
    order.emplace_back( groupKeys.domain.key() );
    usedNames.insert( groupKeys.domain.key() );
  }

  {
    order.emplace_back( groupKeys.eventManager.key() );
    usedNames.insert( groupKeys.eventManager.key() );
  }

  for( auto const & subGroup : this->getSubGroups() )
  {
    if( usedNames.count( subGroup.first ) == 0 )
    {
      order.emplace_back( subGroup.first );
    }
  }
}


void ProblemManager::generateMesh()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition & domain = getDomainPartition();

  MeshManager & meshManager = this->getGroup< MeshManager >( groupKeys.meshManager );


  meshManager.generateMeshes( domain );

  // get all the discretizations from the numerical methods.
  // map< pair< mesh body name, pointer to discretization>, array of region names >
  map< std::pair< string, Group const * const >, arrayView1d< string const > const >
  discretizations = getDiscretizations();

  // setup the base discretizations (hard code this for now)
  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {
    CellBlockManagerABC & cellBlockManager = meshBody.getGroup< CellBlockManagerABC >( keys::cellManager );

    MeshLevel & baseMesh = meshBody.getBaseDiscretization();
    array1d< string > junk;
    this->generateMeshLevel( baseMesh, cellBlockManager, nullptr, junk.toViewConst() );

    ElementRegionManager & elemManager = baseMesh.getElemManager();
    elemManager.generateWells( meshManager, baseMesh );

  } );

  // setup the MeshLevel assocaited with the discretizations
  for( auto const & discretizationPair: discretizations )
  {
    string const & meshBodyName = discretizationPair.first.first;
    MeshBody & meshBody = domain.getMeshBody( meshBodyName );

    if( discretizationPair.first.second!=nullptr ) // this check shouldn't be required
    {
      FiniteElementDiscretization const * const
      feDiscretization = dynamic_cast< FiniteElementDiscretization const * >( discretizationPair.first.second );

      // if the discretization is a finite element discretization
      if( feDiscretization != nullptr )
      {
        int const order = feDiscretization->getOrder();
        string const & discretizationName = feDiscretization->getName();
        arrayView1d< string const > const regionNames = discretizationPair.second;
        CellBlockManagerABC & cellBlockManager = meshBody.getGroup< CellBlockManagerABC >( keys::cellManager );

        // create a high order MeshLevel
        if( order > 1 )
        {
          MeshLevel & mesh = meshBody.createMeshLevel( MeshBody::groupStructKeys::baseDiscretizationString(),
                                                       discretizationName,
                                                       order );

          this->generateMeshLevel( mesh,
                                   cellBlockManager,
                                   feDiscretization,
                                   regionNames );
        }
        // Just create a shallow copy of the base discretization.
        else if( order==1 )
        {
          meshBody.createShallowMeshLevel( MeshBody::groupStructKeys::baseDiscretizationString(),
                                           discretizationName );
        }
      }
      else // this is a finite volume discretization...i hope
      {
        Group const * const discretization = discretizationPair.first.second;

        if( discretization != nullptr ) // ...it is FV if it isn't nullptr
        {
          string const & discretizationName = discretization->getName();
          meshBody.createShallowMeshLevel( MeshBody::groupStructKeys::baseDiscretizationString(),
                                           discretizationName );
        }
      }
    }
  }



  Group const & commandLine = this->getGroup< Group >( groupKeys.commandLine );
  integer const useNonblockingMPI = commandLine.getReference< integer >( viewKeys.useNonblockingMPI );
  domain.setupCommunications( useNonblockingMPI );

  domain.forMeshBodies( [&]( MeshBody & meshBody )
  {

    meshBody.deregisterGroup( keys::cellManager );

    meshBody.forMeshLevels( [&]( MeshLevel & meshLevel )
    {
      FaceManager & faceManager = meshLevel.getFaceManager();
      EdgeManager & edgeManager = meshLevel.getEdgeManager();

      faceManager.setIsExternal();
      edgeManager.setIsExternal( faceManager );
    } );
  } );

}


void ProblemManager::importFields()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition & domain = getDomainPartition();
  MeshManager & meshManager = this->getGroup< MeshManager >( groupKeys.meshManager );
  meshManager.importFields( domain );
}

void ProblemManager::applyNumericalMethods()
{

  DomainPartition & domain  = getDomainPartition();
  ConstitutiveManager & constitutiveManager = domain.getGroup< ConstitutiveManager >( groupKeys.constitutiveManager );
  Group & meshBodies = domain.getMeshBodies();

  // this contains a key tuple< mesh body name, mesh level name, region name, subregion name> with a value of the number of quadrature
  // points.
  map< std::tuple< string, string, string, string >, localIndex > const regionQuadrature = calculateRegionQuadrature( meshBodies );

  setRegionQuadrature( meshBodies, constitutiveManager, regionQuadrature );
}



map< std::pair< string, Group const * const >, arrayView1d< string const > const >
ProblemManager::getDiscretizations() const
{

  map< std::pair< string, Group const * const >, arrayView1d< string const > const > meshDiscretizations;

  NumericalMethodsManager const &
  numericalMethodManager = getGroup< NumericalMethodsManager >( groupKeys.numericalMethodsManager.key() );

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteVolumeManager const &
  fvDiscretizationManager = numericalMethodManager.getFiniteVolumeManager();

  DomainPartition const & domain  = getDomainPartition();
  Group const & meshBodies = domain.getMeshBodies();

  m_physicsSolverManager->forSubGroups< SolverBase >( [&]( SolverBase & solver )
  {

    solver.generateMeshTargetsFromTargetRegions( meshBodies );

    string const discretizationName = solver.getDiscretizationName();


    Group const *
    discretization = feDiscretizationManager.getGroupPointer( discretizationName );

    if( discretization==nullptr )
    {
      discretization = fvDiscretizationManager.getGroupPointer( discretizationName );
    }

    if( discretization!=nullptr )
    {
      solver.forDiscretizationOnMeshTargets( meshBodies,
                                             [&]( string const & meshBodyName,
                                                  MeshLevel const &,
                                                  auto const & regionNames )
      {
        std::pair< string, Group const * const > key = std::make_pair( meshBodyName, discretization );
        meshDiscretizations.insert( { key, regionNames } );
      } );
    }
  } );

  return meshDiscretizations;
}

void ProblemManager::generateMeshLevel( MeshLevel & meshLevel,
                                        CellBlockManagerABC & cellBlockManager,
                                        Group const * const discretization,
                                        arrayView1d< string const > const & )
{
  if( discretization != nullptr )
  {
    auto const * const
    feDisc = dynamic_cast< FiniteElementDiscretization const * >(discretization);

    auto const * const
    fvsDisc = dynamic_cast< FluxApproximationBase const * >(discretization);

    auto const * const
    fvhDisc = dynamic_cast< HybridMimeticDiscretization const * >(discretization);

    if( feDisc==nullptr && fvsDisc==nullptr && fvhDisc==nullptr )
    {
      GEOSX_ERROR( "Group expected to cast to a discretization object." );
    }
  }

  NodeManager & nodeManager = meshLevel.getNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

//  GeometricObjectManager & geometricObjects = this->getGroup< GeometricObjectManager >( groupKeys.geometricObjectManager );

//  MeshUtilities::generateNodesets( geometricObjects, nodeManager );


//  nodeManager.constructGlobalToLocalMap();


  if( meshLevel.getName() == MeshBody::groupStructKeys::baseDiscretizationString() )
  {
    elemManager.generateMesh( cellBlockManager );
    nodeManager.setGeometricalRelations( cellBlockManager, elemManager );
    edgeManager.setGeometricalRelations( cellBlockManager );
    faceManager.setGeometricalRelations( cellBlockManager,
                                         elemManager,
                                         nodeManager );
    nodeManager.constructGlobalToLocalMap( cellBlockManager );
    // Edge, face and element region managers rely on the sets provided by the node manager.
    // This is why `nodeManager.buildSets` is called first.
    nodeManager.buildSets( cellBlockManager, this->getGroup< GeometricObjectManager >( groupKeys.geometricObjectManager ) );
    edgeManager.buildSets( nodeManager );
    faceManager.buildSets( nodeManager );
    elemManager.buildSets( nodeManager );
    // The edge manager do not hold any information related to the regions nor the elements.
    // This is why the element region manager is not provided.
    nodeManager.setupRelatedObjectsInRelations( edgeManager, faceManager, elemManager );
    edgeManager.setupRelatedObjectsInRelations( nodeManager, faceManager );
    faceManager.setupRelatedObjectsInRelations( nodeManager, edgeManager, elemManager );
    // Node and edge managers rely on the boundary information provided by the face manager.
    // This is why `faceManager.setDomainBoundaryObjects` is called first.
    faceManager.setDomainBoundaryObjects();
    nodeManager.setDomainBoundaryObjects( faceManager );
    edgeManager.setDomainBoundaryObjects( faceManager );
  }
  meshLevel.generateSets();


  if( meshLevel.getName() == MeshBody::groupStructKeys::baseDiscretizationString() )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      subRegion.setupRelatedObjectsInRelations( meshLevel );
      subRegion.calculateElementGeometricQuantities( nodeManager, faceManager );
      subRegion.setMaxGlobalIndex();
    } );
    elemManager.setMaxGlobalIndex();
  }
}


map< std::tuple< string, string, string, string >, localIndex > ProblemManager::calculateRegionQuadrature( Group & meshBodies )
{

  NumericalMethodsManager const &
  numericalMethodManager = getGroup< NumericalMethodsManager >( groupKeys.numericalMethodsManager.key() );

  map< std::tuple< string, string, string, string >, localIndex > regionQuadrature;

  for( localIndex solverIndex=0; solverIndex<m_physicsSolverManager->numSubGroups(); ++solverIndex )
  {
    SolverBase const * const solver = m_physicsSolverManager->getGroupPointer< SolverBase >( solverIndex );

    if( solver != nullptr )
    {
      string const discretizationName = solver->getDiscretizationName();

      FiniteElementDiscretizationManager const &
      feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

      FiniteElementDiscretization const * const
      feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( discretizationName );

      solver->forDiscretizationOnMeshTargets( meshBodies,
                                              [&]( string const & meshBodyName,
                                                   MeshLevel & targetMeshLevel,
                                                   auto const & regionNames )
      {
        MeshLevel & baseMeshLevel = meshBodies.getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();

        MeshLevel & meshLevel = targetMeshLevel.isShallowCopyOf( baseMeshLevel ) ? baseMeshLevel : targetMeshLevel;

        NodeManager & nodeManager = meshLevel.getNodeManager();
        ElementRegionManager & elemManager = meshLevel.getElemManager();
        FaceManager const & faceManager = meshLevel.getFaceManager();
        EdgeManager const & edgeManager = meshLevel.getEdgeManager();
        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

        for( auto const & regionName : regionNames )
        {
          if( elemManager.hasRegion( regionName ) )
          {
            ElementRegionBase & elemRegion = elemManager.getRegion( regionName );

            if( feDiscretization != nullptr )
            {
              elemRegion.forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( auto & subRegion )
              {
                std::unique_ptr< finiteElement::FiniteElementBase > newFE = feDiscretization->factory( subRegion.getElementType() );

                finiteElement::FiniteElementBase &
                fe = subRegion.template registerWrapper< finiteElement::FiniteElementBase >( discretizationName, std::move( newFE ) ).
                       setRestartFlags( dataRepository::RestartFlags::NO_WRITE ).reference();
                subRegion.excludeWrappersFromPacking( { discretizationName } );

                finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe,
                                                                                         [&] ( auto & finiteElement )
                {
                  using FE_TYPE = std::remove_const_t< TYPEOFREF( finiteElement ) >;
                  using SUBREGION_TYPE = TYPEOFREF( subRegion );

                  typename FE_TYPE::template MeshData< SUBREGION_TYPE > meshData;
                  finiteElement::FiniteElementBase::initialize< FE_TYPE, SUBREGION_TYPE >( nodeManager,
                                                                                           edgeManager,
                                                                                           faceManager,
                                                                                           subRegion,
                                                                                           meshData );

                  localIndex const numQuadraturePoints = FE_TYPE::numQuadraturePoints;

#if !defined(GEOSX_USE_CUDA)
                  feDiscretization->calculateShapeFunctionGradients< SUBREGION_TYPE, FE_TYPE >( X, &subRegion, meshData, finiteElement );
#endif

                  localIndex & numQuadraturePointsInList = regionQuadrature[ std::make_tuple( meshBodyName,
                                                                                              meshLevel.getName(),
                                                                                              regionName,
                                                                                              subRegion.getName() ) ];

                  numQuadraturePointsInList = std::max( numQuadraturePointsInList, numQuadraturePoints );
                } );
              } );
            }
            else   //if( fvFluxApprox != nullptr )
            {
              elemRegion.forElementSubRegions( [&]( auto & subRegion )
              {
                localIndex & numQuadraturePointsInList = regionQuadrature[ std::make_tuple( meshBodyName,
                                                                                            meshLevel.getName(),
                                                                                            regionName,
                                                                                            subRegion.getName() ) ];
                localIndex const numQuadraturePoints = 1;
                numQuadraturePointsInList = std::max( numQuadraturePointsInList, numQuadraturePoints );
              } );
            }
          }
        }
      } );
    } // if( solver!=nullptr )
  }

  return regionQuadrature;
}


void ProblemManager::setRegionQuadrature( Group & meshBodies,
                                          ConstitutiveManager const & constitutiveManager,
                                          map< std::tuple< string, string, string, string >, localIndex > const & regionQuadratures )
{
  for( auto const & regionQuadrature : regionQuadratures )
  {
    std::tuple< string, string, string, string > const key = regionQuadrature.first;
    localIndex const numQuadraturePoints = regionQuadrature.second;
    string const meshBodyName = std::get< 0 >( key );
    string const meshLevelName = std::get< 1 >( key );
    string const regionName = std::get< 2 >( key );
    string const subRegionName = std::get< 3 >( key );

    GEOSX_LOG_RANK_0( "regionQuadrature: meshBodyName, meshLevelName, regionName, subRegionName = "<<
                      meshBodyName<<", "<<meshLevelName<<", "<<regionName<<", "<<subRegionName );



    MeshBody & meshBody = meshBodies.getGroup< MeshBody >( meshBodyName );
    MeshLevel & meshLevel = meshBody.getMeshLevel( meshLevelName );

//    if( meshLevel.isShallowCopy() )
    {
      ElementRegionManager & elemRegionManager = meshLevel.getElemManager();
      ElementRegionBase & elemRegion = elemRegionManager.getRegion( regionName );
      ElementSubRegionBase & elemSubRegion = elemRegion.getSubRegion( subRegionName );

      string_array const & materialList = elemRegion.getMaterialList();
      for( auto & materialName : materialList )
      {
        constitutiveManager.hangConstitutiveRelation( materialName, &elemSubRegion, numQuadraturePoints );
        GEOSX_LOG_RANK_0( GEOSX_FMT( "{}/{}/{}/{}/{} allocated {} quadrature points",
                                     meshBodyName,
                                     meshLevelName,
                                     regionName,
                                     subRegionName,
                                     materialName,
                                     numQuadraturePoints ) );
      }
    }
  }

  // check that every mesh element entity got a discretization
//  for( localIndex a = 0; a < meshBodies.getSubGroups().size(); ++a )
//  {
//    MeshBody & meshBody = meshBodies.getGroup< MeshBody >( a );
//    meshBody.forMeshLevels( [&] ( MeshLevel & meshLevel )
//    {
//      ElementRegionManager & elemManager = meshLevel.getElemManager();
//
//      elemManager.forElementSubRegionsComplete( [&]( localIndex const,
//                                                     localIndex const,
//                                                     ElementRegionBase & elemRegion,
//                                                     ElementSubRegionBase & elemSubRegion )
//      {
//        string const & regionName = elemRegion.getName();
//        string const & subRegionName = elemSubRegion.getName();
//
//        std::cout<<"looking for: meshBodyName, meshLevelName, regionName, subRegionName = "<<
//            meshBody.getName()<<", "<<meshLevel.getName()<<", "<<regionName<<", "<<subRegionName<<std::endl;
//
//
//        TYPEOFREF( regionQuadratures ) ::const_iterator rqIter = regionQuadratures.find( std::make_tuple( meshBody.getName(),
//                                                                                                        meshLevel.getName(),
//                                                                                                        regionName,
//                                                                                                        subRegionName ) );
//        GEOSX_ERROR_IF( rqIter == regionQuadratures.end(),
//                        GEOSX_FMT( "{}/{}/{}/{} does not have a discretization associated with it.",
//                                   meshBody.getName(),
//                                   meshLevel.getName(),
//                                   regionName,
//                                   subRegionName ) );
//      } );
//    } );
//  }
}

bool ProblemManager::runSimulation()
{
  return m_eventManager->run( getDomainPartition() );
}

DomainPartition & ProblemManager::getDomainPartition()
{
  return getGroup< DomainPartition >( groupKeys.domain );
}

DomainPartition const & ProblemManager::getDomainPartition() const
{
  return getGroup< DomainPartition >( groupKeys.domain );
}

void ProblemManager::applyInitialConditions()
{

  m_fieldSpecificationManager->forSubGroups< FieldSpecificationBase >( [&]( FieldSpecificationBase & fs )
  {
    fs.setMeshObjectPath( getDomainPartition().getMeshBodies() );
  } );

  getDomainPartition().forMeshBodies( [&] ( MeshBody & meshBody )
  {
    meshBody.forMeshLevels( [&] ( MeshLevel & meshLevel )
    {
      if( !meshLevel.isShallowCopy() ) // to avoid messages printed three times
      {
        m_fieldSpecificationManager->validateBoundaryConditions( meshLevel );
      }
      m_fieldSpecificationManager->applyInitialConditions( meshLevel );
    } );
  } );
  initializePostInitialConditions();
}

void ProblemManager::readRestartOverwrite()
{
  this->loadFromConduit();
  this->postRestartInitializationRecursive();
}

} /* namespace geosx */
