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

// Source includes
#include "ProblemManager.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/Path.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/initialization.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/Outputs/OutputManager.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "mesh/MeshBody.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "meshUtilities/MeshUtilities.hpp"
#include "meshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/SolverBase.hpp"

// System includes
#include <vector>
#include <regex>
#include <sys/stat.h>
#include <unistd.h>

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

class CellElementSubRegion;
class FaceElementSubRegion;


ProblemManager::ProblemManager( std::string const & name, conduit::Node & root ):
  dataRepository::Group( name, root ),
  m_physicsSolverManager( nullptr ),
  m_eventManager( nullptr ),
  m_functionManager( nullptr ),
  m_fieldSpecificationManager( nullptr )
{
  // Groups that do not read from the xml
  RegisterGroup< DomainPartition >( groupKeys.domain );
  Group * commandLine = RegisterGroup< Group >( groupKeys.commandLine );
  commandLine->setRestartFlags( RestartFlags::WRITE );

  setInputFlags( InputFlags::PROBLEM_ROOT );

  m_fieldSpecificationManager = RegisterGroup< FieldSpecificationManager >( groupKeys.fieldSpecificationManager );

  m_eventManager = RegisterGroup< EventManager >( groupKeys.eventManager );
  RegisterGroup< NumericalMethodsManager >( groupKeys.numericalMethodsManager );
  RegisterGroup< GeometricObjectManager >( groupKeys.geometricObjectManager );
  RegisterGroup< MeshManager >( groupKeys.meshManager );
  RegisterGroup< OutputManager >( groupKeys.outputManager );
  m_physicsSolverManager = RegisterGroup< PhysicsSolverManager >( groupKeys.physicsSolverManager );
  m_functionManager = RegisterGroup< FunctionManager >( groupKeys.functionManager );

  // Command line entries
  commandLine->registerWrapper< string >( viewKeys.inputFileName.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Name of the input xml file." );

  commandLine->registerWrapper< string >( viewKeys.restartFileName.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Name of the restart file." );

  commandLine->registerWrapper< integer >( viewKeys.beginFromRestart.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Flag to indicate restart run." );

  commandLine->registerWrapper< string >( viewKeys.problemName.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Used in writing the output files, if not specified defaults to the name of the input file." );

  commandLine->registerWrapper< string >( viewKeys.outputDirectory.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Directory in which to put the output files, if not specified defaults to the current directory." );

  commandLine->registerWrapper< integer >( viewKeys.xPartitionsOverride.Key() )->
    setApplyDefaultValue( 1 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Number of partitions in the x-direction" );

  commandLine->registerWrapper< integer >( viewKeys.yPartitionsOverride.Key() )->
    setApplyDefaultValue( 1 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Number of partitions in the y-direction" );

  commandLine->registerWrapper< integer >( viewKeys.zPartitionsOverride.Key() )->
    setApplyDefaultValue( 1 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Number of partitions in the z-direction" );

  commandLine->registerWrapper< integer >( viewKeys.overridePartitionNumbers.Key() )->
    setApplyDefaultValue( 0 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Flag to indicate partition number override" );

  commandLine->registerWrapper< string >( viewKeys.schemaFileName.Key() )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Name of the output schema" );

  commandLine->registerWrapper< integer >( viewKeys.useNonblockingMPI.Key() )->
    setApplyDefaultValue( 0 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Whether to prefer using non-blocking MPI communication where implemented (results in non-deterministic DOF numbering)." );

  commandLine->registerWrapper< integer >( viewKeys.suppressPinned.Key( ) )->
    setApplyDefaultValue( 0 )->
    setRestartFlags( RestartFlags::WRITE )->
    setDescription( "Whether to disallow using pinned memory allocations for MPI communication buffers." );

}

ProblemManager::~ProblemManager()
{}


Group * ProblemManager::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{ return nullptr; }


void ProblemManager::ProblemSetup()
{
  GEOSX_MARK_FUNCTION;
  PostProcessInputRecursive();

  GenerateMesh();

  ApplyNumericalMethods();

  RegisterDataOnMeshRecursive( GetGroup< DomainPartition >( groupKeys.domain )->getMeshBodies() );

  Initialize( this );
}


void ProblemManager::ParseCommandLineInput()
{
  Group * commandLine = GetGroup< Group >( groupKeys.commandLine );

  CommandLineOptions const & opts = getCommandLineOptions();

  commandLine->getReference< std::string >( viewKeys.restartFileName ) = opts.restartFileName;
  commandLine->getReference< integer >( viewKeys.beginFromRestart ) = opts.beginFromRestart;
  commandLine->getReference< integer >( viewKeys.xPartitionsOverride ) = opts.xPartitionsOverride;
  commandLine->getReference< integer >( viewKeys.yPartitionsOverride ) = opts.yPartitionsOverride;
  commandLine->getReference< integer >( viewKeys.zPartitionsOverride ) = opts.zPartitionsOverride;
  commandLine->getReference< integer >( viewKeys.overridePartitionNumbers ) = opts.overridePartitionNumbers;
  commandLine->getReference< integer >( viewKeys.useNonblockingMPI ) = opts.useNonblockingMPI;
  commandLine->getReference< integer >( viewKeys.suppressPinned ) = opts.suppressPinned;

  std::string & inputFileName = commandLine->getReference< std::string >( viewKeys.inputFileName );
  inputFileName = opts.inputFileName;

  std::string & schemaName = commandLine->getReference< std::string >( viewKeys.schemaFileName );
  schemaName = opts.schemaName;

  std::string & problemName = commandLine->getReference< std::string >( viewKeys.problemName );
  problemName = opts.problemName;

  std::string & outputDirectory = commandLine->getReference< std::string >( viewKeys.outputDirectory );
  outputDirectory = opts.outputDirectory;

  if( schemaName.empty())
  {
    getAbsolutePath( inputFileName, inputFileName );
    string xmlFolder;
    string notUsed;
    splitPath( inputFileName, xmlFolder, notUsed );
    Path::pathPrefix() = xmlFolder;

    if( outputDirectory != "." )
    {
      mkdir( outputDirectory.data(), 0755 );
      if( chdir( outputDirectory.data()) != 0 )
      {
        GEOSX_ERROR( "Could not change to the ouput directory: " + outputDirectory );
      }
    }
  }

  if( opts.suppressMoveLogging )
  {
    chai::ArrayManager::getInstance()->disableCallbacks();
  }
}


bool ProblemManager::ParseRestart( std::string & restartFileName )
{
  CommandLineOptions const & opts = getCommandLineOptions();
  bool const beginFromRestart = opts.beginFromRestart;
  restartFileName = opts.restartFileName;

  if( beginFromRestart == 1 )
  {
    std::string dirname;
    std::string basename;
    splitPath( restartFileName, dirname, basename );

    std::vector< std::string > dir_contents;
    readDirectory( dirname, dir_contents );

    GEOSX_ERROR_IF( dir_contents.size() == 0, "Directory gotten from " << restartFileName << " " << dirname << " is empty." );

    std::regex basename_regex( basename );

    std::string min_str( "" );
    std::string & max_match = min_str;
    bool match_found = false;
    for( std::string & s : dir_contents )
    {
      if( std::regex_match( s, basename_regex ))
      {
        match_found = true;
        max_match = (s > max_match)? s : max_match;
      }
    }

    GEOSX_ERROR_IF( !match_found, "No matches found for pattern " << basename << " in directory " << dirname << "." );

    restartFileName = dirname + "/" + max_match;
    getAbsolutePath( restartFileName, restartFileName );
  }

  return beginFromRestart;
}


void ProblemManager::GenerateDocumentation()
{
  // Documentation output
  std::cout << "Trying to generate schema..." << std::endl;
  Group * commandLine = GetGroup< Group >( groupKeys.commandLine );
  std::string const & schemaName = commandLine->getReference< std::string >( viewKeys.schemaFileName );

  if( schemaName.empty() == 0 )
  {
    // Generate an extensive data structure
    GenerateDataStructureSkeleton( 0 );

    MeshManager * meshManager = this->GetGroup< MeshManager >( groupKeys.meshManager );
    DomainPartition * domain  = getDomainPartition();
    meshManager->GenerateMeshLevels( domain );

    RegisterDataOnMeshRecursive( domain->getMeshBodies() );

    // Generate schema
    schemaUtilities::ConvertDocumentationToSchema( schemaName.c_str(), this, 0 );

    // Generate non-schema documentation
    schemaUtilities::ConvertDocumentationToSchema((schemaName + ".other").c_str(), this, 1 );
  }
}


void ProblemManager::SetSchemaDeviations( xmlWrapper::xmlNode schemaRoot,
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
  DomainPartition * domain  = getDomainPartition();

  m_functionManager->GenerateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( m_functionManager, schemaRoot, targetChoiceNode, documentationType );

  m_fieldSpecificationManager->GenerateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( m_fieldSpecificationManager, schemaRoot, targetChoiceNode, documentationType );

  ConstitutiveManager * constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  schemaUtilities::SchemaConstruction( constitutiveManager, schemaRoot, targetChoiceNode, documentationType );

  MeshManager * meshManager = this->GetGroup< MeshManager >( groupKeys.meshManager );
  meshManager->GenerateMeshLevels( domain );
  ElementRegionManager * elementManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  elementManager->GenerateDataStructureSkeleton( 0 );
  schemaUtilities::SchemaConstruction( elementManager, schemaRoot, targetChoiceNode, documentationType );


  // Add entries that are only used in the pre-processor
  Group * IncludedList = this->RegisterGroup< Group >( "Included" );
  IncludedList->setInputFlags( InputFlags::OPTIONAL );

  Group * includedFile = IncludedList->RegisterGroup< Group >( "File" );
  includedFile->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  schemaUtilities::SchemaConstruction( IncludedList, schemaRoot, targetChoiceNode, documentationType );

  Group * parameterList = this->RegisterGroup< Group >( "Parameters" );
  parameterList->setInputFlags( InputFlags::OPTIONAL );

  Group * parameter = parameterList->RegisterGroup< Group >( "Parameter" );
  parameter->setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
  parameter->registerWrapper< string >( "value" )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Input parameter definition for the preprocessor" );

  schemaUtilities::SchemaConstruction( parameterList, schemaRoot, targetChoiceNode, documentationType );

  Group * benchmarks = this->RegisterGroup< Group >( "Benchmarks" );
  benchmarks->setInputFlags( InputFlags::OPTIONAL );

  for( std::string const & machineName : {"quartz", "lassen"} )
  {
    Group * machine = benchmarks->RegisterGroup< Group >( machineName );
    machine->setInputFlags( InputFlags::OPTIONAL );

    Group * run = machine->RegisterGroup< Group >( "Run" );
    run->setInputFlags( InputFlags::OPTIONAL );

    run->registerWrapper< std::string >( "name" )->setInputFlag( InputFlags::REQUIRED )->
      setDescription( "The name of this benchmark." );

    run->registerWrapper< int >( "nodes" )->setInputFlag( InputFlags::REQUIRED )->
      setDescription( "The number of nodes needed to run the benchmark." );

    run->registerWrapper< int >( "tasksPerNode" )->setInputFlag( InputFlags::REQUIRED )->
      setDescription( "The number of tasks per node to run the benchmark with." );

    run->registerWrapper< int >( "threadsPerTask" )->setInputFlag( InputFlags::OPTIONAL )->
      setDescription( "The number of threads per task to run the benchmark with." );

    run->registerWrapper< int >( "timeLimit" )->setInputFlag( InputFlags::OPTIONAL )->
      setDescription( "The time limit of the benchmark." );

    run->registerWrapper< std::string >( "args" )->setInputFlag( InputFlags::OPTIONAL )->
      setDescription( "Any extra command line arguments to pass to GEOSX." );

    run->registerWrapper< std::string >( "autoPartition" )->setInputFlag( InputFlags::OPTIONAL )->
      setDescription( "May be 'Off' or 'On', if 'On' partitioning arguments are created automatically. Default is Off." );

    run->registerWrapper< array1d< int > >( "strongScaling" )->setInputFlag( InputFlags::OPTIONAL )->
      setDescription( "Repeat the benchmark N times, scaling the number of nodes in the benchmark by these values." );
  }

  schemaUtilities::SchemaConstruction( benchmarks, schemaRoot, targetChoiceNode, documentationType );
}


void ProblemManager::ParseInputFile()
{
  DomainPartition * domain  = getDomainPartition();

  Group * commandLine = GetGroup< Group >( groupKeys.commandLine );
  std::string const & inputFileName = commandLine->getReference< std::string >( viewKeys.inputFileName );

  // Load preprocessed xml file and check for errors
  xmlResult = xmlDocument.load_file( inputFileName.c_str());
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  string::size_type const pos=inputFileName.find_last_of( '/' );
  string path = inputFileName.substr( 0, pos + 1 );
  xmlDocument.append_child( xmlWrapper::filePathString ).append_attribute( xmlWrapper::filePathString ) = path.c_str();
  xmlProblemNode = xmlDocument.child( this->getName().c_str());
  ProcessInputFileRecursive( xmlProblemNode );

  // The objects in domain are handled separately for now
  {
    ConstitutiveManager * constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child( constitutiveManager->getName().c_str());
    constitutiveManager->ProcessInputFileRecursive( topLevelNode );

    // Open mesh levels
    MeshManager * meshManager = this->GetGroup< MeshManager >( groupKeys.meshManager );
    meshManager->GenerateMeshLevels( domain );
    ElementRegionManager * elementManager = domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
    topLevelNode = xmlProblemNode.child( elementManager->getName().c_str());
    elementManager->ProcessInputFileRecursive( topLevelNode );

  }
}


void ProblemManager::PostProcessInput()
{
  DomainPartition * domain  = getDomainPartition();

  Group const * commandLine = GetGroup< Group >( groupKeys.commandLine );
  integer const & xparCL = commandLine->getReference< integer >( viewKeys.xPartitionsOverride );
  integer const & yparCL = commandLine->getReference< integer >( viewKeys.yPartitionsOverride );
  integer const & zparCL = commandLine->getReference< integer >( viewKeys.zPartitionsOverride );

  integer const & suppressPinned = commandLine->getReference< integer >( viewKeys.suppressPinned );
  setPreferPinned((suppressPinned == 0));

  PartitionBase & partition = domain->getReference< PartitionBase >( keys::partitionManager );
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
    int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
    // Case : Using MPI domain decomposition and partition are not defined (mainly pamela usage)
    if( mpiSize > 1 && xpar == 1 && ypar == 1 && zpar == 1 )
    {
      //TODO  confirm creates no issues with MPI_Cart_Create
      partition.setPartitions( 1, 1, mpiSize );
    }
  }
}


void ProblemManager::InitializationOrder( string_array & order )
{
  SortedArray< string > usedNames;


  {
    order.emplace_back( groupKeys.numericalMethodsManager.Key() );
    usedNames.insert( groupKeys.numericalMethodsManager.Key() );
  }

  {
    order.emplace_back( groupKeys.domain.Key() );
    usedNames.insert( groupKeys.domain.Key() );
  }

  {
    order.emplace_back( groupKeys.eventManager.Key() );
    usedNames.insert( groupKeys.eventManager.Key() );
  }

  for( auto const & subGroup : this->GetSubGroups() )
  {
    if( usedNames.count( subGroup.first ) == 0 )
    {
      order.emplace_back( subGroup.first );
    }
  }
}


void ProblemManager::GenerateMesh()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition * domain  = getDomainPartition();

  MeshManager * meshManager = this->GetGroup< MeshManager >( groupKeys.meshManager );
  meshManager->GenerateMeshes( domain );
  Group * const cellBlockManager = domain->GetGroup( keys::cellManager );


  Group * const meshBodies = domain->getMeshBodies();

  for( localIndex a=0; a<meshBodies->numSubGroups(); ++a )
  {
    MeshBody * const meshBody = meshBodies->GetGroup< MeshBody >( a );
    for( localIndex b=0; b<meshBody->numSubGroups(); ++b )
    {
      MeshLevel * const meshLevel = meshBody->GetGroup< MeshLevel >( b );

      NodeManager * const nodeManager = meshLevel->getNodeManager();
      EdgeManager * edgeManager = meshLevel->getEdgeManager();
      FaceManager * const faceManager = meshLevel->getFaceManager();
      ElementRegionManager * const elemManager = meshLevel->getElemManager();

      GeometricObjectManager * geometricObjects = this->GetGroup< GeometricObjectManager >( groupKeys.geometricObjectManager );

      MeshUtilities::GenerateNodesets( geometricObjects,
                                       nodeManager );
      nodeManager->ConstructGlobalToLocalMap();

      elemManager->GenerateMesh( cellBlockManager );
      nodeManager->SetElementMaps( meshLevel->getElemManager() );

      faceManager->BuildFaces( nodeManager, elemManager );
      nodeManager->SetFaceMaps( meshLevel->getFaceManager() );

      edgeManager->BuildEdges( faceManager, nodeManager );
      nodeManager->SetEdgeMaps( meshLevel->getEdgeManager() );

      domain->GenerateSets();

      elemManager->forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
      {
        subRegion.setupRelatedObjectsInRelations( meshLevel );
        subRegion.CalculateElementGeometricQuantities( *nodeManager,
                                                       *faceManager );
      } );

      elemManager->GenerateCellToEdgeMaps( faceManager );

      elemManager->GenerateAggregates( faceManager, nodeManager );

      elemManager->GenerateWells( meshManager, meshLevel );
    }
  }

  GEOSX_ERROR_IF_NE( meshBodies->numSubGroups(), 1 );
  MeshBody * const meshBody = meshBodies->GetGroup< MeshBody >( 0 );

  GEOSX_ERROR_IF_NE( meshBody->numSubGroups(), 1 );
  MeshLevel * const meshLevel = meshBody->GetGroup< MeshLevel >( 0 );

  FaceManager * const faceManager = meshLevel->getFaceManager();
  EdgeManager * edgeManager = meshLevel->getEdgeManager();

  Group * commandLine = this->GetGroup< Group >( groupKeys.commandLine );
  integer const & useNonblockingMPI = commandLine->getReference< integer >( viewKeys.useNonblockingMPI );
  domain->SetupCommunications( useNonblockingMPI );
  faceManager->SetIsExternal();
  edgeManager->SetIsExternal( faceManager );
}


void ProblemManager::ApplyNumericalMethods()
{
  NumericalMethodsManager const * const
  numericalMethodManager = GetGroup< NumericalMethodsManager >( groupKeys.numericalMethodsManager.Key() );

  DomainPartition * domain  = getDomainPartition();
  ConstitutiveManager const * constitutiveManager = domain->GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );
  Group * const meshBodies = domain->getMeshBodies();

  map< string, localIndex > regionQuadrature;
  for( localIndex solverIndex=0; solverIndex<m_physicsSolverManager->numSubGroups(); ++solverIndex )
  {
    SolverBase const * const solver = m_physicsSolverManager->GetGroup< SolverBase >( solverIndex );

    string const numericalMethodName = solver->getDiscretization();
    arrayView1d< string const > const & targetRegions = solver->targetRegionNames();

    FiniteElementDiscretizationManager const &
    feDiscretizationManager = numericalMethodManager->getFiniteElementDiscretizationManager();

    FiniteElementDiscretization const *
      feDiscretization = feDiscretizationManager.GetGroup< FiniteElementDiscretization >( numericalMethodName );

    for( localIndex a=0; a<meshBodies->GetSubGroups().size(); ++a )
    {
      MeshBody * const meshBody = meshBodies->GetGroup< MeshBody >( a );
      for( localIndex b=0; b<meshBody->numSubGroups(); ++b )
      {
        MeshLevel * const meshLevel = meshBody->GetGroup< MeshLevel >( b );
        NodeManager * const nodeManager = meshLevel->getNodeManager();
        ElementRegionManager * const elemManager = meshLevel->getElemManager();
        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager->referencePosition();

        for( auto const & regionName : targetRegions )
        {
          ElementRegionBase * const elemRegion = elemManager->GetRegion( regionName );
          localIndex const quadratureSize = feDiscretization == nullptr ? 1 : feDiscretization->getNumberOfQuadraturePoints();
          if( quadratureSize > regionQuadrature[regionName] )
          {
            regionQuadrature[regionName] = quadratureSize;
          }
          elemRegion->forElementSubRegions< CellElementSubRegion,
                                            FaceElementSubRegion >( [&]( auto & subRegion )
          {
            if( feDiscretization != nullptr )
            {
              feDiscretization->CalculateShapeFunctionGradients( X, &subRegion );
            }
          } );
        }
      }
    }
  }

  for( localIndex a=0; a<meshBodies->GetSubGroups().size(); ++a )
  {
    MeshBody * const meshBody = meshBodies->GetGroup< MeshBody >( a );
    for( localIndex b=0; b<meshBody->numSubGroups(); ++b )
    {
      MeshLevel * const meshLevel = meshBody->GetGroup< MeshLevel >( b );
      ElementRegionManager * const elemManager = meshLevel->getElemManager();

      for( map< string, localIndex >::iterator iter=regionQuadrature.begin(); iter!=regionQuadrature.end(); ++iter )
      {
        string const regionName = iter->first;
        localIndex const quadratureSize = iter->second;

        ElementRegionBase * const elemRegion = elemManager->GetRegion( regionName );
        if( elemRegion != nullptr )
        {
          string_array const & materialList = elemRegion->getMaterialList();
          elemRegion->forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
          {
            for( auto & materialName : materialList )
            {
              constitutiveManager->HangConstitutiveRelation( materialName, &subRegion, quadratureSize );
            }
          } );
        }
      }
    }
  }
}

bool ProblemManager::RunSimulation()
{
  DomainPartition * domain  = getDomainPartition();
  return m_eventManager->Run( domain );
}

DomainPartition * ProblemManager::getDomainPartition()
{
  return GetGroup< DomainPartition >( keys::domain );
}

DomainPartition const * ProblemManager::getDomainPartition() const
{
  return GetGroup< DomainPartition >( keys::domain );
}

void ProblemManager::ApplyInitialConditions()
{
  DomainPartition * domain = GetGroup< DomainPartition >( keys::domain );
  m_fieldSpecificationManager->ApplyInitialConditions( domain );
  InitializePostInitialConditions( this );
}

void ProblemManager::ReadRestartOverwrite()
{
  this->loadFromConduit();
  this->postRestartInitializationRecursive( GetGroup< DomainPartition >( keys::domain ) );
}

} /* namespace geosx */
