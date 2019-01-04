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
 * ProblemManager.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#include "ProblemManager.hpp"

#include <vector>

#include "optionparser.h"

#include "DomainPartition.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "NumericalMethodsManager.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/Outputs/OutputManager.hpp"
#include "fileIO/utils/utils.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"
#include "MPI_Communications/SpatialPartition.hpp"
#include "meshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "dataRepository/SidreWrapper.hpp"
#include "dataRepository/RestartFlags.hpp"

#include "mesh/MeshBody.hpp"
#include "meshUtilities/MeshUtilities.hpp"
#include "common/TimingMacros.hpp"
// #include "managers/MeshLevel.hpp"
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;


struct Arg : public option::Arg
{
  static option::ArgStatus Unknown(const option::Option& option, bool /*error*/)
  {
    GEOS_LOG_RANK("Unknown option: " << option.name);
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus NonEmpty(const option::Option& option, bool /*error*/)
  {
    if ((option.arg != nullptr) && (option.arg[0] != 0))
    {
      return option::ARG_OK;
    }

    GEOS_LOG_RANK("Error: " << option.name << " requires a non-empty argument!");
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus Numeric(const option::Option& option, bool /*error*/)
  {
    char* endptr = nullptr;
    if ((option.arg != nullptr) && strtol(option.arg, &endptr, 10)) {};
    if ((endptr != option.arg) && (*endptr == 0))
    {
      return option::ARG_OK;
    }

    GEOS_LOG_RANK("Error: " << option.name << " requires a long-int argument!");
    return option::ARG_ILLEGAL;
  }

};


ProblemManager::ProblemManager( const std::string& name,
                                ManagedGroup * const parent ):
  ObjectManagerBase(name, parent),
  m_physicsSolverManager(nullptr),
  m_eventManager(nullptr),
  m_functionManager(nullptr)
{
  // Groups that do not read from the xml
  // RegisterGroup<DomainPartition>(groupKeys.domain)->BuildDataStructure(nullptr);
  RegisterGroup<DomainPartition>(groupKeys.domain);
  ManagedGroup * commandLine = RegisterGroup<ManagedGroup>(groupKeys.commandLine);
  commandLine->setRestartFlags(RestartFlags::WRITE);

  // Mandatory groups that read from the xml
  //RegisterGroup<BoundaryConditionManager>(groupKeys.boundaryConditionManager);
  // RegisterGroup<ConstitutiveManager>(groupKeys.constitutiveManager);
  // RegisterGroup<ElementRegionManager>(groupKeys.elementRegionManager);
  m_eventManager = RegisterGroup<EventManager>(groupKeys.eventManager);
  RegisterGroup<NumericalMethodsManager>(groupKeys.numericalMethodsManager);
  RegisterGroup<GeometricObjectManager>(groupKeys.geometricObjectManager);
  RegisterGroup<MeshManager>(groupKeys.meshManager);
  RegisterGroup<OutputManager>(groupKeys.outputManager);
  m_physicsSolverManager = RegisterGroup<PhysicsSolverManager>(groupKeys.physicsSolverManager);

  // The function manager is handled separately
  m_functionManager = NewFunctionManager::Instance();


  commandLine->RegisterViewWrapper<string>( viewKeys.inputFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the input xml file.");

  commandLine->RegisterViewWrapper<string>( viewKeys.restartFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the restart file.");

  commandLine->RegisterViewWrapper<integer>( viewKeys.beginFromRestart.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Flag to indicate restart run.");

  commandLine->RegisterViewWrapper<string>( viewKeys.problemName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Used in writing the output files, if not specified defaults to the name of the input file..");

  commandLine->RegisterViewWrapper<string>( viewKeys.outputDirectory.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Directory in which to put the output files, if not specified defaults to the current directory.");

  commandLine->RegisterViewWrapper<integer>( viewKeys.xPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the x-direction");

  commandLine->RegisterViewWrapper<integer>( viewKeys.yPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the y-direction");

  commandLine->RegisterViewWrapper<integer>( viewKeys.zPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the z-direction");

  commandLine->RegisterViewWrapper<integer>( viewKeys.overridePartitionNumbers.Key() )->
    setApplyDefaultValue(0)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Flag to indicate partition number override");

  commandLine->RegisterViewWrapper<string>( viewKeys.schemaFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the output schema");

  commandLine->RegisterViewWrapper<integer>( viewKeys.schemaLevel.Key() )->
    setApplyDefaultValue(0)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Schema verbosity level (0=default, 1=development, 2=all)");
}


ProblemManager::~ProblemManager()
{}


ManagedGroup * ProblemManager::CreateChild( string const & childKey, string const & childName )
{ return nullptr; }

void ProblemManager::RegisterDataOnMeshRecursive( ManagedGroup * const )
{
  ManagedGroup::RegisterDataOnMeshRecursive( GetGroup<DomainPartition>(groupKeys.domain)->getMeshBodies() );
}

void ProblemManager::ParseCommandLineInput( int argc, char** argv)
{
  ManagedGroup * commandLine = GetGroup<ManagedGroup>(groupKeys.commandLine);

  std::string& inputFileName = commandLine->getReference<std::string>(viewKeys.inputFileName);
  std::string& restartFileName = commandLine->getReference<std::string>(viewKeys.restartFileName);
  integer& beginFromRestart = commandLine->getReference<integer>(viewKeys.beginFromRestart);
  integer& xPartitionsOverride = commandLine->getReference<integer>(viewKeys.xPartitionsOverride);
  integer& yPartitionsOverride = commandLine->getReference<integer>(viewKeys.yPartitionsOverride);
  integer& zPartitionsOverride = commandLine->getReference<integer>(viewKeys.zPartitionsOverride);
  integer& overridePartitionNumbers = commandLine->getReference<integer>(viewKeys.overridePartitionNumbers);
  std::string&  schemaName = commandLine->getReference<std::string>(viewKeys.schemaFileName);
  integer& schemaLevel = commandLine->getReference<integer>(viewKeys.schemaLevel);
  schemaLevel = 0;
  std::string& problemName = commandLine->getReference<std::string>(viewKeys.problemName);
  std::string& outputDirectory = commandLine->getReference<std::string>(viewKeys.outputDirectory);
  outputDirectory = ".";
  problemName = "";


  // Set the options structs and parse
  enum optionIndex {UNKNOWN, HELP, INPUT, RESTART, XPAR, YPAR, ZPAR, SCHEMA, SCHEMALEVEL, PROBLEMNAME, OUTPUTDIR};
  const option::Descriptor usage[] =
  {
    {UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: geosx -i input.xml [options]\n\nOptions:"},
    {HELP, 0, "?", "help", Arg::None, "\t-?, --help"},
    {INPUT, 0, "i", "input", Arg::NonEmpty, "\t-i, --input, \t Input xml filename (required)"},
    {RESTART, 0, "r", "restart", Arg::NonEmpty, "\t-r, --restart, \t Target restart filename"},
    {XPAR, 0, "x", "xpartitions", Arg::Numeric, "\t-x, --x-partitions, \t Number of partitions in the x-direction"},
    {YPAR, 0, "y", "ypartitions", Arg::Numeric, "\t-y, --y-partitions, \t Number of partitions in the y-direction"},
    {ZPAR, 0, "z", "zpartitions", Arg::Numeric, "\t-z, --z-partitions, \t Number of partitions in the z-direction"},
    {SCHEMA, 0, "s", "schema", Arg::NonEmpty, "\t-s, --schema, \t Name of the output schema"},
    {SCHEMALEVEL, 0, "s", "schema_level", Arg::NonEmpty, "\t-s, --schema_level, \t Verbosity level of output schema (default=0)"},
    {PROBLEMNAME, 0, "n", "name", Arg::NonEmpty, "\t-n, --name, \t Name of the problem, used for output"},
    {OUTPUTDIR, 0, "o", "output", Arg::NonEmpty, "\t-o, --output, \t Directory to put the output files"},
    { 0, 0, nullptr, nullptr, nullptr, nullptr}
  };

  argc -= (argc>0);
  argv += (argc>0);
  option::Stats stats(usage, argc, argv);
  option::Option options[100];//stats.options_max];
  option::Option buffer[100];//stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);


  // Handle special cases
  if (parse.error())
  {
    GEOS_ERROR("Bad input arguments");
  }

  if (options[HELP] || (argc == 0))
  {
    int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
    option::printUsage(fwrite, stdout, usage, columns);
    exit(0);
  }

  if (options[INPUT].count() == 0)
  {
    if (options[SCHEMA].count() == 0)
    {
      GEOS_ERROR("An input xml must be specified!");
    }
  }


  // Iterate over the remaining inputs
  for (int ii=0 ; ii<parse.optionsCount() ; ++ii)
  {
    option::Option& opt = buffer[ii];
    switch (opt.index())
    {
    case UNKNOWN:
      // This should have thrown an error
      break;
    case HELP:
      // This is already handled above
      break;
    case INPUT:
      inputFileName = opt.arg;
      break;
    case RESTART:
      restartFileName = opt.arg;
      beginFromRestart = 1;
      break;
    case XPAR:
      xPartitionsOverride = std::stoi(opt.arg);
      overridePartitionNumbers = 1;
      break;
    case YPAR:
      yPartitionsOverride = std::stoi(opt.arg);
      overridePartitionNumbers = 1;
      break;
    case ZPAR:
      zPartitionsOverride = std::stoi(opt.arg);
      overridePartitionNumbers = 1;
      break;
    case SCHEMA:
      schemaName = opt.arg;
      break;
    case SCHEMALEVEL:
      schemaLevel = std::stoi(opt.arg);
      break;
    case PROBLEMNAME:
      problemName = opt.arg;
      break;
    case OUTPUTDIR:
      outputDirectory = opt.arg;
      break;
    }
  }

  if (schemaName.empty())
  {
    getAbsolutePath(inputFileName, inputFileName);

    if (problemName == "") 
    {
      if (inputFileName.length() > 4 && inputFileName.substr(inputFileName.length() - 4, 4) == ".xml")
      {
        string::size_type start = inputFileName.find_last_of('/') + 1;
        if (start >= inputFileName.length())
        {
          start = 0;
        }
        problemName.assign(inputFileName, start, inputFileName.length() - 4 - start);
      }
      else {
        problemName.assign(inputFileName);
      }
    }

    if (outputDirectory != ".")
    {
      mkdir(outputDirectory.data(), 0755);
      if (chdir(outputDirectory.data()) != 0)
      {
        GEOS_ERROR("Could not change to the ouput directory: " + outputDirectory);
      }
    }
  }
}


bool ProblemManager::ParseRestart( int argc, char** argv, std::string& restartFileName )
{
  // Set the options structs and parse
  enum optionIndex {UNKNOWN, HELP, INPUT, RESTART, XPAR, YPAR, ZPAR, SCHEMA, SCHEMALEVEL, PROBLEMNAME, OUTPUTDIR};
  const option::Descriptor usage[] =
  {
    {UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: geosx -i input.xml [options]\n\nOptions:"},
    {HELP, 0, "?", "help", Arg::None, "\t-?, --help"},
    {INPUT, 0, "i", "input", Arg::NonEmpty, "\t-i, --input, \t Input xml filename (required)"},
    {RESTART, 0, "r", "restart", Arg::NonEmpty, "\t-r, --restart, \t Target restart filename"},
    {XPAR, 0, "x", "xpartitions", Arg::Numeric, "\t-x, --x-partitions, \t Number of partitions in the x-direction"},
    {YPAR, 0, "y", "ypartitions", Arg::Numeric, "\t-y, --y-partitions, \t Number of partitions in the y-direction"},
    {ZPAR, 0, "z", "zpartitions", Arg::Numeric, "\t-z, --z-partitions, \t Number of partitions in the z-direction"},
    {SCHEMA, 0, "s", "schema", Arg::NonEmpty, "\t-s, --schema, \t Name of the output schema"},
    {SCHEMALEVEL, 0, "l", "schema_level", Arg::NonEmpty, "\t-l, --schema_level, \t Verbosity level of output schema (default=0)"},
    {PROBLEMNAME, 0, "n", "name", Arg::NonEmpty, "\t-n, --name, \t Name of the problem, used for output"},
    {OUTPUTDIR, 0, "o", "output", Arg::NonEmpty, "\t-o, --output, \t Directory to put the output files"},
    { 0, 0, nullptr, nullptr, nullptr, nullptr}
  };

  argc -= (argc>0); 
  argv += (argc>0);
  option::Stats stats(usage, argc, argv);
  option::Option options[100];//stats.options_max];
  option::Option buffer[100];//stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  
  // Handle special cases
  if (parse.error())
  {
    GEOS_ERROR("Bad input arguments");
  }

  if (options[HELP] || (argc == 0))
  {
    int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
    option::printUsage(fwrite, stdout, usage, columns);
    exit(0);
  }

  if (options[INPUT].count() == 0)
  {
    if (options[SCHEMA].count() == 0)
    {
      GEOS_ERROR("An input xml must be specified!");
    }
  }

  // Iterate over the remaining inputs
  bool beginFromRestart = false;
  for (int ii=0; ii<parse.optionsCount(); ++ii)
  {
    option::Option& opt = buffer[ii];
    switch (opt.index())
    {
      case UNKNOWN:
        break;
      case HELP:
        break;
      case INPUT:
        break;
      case RESTART:
        restartFileName = opt.arg;
        beginFromRestart = 1;
        break;
      case XPAR:
        break;
      case YPAR:
        break;
      case ZPAR:
        break;
      case SCHEMA:
        break;
      case SCHEMALEVEL:
        break;
      case PROBLEMNAME:
        break;
      case OUTPUTDIR:
        break;
    }
  }

  if (beginFromRestart == 1)
  {
    std::string dirname;
    std::string basename;
    splitPath(restartFileName, dirname, basename);

    std::vector<std::string> dir_contents;
    readDirectory(dirname, dir_contents);

    if (dir_contents.size() == 0)
    {
      GEOS_ERROR("Directory gotten from " << restartFileName << " " << dirname << " is empty.");
    }

    std::regex basename_regex(basename);

    std::string min_str("");
    std::string & max_match = min_str;
    bool match_found = false;
    for (std::string & s : dir_contents)
    {
      if (std::regex_match(s, basename_regex))
      {
        match_found = true;
        max_match = (s > max_match)? s : max_match;
      }
    }

    if (!match_found) {
      GEOS_ERROR("No matches found for pattern " << basename << " in directory " << dirname << ".");
    }

    restartFileName = dirname + "/" + max_match;
    getAbsolutePath(restartFileName, restartFileName);
  }

  return beginFromRestart;
}


void ProblemManager::InitializePythonInterpreter()
{  
#ifdef GEOSX_USE_PYTHON
  // Initialize python and numpy
  GEOS_LOG_RANK_0("Loading python interpreter");

  // Check to make sure the appropriate environment variables are set
  if (getenv("GPAC_SCHEMA") == NULL)
  {
    GEOS_ERROR("GPAC_SCHEMA must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONPATH") == NULL)
  {
    GEOS_ERROR("GEOS_PYTHONPATH must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONHOME") == NULL)
  {
    GEOS_ERROR("GEOS_PYTHONHOME must be defined to use the new preprocessor!");
  }

  setenv("PYTHONPATH", getenv("GEOS_PYTHONPATH"), 1);
  Py_SetPythonHome(getenv("GEOS_PYTHONHOME"));
  Py_Initialize();
  import_array();
#endif
}


void ProblemManager::ClosePythonInterpreter()
{
#ifdef GEOSX_USE_PYTHON
  // Add any other cleanup here
  GEOS_LOG_RANK_0("Closing python interpreter");
  Py_Finalize();
#endif
}


void ProblemManager::GenerateDocumentation()
{
  // Documentation output
  std::cout << "Trying to generate schema..." << std::endl;
  ManagedGroup * commandLine = GetGroup<ManagedGroup>(groupKeys.commandLine);
  std::string const & schemaName = commandLine->getReference<std::string>(viewKeys.schemaFileName);
  
  if (schemaName.empty() == 0)
  {
    GenerateDataStructureSkeleton(0);
  }
}


void ProblemManager::ParseInputFile()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition * domain  = getDomainPartition();

  ManagedGroup * commandLine = GetGroup<ManagedGroup>(groupKeys.commandLine);
  std::string const& inputFileName = commandLine->getReference<std::string>(viewKeys.inputFileName);


#ifdef GEOSX_USE_PYTHON
  // Load the pygeos module
  PyObject *pModule = PyImport_ImportModule("pygeos");
  if (pModule == NULL)
  {
    PyErr_Print();
    GEOS_ERROR("Could not find the pygeos module in GEOS_PYTHONPATH!");
  }

  // Call the xml preprocessor
  PyObject *pPreprocessorFunction = PyObject_GetAttrString(pModule, "PreprocessGEOSXML");
  PyObject *pPreprocessorInputStr = Py_BuildValue("(s)", inputFileName.c_str());
  PyObject *pKeywordDict = Py_BuildValue("{s:s}", "schema", getenv("GPAC_SCHEMA"));
  PyObject *pPreprocessorResult = PyObject_Call(pPreprocessorFunction, pPreprocessorInputStr, pKeywordDict);
  inputFileName = PyString_AsString(pPreprocessorResult);

  // Cleanup
  Py_DECREF(pPreprocessorResult);
  Py_DECREF(pKeywordDict);
  Py_DECREF(pPreprocessorInputStr);
  Py_DECREF(pPreprocessorFunction);
  Py_DECREF(pModule);

#else
  GEOS_LOG_RANK_0("GEOS must be configured to use Python to use parameters, symbolic math, etc. in input files");
#endif


  // Load preprocessed xml file and check for errors
  xmlResult = xmlDocument.load_file(inputFileName.c_str());
  if (!xmlResult)
  {
    GEOS_LOG_RANK_0("XML parsed with errors!");
    GEOS_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOS_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }
  xmlProblemNode = xmlDocument.child("Problem");

  ProcessInputFileRecursive( xmlProblemNode );
  ProcessInputFileRecursive_PostProcess();


  // The function manager is handled separately
  {
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child("Functions");
    m_functionManager->ProcessInputFileRecursive( topLevelNode );
    m_functionManager->ProcessInputFileRecursive_PostProcess(  );
  }

  {
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child("BoundaryConditions");
    BoundaryConditionManager * const bcManager = BoundaryConditionManager::get();
    bcManager->ProcessInputFileRecursive( topLevelNode );
    bcManager->ProcessInputFileRecursive_PostProcess();

  }

  // The objects in domain are handled separately for now
  {
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child("Constitutive");
    ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
    constitutiveManager->ProcessInputFileRecursive( topLevelNode );
    constitutiveManager->ProcessInputFileRecursive_PostProcess();

    // Open mesh levels
    MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
//    meshManager->ProcessInputFileRecursive()
    meshManager->GenerateMeshLevels(domain);

    topLevelNode = xmlProblemNode.child("ElementRegions");
    ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->ProcessInputFileRecursive_PostProcess();
  }

  // Documentation output
  // std::string const & schemaName = commandLine->getReference<std::string>(viewKeys.schemaFileName);
  // if (schemaName.empty() == 0)
  // {
  //   integer& schemaLevel = commandLine->getReference<integer>(viewKeys.schemaLevel);
  //   ConvertDocumentationToSchema(schemaName.c_str(), *(getDocumentationNode()), schemaLevel);
  // }
}


void ProblemManager::InitializationOrder( string_array & order )
{
  set<string> usedNames;


  {
    order.push_back(groupKeys.numericalMethodsManager.Key());
    usedNames.insert(groupKeys.numericalMethodsManager.Key());
  }

  {
    order.push_back(groupKeys.domain.Key());
    usedNames.insert(groupKeys.domain.Key());
  }

  {
    order.push_back(groupKeys.eventManager.Key());
    usedNames.insert(groupKeys.eventManager.Key());
  }

  for( auto const & subGroup : this->GetSubGroups() )
  {
    if( usedNames.count(subGroup.first) == 0 )
    {
      order.push_back(subGroup.first);
    }
  }
}



void ProblemManager::InitializePreSubGroups( ManagedGroup * const group )
{
  DomainPartition * domain  = getDomainPartition();

  ManagedGroup const * commandLine = GetGroup<ManagedGroup>(groupKeys.commandLine);
  integer const & xparCL = commandLine->getReference<integer>(viewKeys.xPartitionsOverride);
  integer const & yparCL = commandLine->getReference<integer>(viewKeys.yPartitionsOverride);
  integer const & zparCL = commandLine->getReference<integer>(viewKeys.zPartitionsOverride);

  PartitionBase & partition = domain->getReference<PartitionBase>(keys::partitionManager);
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
    partition.setPartitions( xpar,  ypar, zpar );
  }

  MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
  meshManager->GenerateMeshes(domain);

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    NodeManager * const nodeManager = (*mesh.second).group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();

    GeometricObjectManager * geometricObjects = this->GetGroup<GeometricObjectManager>(groupKeys.geometricObjectManager);

    MeshUtilities::GenerateNodesets( geometricObjects,
                                     nodeManager );
    nodeManager->ConstructGlobalToLocalMap();

  }
}


void ProblemManager::InitializePostSubGroups( ManagedGroup * const group )
{

  ObjectManagerBase::InitializePostSubGroups(nullptr);

  DomainPartition * domain  = getDomainPartition();

  ManagedGroup * const meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  FaceManager * const faceManager = meshLevel->getFaceManager();

  ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();

  NodeManager * nodeManager = meshLevel->getNodeManager();
  faceManager->BuildFaces( nodeManager, elementManager );

  EdgeManager * edgeManager = meshLevel->getEdgeManager();
  edgeManager->BuildEdges(faceManager, nodeManager );

  nodeManager->SetEdgeMaps( meshLevel->getEdgeManager() );
  nodeManager->SetFaceMaps( meshLevel->getFaceManager() );
  nodeManager->SetElementMaps( meshLevel->getElemManager() );

  domain->GenerateSets();
  domain->SetupCommunications();

  faceManager->SetIsExternal();
  edgeManager->SetIsExternal( faceManager );

}

void ProblemManager::RunSimulation()
{
  DomainPartition * domain  = getDomainPartition();
  m_eventManager->Run(domain);
}


void ProblemManager::ApplySchedulerEvent()
{}


DomainPartition * ProblemManager::getDomainPartition()
{
  return GetGroup<DomainPartition>(keys::domain);
}

DomainPartition const * ProblemManager::getDomainPartition() const
{
  return GetGroup<DomainPartition>(keys::domain);
}

void ProblemManager::ApplyInitialConditions()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition * domain = GetGroup<DomainPartition>(keys::domain);

  BoundaryConditionManager const * boundaryConditionManager = BoundaryConditionManager::get();

  boundaryConditionManager->ApplyInitialConditions( domain );

}

void ProblemManager::ReadRestartOverwrite( const std::string& restartFileName )
{
#ifdef GEOSX_USE_ATK
  this->prepareToRead();
  m_functionManager->prepareToRead();
  BoundaryConditionManager::get()->prepareToRead();
  SidreWrapper::loadExternalData(restartFileName, MPI_COMM_GEOSX);
  this->finishReading();
  m_functionManager->finishReading();
  BoundaryConditionManager::get()->finishReading();
#endif
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, ProblemManager, string const &, ManagedGroup * const )

} /* namespace geosx */
