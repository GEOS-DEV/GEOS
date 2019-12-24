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


#include "ProblemManager.hpp"

#include <vector>
#include <regex>

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
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
#include "common/Path.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "dataRepository/ConduitRestart.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "mesh/MeshBody.hpp"
#include "wells/InternalWellGenerator.hpp"
#include "wells/WellElementRegion.hpp"
#include "meshUtilities/MeshUtilities.hpp"
#include "common/TimingMacros.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
// #include "managers/MeshLevel.hpp"
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

class CellElementSubRegion;
class FaceElementSubRegion;

struct Arg : public option::Arg
{
  static option::ArgStatus Unknown(const option::Option& option, bool /*error*/)
  {
    GEOSX_LOG_RANK("Unknown option: " << option.name);
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus NonEmpty(const option::Option& option, bool /*error*/)
  {
    if ((option.arg != nullptr) && (option.arg[0] != 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK("Error: " << option.name << " requires a non-empty argument!");
    return option::ARG_ILLEGAL;
  }


  static option::ArgStatus Numeric(const option::Option& option, bool /*error*/)
  {
    char* endptr = nullptr;
    if ((option.arg != nullptr) && strtol(option.arg, &endptr, 10)) {}
    if ((endptr != option.arg) && (*endptr == 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK("Error: " << option.name << " requires a long-int argument!");
    return option::ARG_ILLEGAL;
  }

};


ProblemManager::ProblemManager( const std::string& name,
                                Group * const parent ):
  ObjectManagerBase(name, parent),
  m_physicsSolverManager(nullptr),
  m_eventManager(nullptr),
  m_functionManager(nullptr)
{
  // Groups that do not read from the xml
  RegisterGroup<DomainPartition>(groupKeys.domain);
  Group * commandLine = RegisterGroup<Group>(groupKeys.commandLine);
  commandLine->setRestartFlags(RestartFlags::WRITE);

  setInputFlags(InputFlags::PROBLEM_ROOT);

  // Mandatory groups that read from the xml
  RegisterGroup<FieldSpecificationManager>( groupKeys.fieldSpecificationManager.Key(),
                                            &FieldSpecificationManager::get(),
                                            false );//->setRestartFlags(RestartFlags::NO_WRITE);


  // RegisterGroup<ConstitutiveManager>(groupKeys.constitutiveManager);
  // RegisterGroup<ElementRegionManager>(groupKeys.elementRegionManager);
  m_eventManager = RegisterGroup<EventManager>(groupKeys.eventManager);
  RegisterGroup<NumericalMethodsManager>(groupKeys.numericalMethodsManager);
  RegisterGroup<GeometricObjectManager>(groupKeys.geometricObjectManager);
  RegisterGroup<MeshManager>(groupKeys.meshManager);
  RegisterGroup<OutputManager>(groupKeys.outputManager);
  m_physicsSolverManager = RegisterGroup<PhysicsSolverManager>(groupKeys.physicsSolverManager);

  // The function manager is handled separately
  m_functionManager = &FunctionManager::Instance();
  // Mandatory groups that read from the xml
  RegisterGroup<FunctionManager>( groupKeys.functionManager.Key(),
                                  m_functionManager,
                                  false );

  // Command line entries
  commandLine->registerWrapper<string>( viewKeys.inputFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the input xml file.");

  commandLine->registerWrapper<string>( viewKeys.restartFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the restart file.");

  commandLine->registerWrapper<integer>( viewKeys.beginFromRestart.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Flag to indicate restart run.");

  commandLine->registerWrapper<string>( viewKeys.problemName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Used in writing the output files, if not specified defaults to the name of the input file..");

  commandLine->registerWrapper<string>( viewKeys.outputDirectory.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Directory in which to put the output files, if not specified defaults to the current directory.");

  commandLine->registerWrapper<integer>( viewKeys.xPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the x-direction");

  commandLine->registerWrapper<integer>( viewKeys.yPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the y-direction");

  commandLine->registerWrapper<integer>( viewKeys.zPartitionsOverride.Key() )->
    setApplyDefaultValue(1)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Number of partitions in the z-direction");

  commandLine->registerWrapper<integer>( viewKeys.overridePartitionNumbers.Key() )->
    setApplyDefaultValue(0)->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Flag to indicate partition number override");

  commandLine->registerWrapper<string>( viewKeys.schemaFileName.Key() )->
    setRestartFlags(RestartFlags::WRITE)->
    setDescription("Name of the output schema");
}


ProblemManager::~ProblemManager()
{}


Group * ProblemManager::CreateChild( string const & GEOSX_UNUSED_ARG( childKey ), string const & GEOSX_UNUSED_ARG( childName ) )
{ return nullptr; }


void ProblemManager::ProblemSetup()
{
  GEOSX_MARK_FUNCTION;
  PostProcessInputRecursive();

  GenerateMesh();

  ApplyNumericalMethods();

  RegisterDataOnMeshRecursive( GetGroup<DomainPartition>(groupKeys.domain)->getMeshBodies() );

  Initialize( this );

  ApplyInitialConditions();

  InitializePostInitialConditions( this );
}


void ProblemManager::ParseCommandLineInput( int argc, char** argv)
{
  Group * commandLine = GetGroup<Group>(groupKeys.commandLine);

  std::string& inputFileName = commandLine->getReference<std::string>(viewKeys.inputFileName);
  std::string& restartFileName = commandLine->getReference<std::string>(viewKeys.restartFileName);
  integer& beginFromRestart = commandLine->getReference<integer>(viewKeys.beginFromRestart);
  integer& xPartitionsOverride = commandLine->getReference<integer>(viewKeys.xPartitionsOverride);
  integer& yPartitionsOverride = commandLine->getReference<integer>(viewKeys.yPartitionsOverride);
  integer& zPartitionsOverride = commandLine->getReference<integer>(viewKeys.zPartitionsOverride);
  integer& overridePartitionNumbers = commandLine->getReference<integer>(viewKeys.overridePartitionNumbers);
  std::string& schemaName = commandLine->getReference<std::string>(viewKeys.schemaFileName);
  std::string& problemName = commandLine->getReference<std::string>(viewKeys.problemName);
  std::string& outputDirectory = commandLine->getReference<std::string>(viewKeys.outputDirectory);
  outputDirectory = ".";
  problemName = "";


  // Set the options structs and parse
  enum optionIndex {UNKNOWN, HELP, INPUT, RESTART, XPAR, YPAR, ZPAR, SCHEMA, PROBLEMNAME, OUTPUTDIR};
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
    GEOSX_ERROR("Bad input arguments");
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
      GEOSX_ERROR("An input xml must be specified!");
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
    string xmlFolder;
    string notUsed;
    splitPath( inputFileName, xmlFolder, notUsed );
    Path::pathPrefix() = xmlFolder;

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
        GEOSX_ERROR("Could not change to the ouput directory: " + outputDirectory);
      }
    }
  }
}


bool ProblemManager::ParseRestart( int argc, char** argv, std::string& restartFileName )
{
  // Set the options structs and parse
  enum optionIndex {UNKNOWN, HELP, INPUT, RESTART, XPAR, YPAR, ZPAR, SCHEMA, PROBLEMNAME, OUTPUTDIR};
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
    GEOSX_ERROR("Bad input arguments");
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
      GEOSX_ERROR("An input xml must be specified!");
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
      GEOSX_ERROR("Directory gotten from " << restartFileName << " " << dirname << " is empty.");
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
      GEOSX_ERROR("No matches found for pattern " << basename << " in directory " << dirname << ".");
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
  GEOSX_LOG_RANK_0("Loading python interpreter");

  // Check to make sure the appropriate environment variables are set
  if (getenv("GPAC_SCHEMA") == NULL)
  {
    GEOSX_ERROR("GPAC_SCHEMA must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONPATH") == NULL)
  {
    GEOSX_ERROR("GEOS_PYTHONPATH must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONHOME") == NULL)
  {
    GEOSX_ERROR("GEOS_PYTHONHOME must be defined to use the new preprocessor!");
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
  GEOSX_LOG_RANK_0("Closing python interpreter");
  Py_Finalize();
#endif
}


void ProblemManager::GenerateDocumentation()
{
  // Documentation output
  std::cout << "Trying to generate schema..." << std::endl;
  Group * commandLine = GetGroup<Group>(groupKeys.commandLine);
  std::string const & schemaName = commandLine->getReference<std::string>(viewKeys.schemaFileName);
  
  if (schemaName.empty() == 0)
  {
    // Generate an extensive data structure
    GenerateDataStructureSkeleton(0);

    MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
    DomainPartition * domain  = getDomainPartition();
    meshManager->GenerateMeshLevels(domain);

    RegisterDataOnMeshRecursive( domain->getMeshBodies() );

    // Generate schema
    SchemaUtilities::ConvertDocumentationToSchema(schemaName.c_str(), this, 0);

    // Generate non-schema documentation
    SchemaUtilities::ConvertDocumentationToSchema((schemaName + ".other").c_str(), this, 1);
  }
}


void ProblemManager::SetSchemaDeviations(xmlWrapper::xmlNode schemaRoot,
                                         xmlWrapper::xmlNode schemaParent,
                                         integer documentationType)
{
  xmlWrapper::xmlNode targetChoiceNode = schemaParent.child("xsd:choice");
  if( targetChoiceNode.empty() )
  {
    targetChoiceNode = schemaParent.prepend_child("xsd:choice");
    targetChoiceNode.append_attribute("minOccurs") = "0";
    targetChoiceNode.append_attribute("maxOccurs") = "unbounded";
  }

  // These objects are handled differently during the xml read step,
  // so we need to explicitly add them into the schema structure
  DomainPartition * domain  = getDomainPartition();

  m_functionManager->GenerateDataStructureSkeleton(0);
  SchemaUtilities::SchemaConstruction(m_functionManager, schemaRoot, targetChoiceNode, documentationType);

  FieldSpecificationManager & bcManager = FieldSpecificationManager::get();
  bcManager.GenerateDataStructureSkeleton(0);
  SchemaUtilities::SchemaConstruction(&bcManager, schemaRoot, targetChoiceNode, documentationType);

  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  SchemaUtilities::SchemaConstruction(constitutiveManager, schemaRoot, targetChoiceNode, documentationType);

  MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
  meshManager->GenerateMeshLevels(domain);
  ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
  elementManager->GenerateDataStructureSkeleton(0);
  SchemaUtilities::SchemaConstruction(elementManager, schemaRoot, targetChoiceNode, documentationType);


  // Add entries that are only used in the pre-processor
  Group * IncludedList = this->RegisterGroup<Group>("Included");
  IncludedList->setInputFlags(InputFlags::OPTIONAL);

  Group * includedFile = IncludedList->RegisterGroup<Group>("File");
  includedFile->setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);
  
  Group * parameterList = this->RegisterGroup<Group>("Parameters");
  parameterList->setInputFlags(InputFlags::OPTIONAL);

  Group * parameter = parameterList->RegisterGroup<Group>("Parameter");
  parameter->setInputFlags(InputFlags::OPTIONAL_NONUNIQUE);
  parameter->registerWrapper<string>("value")->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Input parameter definition for the preprocessor");

  SchemaUtilities::SchemaConstruction(IncludedList, schemaRoot, targetChoiceNode, documentationType);
  SchemaUtilities::SchemaConstruction(parameterList, schemaRoot, targetChoiceNode, documentationType);
}


void ProblemManager::ParseInputFile()
{
  DomainPartition * domain  = getDomainPartition();

  Group * commandLine = GetGroup<Group>(groupKeys.commandLine);
  std::string const& inputFileName = commandLine->getReference<std::string>(viewKeys.inputFileName);


#ifdef GEOSX_USE_PYTHON
  // Load the pygeos module
  PyObject *pModule = PyImport_ImportModule("pygeos");
  if (pModule == NULL)
  {
    PyErr_Print();
    GEOSX_ERROR("Could not find the pygeos module in GEOS_PYTHONPATH!");
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
  GEOSX_LOG_RANK_0("GEOS must be configured to use Python to use parameters, symbolic math, etc. in input files");
#endif


  // Load preprocessed xml file and check for errors
  xmlResult = xmlDocument.load_file(inputFileName.c_str());
  if (!xmlResult)
  {
    GEOSX_LOG_RANK_0("XML parsed with errors!");
    GEOSX_LOG_RANK_0("Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0("Error offset: " << xmlResult.offset);
  }

  string::size_type const pos=inputFileName.find_last_of('/');
  string path = inputFileName.substr( 0, pos + 1 );
  xmlDocument.append_child(xmlWrapper::filePathString).append_attribute(xmlWrapper::filePathString) = path.c_str();
  xmlProblemNode = xmlDocument.child(this->getName().c_str());
  ProcessInputFileRecursive( xmlProblemNode );

  // The objects in domain are handled separately for now
  {
    ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
    xmlWrapper::xmlNode topLevelNode = xmlProblemNode.child(constitutiveManager->getName().c_str());
    constitutiveManager->ProcessInputFileRecursive( topLevelNode );
    constitutiveManager->PostProcessInputRecursive();

    // Open mesh levels
    MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
    meshManager->GenerateMeshLevels(domain);
    ElementRegionManager * elementManager = domain->getMeshBody(0)->getMeshLevel(0)->getElemManager();
    topLevelNode = xmlProblemNode.child(elementManager->getName().c_str());
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->PostProcessInputRecursive();

  }
}


void ProblemManager::PostProcessInput()
{
  DomainPartition * domain  = getDomainPartition();

  Group const * commandLine = GetGroup<Group>(groupKeys.commandLine);
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
    partition.setPartitions( xpar, ypar, zpar );
    int const mpiSize = MpiWrapper::Comm_size(MPI_COMM_GEOSX) ;
    // Case : Using MPI domain decomposition and partition are not defined (mainly pamela usage)
    if( mpiSize > 1 && xpar == 1 && ypar == 1 && zpar == 1)
    {
      //TODO  confirm creates no issues with MPI_Cart_Create
      partition.setPartitions( 1,  1, mpiSize );
    }
  }
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


void ProblemManager::GenerateMesh()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition * domain  = getDomainPartition();

  MeshManager * meshManager = this->GetGroup<MeshManager>(groupKeys.meshManager);
  meshManager->GenerateMeshes(domain);
  Group const * const cellBlockManager = domain->GetGroup(keys::cellManager);


  Group * const meshBodies = domain->getMeshBodies();

  for( localIndex a=0; a<meshBodies->GetSubGroups().size() ; ++a )
  {
    MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(a);
    for( localIndex b=0 ; b<meshBody->numSubGroups() ; ++b )
    {
      MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(b);

      NodeManager * const nodeManager = meshLevel->getNodeManager();
      EdgeManager * edgeManager = meshLevel->getEdgeManager();
      FaceManager * const faceManager = meshLevel->getFaceManager();
      ElementRegionManager * const elemManager = meshLevel->getElemManager();

      GeometricObjectManager * geometricObjects = this->GetGroup<GeometricObjectManager>(groupKeys.geometricObjectManager);

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

      elemManager->forElementRegions( [&](ElementRegionBase * const region )->void
      {
        Group * subRegions = region->GetGroup(ElementRegionBase::viewKeyStruct::elementSubRegions);
        subRegions->forSubGroups<ElementSubRegionBase>( [&]( ElementSubRegionBase * const subRegion ) -> void
        {
          subRegion->setupRelatedObjectsInRelations( meshLevel );
          subRegion->CalculateElementGeometricQuantities( *nodeManager,
                                                          *faceManager );
        });

      });

      elemManager->GenerateAggregates( faceManager, nodeManager );

      elemManager->GenerateWells( meshManager, meshLevel );

    }
  }
}


void ProblemManager::ApplyNumericalMethods()
{
  NumericalMethodsManager const * const
  numericalMethodManager = GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  DomainPartition * domain  = getDomainPartition();
  ConstitutiveManager const * constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);
  Group * const meshBodies = domain->getMeshBodies();

  map<string,localIndex> regionQuadrature;
  for( localIndex solverIndex=0 ; solverIndex<m_physicsSolverManager->numSubGroups() ; ++solverIndex )
  {
    SolverBase const * const solver = m_physicsSolverManager->GetGroup<SolverBase>(solverIndex);

    string const numericalMethodName = solver->getDiscretization();
    string_array const & targetRegions = solver->getTargetRegions();

    FiniteElementDiscretizationManager const *
    feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(numericalMethodName);

    for( localIndex a=0; a<meshBodies->GetSubGroups().size() ; ++a )
    {
      MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(a);
      for( localIndex b=0 ; b<meshBody->numSubGroups() ; ++b )
      {
        MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(b);
        NodeManager * const nodeManager = meshLevel->getNodeManager();
        ElementRegionManager * const elemManager = meshLevel->getElemManager();
        arrayView1d<R1Tensor> const & X = nodeManager->referencePosition();

        for( auto const & regionName : targetRegions )
        {
          ElementRegionBase * const elemRegion = elemManager->GetRegion( regionName );
          localIndex const quadratureSize = feDiscretization == nullptr ? 1 : feDiscretization->getNumberOfQuadraturePoints();
          if( quadratureSize > regionQuadrature[regionName] )
          {
            regionQuadrature[regionName] = quadratureSize;
          }
          elemRegion->forElementSubRegions<CellElementSubRegion,
                                           FaceElementSubRegion>([&]( auto * const subRegion )->void
          {
            if( feDiscretization != nullptr )
            {
              feDiscretization->ApplySpaceToTargetCells(subRegion);
              feDiscretization->CalculateShapeFunctionGradients( X, subRegion);
            }
          });
        }
      }
    }
  }

  for( localIndex a=0; a<meshBodies->GetSubGroups().size() ; ++a )
  {
    MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(a);
    for( localIndex b=0 ; b<meshBody->numSubGroups() ; ++b )
    {
      MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(b);
      ElementRegionManager * const elemManager = meshLevel->getElemManager();

      for( map<string,localIndex>::iterator iter=regionQuadrature.begin() ; iter!=regionQuadrature.end() ; ++iter )
      {
        string const regionName = iter->first;
        localIndex const quadratureSize = iter->second;

        ElementRegionBase * const elemRegion = elemManager->GetRegion( regionName );
        if( elemRegion != nullptr )
        {
          string_array const & materialList = elemRegion->getMaterialList();
          elemRegion->forElementSubRegions([&]( auto * const subRegion )->void
          {
            for( auto & materialName : materialList )
            {
              constitutiveManager->HangConstitutiveRelation( materialName, subRegion, quadratureSize );
            }
          });
        }
      }
    }
  }
}


void ProblemManager::InitializePostSubGroups( Group * const GEOSX_UNUSED_ARG( group ) )
{

//  ObjectManagerBase::InitializePostSubGroups(nullptr);
//
  DomainPartition * domain  = getDomainPartition();

  Group * const meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);

  FaceManager * const faceManager = meshLevel->getFaceManager();
  EdgeManager * edgeManager = meshLevel->getEdgeManager();

  domain->SetupCommunications();
  faceManager->SetIsExternal();
  edgeManager->SetIsExternal( faceManager );
}

void ProblemManager::RunSimulation()
{
  GEOSX_MARK_FUNCTION;
  DomainPartition * domain  = getDomainPartition();
  m_eventManager->Run(domain);
}

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
  DomainPartition * domain = GetGroup<DomainPartition>(keys::domain);
  FieldSpecificationManager::get().ApplyInitialConditions( domain );
}

void ProblemManager::ReadRestartOverwrite()
{
  this->loadFromConduit();
  this->postRestartInitializationRecursive( GetGroup< DomainPartition >( keys::domain ) );
}



REGISTER_CATALOG_ENTRY( ObjectManagerBase, ProblemManager, string const &, Group * const )

} /* namespace geosx */
