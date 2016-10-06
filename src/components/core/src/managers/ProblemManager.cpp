/*
 * ProblemManager.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#include "ProblemManager.hpp"

#include "DomainPartition.hpp"
#include "PhysicsSolvers/SolverBase.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "finiteElement/FiniteElementSpace.hpp"
#include <stdexcept>


namespace geosx
{

namespace dataRepository
{
  namespace keys
  {
    std::string const commandLine = "commandLine";
    std::string const inputFileName = "inputFileName";
    std::string const restartFileName = "restartFileName";
    std::string const beginFromRestart = "beginFromRestart";
    std::string const xPartitionsOverride = "xPartitionsOverride";
    std::string const yPartitionsOverride = "yPartitionsOverride";
    std::string const zPartitionsOverride = "zPartitionsOverride";
    std::string const overridePartitionNumbers = "overridePartitionNumbers";
    std::string const solverApplications = "solverApplications";
    std::string const solverApplicationNames = "solverApplicationNames";

    std::string const K = "K";

  }
}

using namespace dataRepository;

ProblemManager::ProblemManager( const std::string& name,
                                ObjectManagerBase * const parent ) :
  ObjectManagerBase( name, parent ),
  m_physicsSolverManager(nullptr)
{
  m_physicsSolverManager = &(RegisterGroup<PhysicsSolverManager>("PhysicsSolverManager" ) ) ;
}

ProblemManager::ProblemManager( const std::string& name,
                                ObjectManagerBase * const parent,
                                cxx_utilities::DocumentationNode * docNode ) :
  ObjectManagerBase( name, parent, docNode ),
  m_physicsSolverManager(nullptr)
{
//  allocateDocumentationNode( "ProblemManager",
//                             "ProblemManager",
//                             0,
//                             "DocumentationNode",
//                             "UniqueNode",
//                             "This is the top level node in the input structure.",
//                             "This is the top level node in the input structure.",
//                             "",
//                             "ProblemManager",
//                             0,
//                             0,
//                             0,
//                             nullptr );
  m_physicsSolverManager = &(RegisterGroup<PhysicsSolverManager>("PhysicsSolverManager" ) ) ;

  m_eventManager = &(RegisterGroup<EventManager>("EventManager" ) ) ;
}

ProblemManager::~ProblemManager()
{
}



void ProblemManager::BuildDataStructure( dataRepository::ManagedGroup * const )
{

//  cxx_utilities::DocumentationNode newNode;

//  newNode.m_name = keys::inputFileName;

//  getDocumentationNode()->

  RegisterGroup<DomainPartition>(keys::domain);

  ManagedGroup& solverApplications = RegisterGroup<ManagedGroup>(keys::solverApplications);


  solverApplications.RegisterViewWrapper<string_array>(keys::solverApplicationNames);

  ManagedGroup& commandLine = RegisterGroup<ManagedGroup >(keys::commandLine);
//  commandLine.RegisterViewWrapper<std::string>(keys::inputFileName);
//  commandLine.RegisterViewWrapper<std::string>(keys::restartFileName);
//  commandLine.RegisterViewWrapper<bool>(keys::beginFromRestart);
//  commandLine.RegisterViewWrapper<int32>(keys::xPartitionsOverride);
//  commandLine.RegisterViewWrapper<int32>(keys::yPartitionsOverride);
//  commandLine.RegisterViewWrapper<int32>(keys::zPartitionsOverride);
//  commandLine.RegisterViewWrapper<bool>(keys::overridePartitionNumbers);
//
//  commandLine.RegisterViewWrapper<int32>(keys::K);

}


void ProblemManager::FillDocumentationNode( dataRepository::ManagedGroup * const group )
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  ObjectManagerBase::FillDocumentationNode( group );

  docNode->AllocateChildNode( keys::inputFileName,
                              keys::inputFileName,
                              -1,
                              "string",
                              "string",
                              "Name of the input xml file.",
                              "Name of the input xml file.",
                              "input.xml",
                              "ProblemManager",
                              0,
                              0 );

  ManagedGroup& commandLine = RegisterGroup<ManagedGroup >(keys::commandLine);
  commandLine.RegisterViewWrapper<std::string>(keys::inputFileName);
  commandLine.RegisterViewWrapper<std::string>(keys::restartFileName);
  commandLine.RegisterViewWrapper<bool>(keys::beginFromRestart);
  commandLine.RegisterViewWrapper<int32>(keys::xPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::yPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::zPartitionsOverride);
  commandLine.RegisterViewWrapper<bool>(keys::overridePartitionNumbers);
  commandLine.RegisterViewWrapper<std::string>(keys::schema);

  commandLine.RegisterViewWrapper<int32>(keys::K);

}

void ProblemManager::ParseCommandLineInput( int & argc, char* argv[])
{
  dataRepository::ManagedGroup& commandLine = GetGroup<ManagedGroup>(keys::commandLine);
  
  ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);
  ViewWrapper<std::string>::rtype  restartFileName = commandLine.getData<std::string>(keys::restartFileName);
  bool&         beginFromRestart = *(commandLine.getData<bool>(keys::beginFromRestart));
  int32&        xPartitionsOverride = *(commandLine.getData<int32>(keys::xPartitionsOverride));
  int32&        yPartitionsOverride = *(commandLine.getData<int32>(keys::yPartitionsOverride));
  int32&        zPartitionsOverride = *(commandLine.getData<int32>(keys::zPartitionsOverride));
  bool&         overridePartitionNumbers = *(commandLine.getData<bool>(keys::overridePartitionNumbers));
  ViewWrapper<std::string>::rtype  schemaName = commandLine.getData<std::string>(keys::schema);

  // Set the options structs and parse
  enum optionIndex {UNKNOWN, HELP, INPUT, RESTART, XPAR, YPAR, ZPAR, SCHEMA};
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
    { 0, 0, 0, 0, 0, 0}
  };

  argc-=(argc>0); 
  argv+=(argc>0);
  option::Stats stats(usage, argc, argv);
  option::Option options[100];//stats.options_max];
  option::Option buffer[100];//stats.buffer_max];
  option::Parser parse(usage, argc, argv, options, buffer);

  
  // Handle special cases
  if (parse.error())
  {
    throw std::invalid_argument("Bad input arguments");
  }

  if (options[HELP] || (argc == 0))
  {
    int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
    option::printUsage(fwrite, stdout, usage, columns);
    exit(0);
  }

  if (options[INPUT] == 0)
  {
    std::cout << "An input xml must be specified!  Exiting..." << std::endl;
    exit(1);
  }


  // Iterate over the remaining inputs
  for (int ii=0; ii<parse.optionsCount(); ++ii)
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
        beginFromRestart = true;
        break;
      case XPAR:
        xPartitionsOverride = std::stoi(opt.arg);
        overridePartitionNumbers = true;
        break;
      case YPAR:
        yPartitionsOverride = std::stoi(opt.arg);
        overridePartitionNumbers = true;
        break;
      case ZPAR:
        zPartitionsOverride = std::stoi(opt.arg);
        overridePartitionNumbers = true;
        break;
      case SCHEMA:
        schemaName = opt.arg;
        break;
    }
  }
}


void ProblemManager::InitializePythonInterpreter()
{
#if USE_PYTHON==1
  // Initialize python and numpy
  std::cout << "Loading python interpreter" << std::endl;

  // Check to make sure the appropriate environment variables are set
  if (getenv("GPAC_SCHEMA") == NULL)
  {
    throw std::invalid_argument("GPAC_SCHEMA must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONPATH") == NULL)
  {
    throw std::invalid_argument("GEOS_PYTHONPATH must be defined to use the new preprocessor!");
  }
  if (getenv("GEOS_PYTHONHOME") == NULL)
  {
    throw std::invalid_argument("GEOS_PYTHONHOME must be defined to use the new preprocessor!");
  }

  setenv("PYTHONPATH", getenv("GEOS_PYTHONPATH"), 1);
  Py_SetPythonHome(getenv("GEOS_PYTHONHOME"));
  Py_Initialize();
  import_array();
#endif
}


void ProblemManager::ClosePythonInterpreter()
{
#if USE_PYTHON==1
  // Add any other cleanup here
  std::cout << "Closing python interpreter" << std::endl;
  Py_Finalize();
#endif
}


void ProblemManager::ParseInputFile()
{
  DomainPartition& domain  = getDomainPartition();

  dataRepository::ManagedGroup& commandLine = GetGroup<ManagedGroup>(keys::commandLine);
  ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);


#if USE_PYTHON==1
  // Load the pygeos module
  PyObject *pModule = PyImport_ImportModule("pygeos");
  if (pModule == NULL)
  {
    PyErr_Print();
    throw std::invalid_argument("Could not find the pygeos module in GEOS_PYTHONPATH!");
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
  std::cout << "Warning: GEOS must be configured to use Python to use parameters, symbolic math, etc. in input files" << std::endl;
#endif


  // Load preprocessed xml file and check for errors
  xmlResult = xmlDocument.load_file(inputFileName.c_str());
  if (!xmlResult)
  {
    std::cout << "XML parsed with errors!" << std::endl;
    std::cout << "Error description: " << xmlResult.description() << std::endl;
    std::cout << "Error offset: " << xmlResult.offset << std::endl;
  }
  xmlProblemNode = xmlDocument.child("Problem");
  pugi::xml_node topLevelNode;

  // Call manager readXML methods:
  this->m_physicsSolverManager->ReadXML(domain, xmlProblemNode );
  this->m_eventManager->ReadXML( xmlProblemNode );
  

  // Documentation output
  ViewWrapper<std::string>::rtype  schemaName = commandLine.getData<std::string>(keys::schema);

  std::cout << schemaName << ", " << schemaName.empty() << ", " << schemaName.size() << std::endl;

  if (schemaName.empty() == 0)
  {
    // m_inputDocumentationHead.Write("test_output.xml");
    ConvertDocumentationToSchema(schemaName.c_str(), *(getDocumentationNode())) ;
    getDocumentationNode()->Print();
  }
}


void ProblemManager::InitializeObjects()
{
  DomainPartition& domain  = getDomainPartition();

  // Initialize solvers
  ViewWrapper<string_array>::rtype  solverNames = this->m_physicsSolverManager->getData<string_array>(keys::solverNames);
  for (auto ii=0; ii<this->m_physicsSolverManager->size(); ++ii)
  {
    SolverBase& currentSolver = this->m_physicsSolverManager->GetGroup<SolverBase>( solverNames[ii] );
    currentSolver.Initialize( domain );
  }
}


void ProblemManager::RunSimulation()
{
  DomainPartition& domain  = getDomainPartition();

  double time = 0.0;
  int cycle = 0;
  real64 dt = 0.0;

  cxx_utilities::DocumentationNode * const eventDocNode = m_eventManager->getDocumentationNode();
  for( auto const & subEventDocNode : eventDocNode->m_child )
  {
    dataRepository::ManagedGroup& currentApplication = m_eventManager->GetGroup( subEventDocNode.first );

    ViewWrapper<string_array>::rtype solverList = currentApplication.getData<string_array>(keys::solvers);
    real64& appDt = *(currentApplication.getData<real64>(keys::dt));
    real64& endTime = *(currentApplication.getData<real64>(keys::endTime));


    bool lockDt = (appDt > 0.0);
    if (lockDt)
    {
      dt = appDt;
    }

    while( time < endTime )
    {
      std::cout << "Time: " << time << "s, dt:" << dt << "s, Cycle: " << cycle << std::endl;
      real64 nextDt = std::numeric_limits<real64>::max();

      for ( auto jj=0; jj<currentApplication.size(); ++jj)
      {
        SolverBase& currentSolver = this->m_physicsSolverManager->GetGroup<SolverBase>( solverList[jj] );
        currentSolver.TimeStep( time, dt, cycle, domain );
        nextDt = std::min(nextDt, *(currentSolver.getData<real64>(keys::maxDt)));
      }

      // Update time, cycle, timestep
      time += dt;
      cycle ++;
      dt = (lockDt)?(dt):(nextDt);
      dt = (endTime - time < dt)?(endTime-time):(dt);
    } 
  }
}




void ProblemManager::ApplySchedulerEvent()
{}


DomainPartition & ProblemManager::getDomainPartition()
{
  return GetGroup<DomainPartition>(keys::domain);
}

DomainPartition const & ProblemManager::getDomainPartition() const
{
  return GetGroup<DomainPartition>(keys::domain);
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, ProblemManager, string const &, ObjectManagerBase * const )

} /* namespace geosx */
