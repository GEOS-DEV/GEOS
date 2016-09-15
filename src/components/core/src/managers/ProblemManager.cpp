/*
 * ProblemManager.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: rrsettgast
 */

#include "ProblemManager.hpp"

#include <getopt.h>
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
    std::string const beginTime = "beginTime";
    std::string const endTime = "endTime";
    std::string const dt = "dt";
    std::string const solverList = "solverList";
  }
}

using namespace dataRepository;

ProblemManager::ProblemManager( const std::string& name,
                                ManagedGroup * const parent ) :
  ManagedGroup( name, parent ),
  m_physicsSolverManager(this->RegisterGroup<PhysicsSolverManager>("PhysicsSolverManager" ) )
{}

ProblemManager::~ProblemManager()
{}

void ProblemManager::Registration( dataRepository::ManagedGroup * const )
{
  RegisterGroup<DomainPartition>(keys::domain);

  ManagedGroup& solverApplications = RegisterGroup<ManagedGroup>(keys::solverApplications);
  solverApplications.RegisterViewWrapper<string_array>(keys::solverApplicationNames);

  ManagedGroup& commandLine = RegisterGroup<ManagedGroup >(keys::commandLine);
  commandLine.RegisterViewWrapper<std::string>(keys::inputFileName);
  commandLine.RegisterViewWrapper<std::string>(keys::restartFileName);
  commandLine.RegisterViewWrapper<bool>(keys::beginFromRestart);
  commandLine.RegisterViewWrapper<int32>(keys::xPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::yPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::zPartitionsOverride);
  commandLine.RegisterViewWrapper<bool>(keys::overridePartitionNumbers);
}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{
  dataRepository::ManagedGroup& commandLine = GetGroup<ManagedGroup>(keys::commandLine);
  
  ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);
  ViewWrapper<std::string>::rtype  restartFileName = commandLine.getData<std::string>(keys::restartFileName);
  bool&         beginFromRestart = *(commandLine.getData<bool>(keys::beginFromRestart));
  int32&        xPartitionsOverride = *(commandLine.getData<int32>(keys::xPartitionsOverride));
  int32&        yPartitionsOverride = *(commandLine.getData<int32>(keys::yPartitionsOverride));
  int32&        zPartitionsOverride = *(commandLine.getData<int32>(keys::zPartitionsOverride));
  bool&         overridePartitionNumbers = *(commandLine.getData<bool>(keys::overridePartitionNumbers));


  // Get command line input
  while (true)
  {
    static struct option long_options[] =
    {
      { "help", no_argument, 0, 'h' },
      { "version", no_argument, 0, 'v' },
      { "xpar", required_argument, 0, 0 },
      { "ypar", required_argument, 0, 0 },
      { "zpar", required_argument, 0, 0 },
      { "include", required_argument, 0, 0 },
      { 0, 0, 0, 0 } };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    int c = getopt_long_only(argc, argv, "ahvi:r:", long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
    {
      /* If option sets a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;

      /* long options without a short arg */
      if( stringutilities::streq( std::string("xpar"), long_options[option_index].name ) )
      {
        xPartitionsOverride = std::stoi(optarg);
        overridePartitionNumbers = true;
      }
      else if( stringutilities::streq( std::string("ypar"), long_options[option_index].name ) )
      {
        yPartitionsOverride = std::stoi(optarg);
        overridePartitionNumbers = true;
      }
      else if( stringutilities::streq( std::string("zpar"), long_options[option_index].name ) )
      {
        zPartitionsOverride = std::stoi(optarg);
        overridePartitionNumbers = true;
      }
    }
    break;
    case 'a':   // Leave Empty: Included for totalview - does nothing
      break;

    case 'i':   // Record input file
    {
      inputFileName = optarg;
    }
    break;

    case 'r':   // From restart
    {
      beginFromRestart = true;
      restartFileName = optarg;
    }
    break;
    
    case 'h':   // help
//      DisplayUsage();   // print help
      exit(0);

    case 'v':   // version
//      DisplayVersion();
      exit(0);

    case '?':
      /* getopt_long has already printed an error message. */
      break;

    default:
      abort();
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

//  dataRepository::SynchronizedGroup& solvers = GetGroup<dataRepository::ManagedGroup>(keys::solvers);
  dataRepository::ManagedGroup& solverApplications = GetGroup<dataRepository::ManagedGroup>(keys::solverApplications);


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

  m_inputDocumentationHead.m_varName = "Problem";
  m_inputDocumentationHead.m_varType = "UniqueNode";
  m_inputDocumentationHead.m_varDescription = "This is the top level node in the input structure.";

  cxx_utilities::InputDocumentation temp;
  temp.m_level   = 1;
  temp.m_varName = "SolverNode";
  temp.m_varType = "";
  temp.m_varDescription = "";

  m_inputDocumentationHead.m_child.insert( { "SolverNode", temp } );

  this->m_physicsSolverManager.ReadXML(domain, xmlProblemNode, m_inputDocumentationHead.m_child["SolverNode"] );


  // Applications
  topLevelNode = xmlProblemNode.child("SolverApplications");
  if (topLevelNode == NULL)
  {
    throw std::invalid_argument("SolverApplications block not present in input xml file!");
  }
  else
  {
    // Determine the number of solver applications, resize the application collection 
    long nApp = std::distance(topLevelNode.children().begin(), topLevelNode.children().end());
    solverApplications.resize(nApp);
    ViewWrapper<string_array>::rtype  solverApplicationNames = solverApplications.getData<string_array>(keys::solverApplicationNames);
    int ii = 0;

    for (pugi::xml_node applicationNode=topLevelNode.first_child(); applicationNode; applicationNode=applicationNode.next_sibling())
    {
      // Register a new solver application (Note: these must be identified by a unique name)
      std::string applicationName = applicationNode.attribute("name").value();
      dataRepository::ManagedGroup& newApplication = solverApplications.RegisterGroup<ManagedGroup>(applicationName);
      newApplication.RegisterViewWrapper<real64>(keys::beginTime);
      newApplication.RegisterViewWrapper<real64>(keys::endTime);
      newApplication.RegisterViewWrapper<real64>(keys::dt);
      newApplication.RegisterViewWrapper<string_array>(keys::solverList);

      // Read application values from the xml
      *(newApplication.getData<real64>(keys::beginTime)) = applicationNode.attribute("beginTime").as_double(0.0);
      *(newApplication.getData<real64>(keys::endTime)) = applicationNode.attribute("endTime").as_double(0.0);
      *(newApplication.getData<real64>(keys::dt)) = applicationNode.attribute("dt").as_double(-1.0);
      
      // Store the application name
      solverApplicationNames[ii] = applicationName;
      ii++;

      // Store the solver list in this application
      std::vector<std::string> newApplicationSolvers;
      applicationNode.attribute("solvers").load_string_array(newApplicationSolvers);
      newApplication.resize(newApplicationSolvers.size());
      ViewWrapper<string_array>::rtype solverList = newApplication.getData<string_array>(keys::solverList);
      for (uint jj=0; jj<newApplicationSolvers.size(); ++jj)
      {
        solverList[jj] = newApplicationSolvers[jj];
      }
    }

    // Test to make sure the applications are valid
    for (auto jj=0; jj<solverApplications.size()-1; ++jj)
    {
      dataRepository::ManagedGroup& applicationA = solverApplications.GetGroup(solverApplicationNames[jj]);
      dataRepository::ManagedGroup& applicationB = solverApplications.GetGroup(solverApplicationNames[jj+1]);

      ViewWrapper<real64>::rtype endTime = applicationA.getData<real64>(keys::endTime);
      ViewWrapper<real64>::rtype beginTime = applicationB.getData<real64>(keys::beginTime);
      if (fabs(*(beginTime) - *(endTime)) > 1e-6)
      {
        std::cout << "Error in solver application times: " << solverApplicationNames[jj] << std::endl;
        throw std::invalid_argument("Solver application times must be contiguous!");
      }
    }
  }

  // m_inputDocumentationHead.Write("test_output.xml");
  ConvertDocumentationToSchema("test_output.xsd", m_inputDocumentationHead);
}


void ProblemManager::InitializeObjects()
{
  DomainPartition& domain  = getDomainPartition();

  // Initialize solvers
  ViewWrapper<string_array>::rtype  solverNames = this->m_physicsSolverManager.getData<string_array>(keys::solverNames);
  for (auto ii=0; ii<this->m_physicsSolverManager.size(); ++ii)
  {
    SolverBase& currentSolver = this->m_physicsSolverManager.GetGroup<SolverBase>( solverNames[ii] );
    currentSolver.Initialize( domain );
  }
}


void ProblemManager::RunSimulation()
{
  DomainPartition& domain  = getDomainPartition();
  dataRepository::ManagedGroup& solverApplications = GetGroup<dataRepository::ManagedGroup>(keys::solverApplications);
  ViewWrapper<string_array>::rtype  solverApplicationNames = solverApplications.getData<string_array>(keys::solverApplicationNames);

  double time = 0.0;
  int cycle = 0;
  real64 dt = 0.0;

  for( auto ii=0; ii<solverApplications.size(); ++ii)
  {
    dataRepository::ManagedGroup& currentApplication = solverApplications.GetGroup( solverApplicationNames[ii] );
    ViewWrapper<string_array>::rtype solverList = currentApplication.getData<string_array>(keys::solverList);
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
        SolverBase& currentSolver = this->m_physicsSolverManager.GetGroup<SolverBase>( solverList[jj] );
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

REGISTER_CATALOG_ENTRY( ManagedGroup, ProblemManager, std::string const &, ManagedGroup * const )

} /* namespace geosx */
