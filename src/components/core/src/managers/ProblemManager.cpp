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
                                SynchronizedGroup * const parent ) :
  SynchronizedGroup( name, parent ),
  m_applicationNames(),
  m_activeSolvers(),
  m_applicationSolvers()
{}

ProblemManager::~ProblemManager()
{}

void ProblemManager::Registration( dataRepository::SynchronizedGroup * const )
{
  RegisterGroup<DomainPartition>(keys::domain);

  SynchronizedGroup& solvers = RegisterGroup<SynchronizedGroup>(keys::solvers);
  solvers.RegisterGroup<SynchronizedGroup>(keys::solverApplications);
  solvers.RegisterViewWrapper<string_array>(keys::solverApplicationNames);

  SynchronizedGroup& commandLine = RegisterGroup<SynchronizedGroup >(keys::commandLine);
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
  dataRepository::SynchronizedGroup& commandLine = GetGroup<SynchronizedGroup>(keys::commandLine);
  
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
      break;

    case 'v':   // version
//      DisplayVersion();
      exit(0);
      break;

    case '?':
      /* getopt_long has already printed an error message. */
      break;

    default:
      abort();
      break;
    }
  }

}


void ProblemManager::InitializePythonInterpreter()
{
#ifdef USE_PYTHON
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
#ifdef USE_PYTHON
  // Add any other cleanup here
  std::cout << "Closing python interpreter" << std::endl;
  Py_Finalize();
#endif
}


void ProblemManager::ParseInputFile()
{
  DomainPartition& domain  = getDomainPartition();

  dataRepository::SynchronizedGroup& commandLine = GetGroup<SynchronizedGroup>(keys::commandLine);
  ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);

  dataRepository::SynchronizedGroup& solvers = GetGroup<dataRepository::SynchronizedGroup>(keys::solvers);
  dataRepository::SynchronizedGroup& solverApplications = solvers.GetGroup<dataRepository::SynchronizedGroup>(keys::solverApplications);
  ViewWrapper<string_array>::rtype  solverApplicationNames = solvers.getData<string_array>(keys::solverApplicationNames);


#ifdef USE_PYTHON
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


  // Solvers
  topLevelNode = xmlProblemNode.child("Solvers");
  std::cout << "Solvers:" << std::endl;
  if (topLevelNode == NULL)
  {
    throw std::invalid_argument("Solver block not present in input xml file!");
  }
  else
  {
    for (pugi::xml_node solverNode=topLevelNode.first_child(); solverNode; solverNode=solverNode.next_sibling())
    {
      std::cout << "   " << solverNode.name() << std::endl;

      // Register the new solver
      std::string solverID = solverNode.attribute("name").value();
      std::unique_ptr<SolverBase> solverPtr = SolverBase::CatalogInterface::Factory(solverNode.name(), solverID, &domain );
      SolverBase& newSolver = solvers.RegisterGroup( solverID, std::move(solverPtr) );
      
      // Register fields in the solver and parse options
      newSolver.Registration( &domain );
      newSolver.ReadXML(solverNode);
      m_activeSolvers.push_back(solverID);
    }
  }


  // Applications
  topLevelNode = xmlProblemNode.child("SolverApplications");
  if (topLevelNode == NULL)
  {
    throw std::invalid_argument("SolverApplications block not present in input xml file!");
  }
  else
  {
    for (pugi::xml_node applicationNode=topLevelNode.first_child(); applicationNode; applicationNode=applicationNode.next_sibling())
    {
      // Register a new solver application (Note: these must be identified by a unique name)
      std::string applicationName = applicationNode.attribute("name").value();
      dataRepository::SynchronizedGroup& newApplication = solverApplications.RegisterGroup<SynchronizedGroup>(applicationName);
      newApplication.RegisterViewWrapper<real64>(keys::beginTime);
      newApplication.RegisterViewWrapper<real64>(keys::endTime);
      newApplication.RegisterViewWrapper<real64>(keys::dt);
      newApplication.RegisterViewWrapper<string_array>(keys::solverList);

      // Read application values from the xml
      *(newApplication.getData<real64>(keys::beginTime)) = applicationNode.attribute("beginTime").as_double(0.0);
      *(newApplication.getData<real64>(keys::endTime)) = applicationNode.attribute("endTime").as_double(0.0);
      *(newApplication.getData<real64>(keys::dt)) = applicationNode.attribute("dt").as_double(-1.0);
      // ViewWrapper<string_array>::rtype solverList = newApplication.getData<string_array>(keys::solverList);
      // applicationNode.attribute("solvers").load_string_array(solverList);
      
      // Having difficulty storing a string_array, so do this for now:
      m_applicationNames.push_back(applicationName);
      std::vector<std::string> newApplicationSolvers;
      applicationNode.attribute("solvers").load_string_array(newApplicationSolvers);
      m_applicationSolvers.insert(applicationSet(applicationName, newApplicationSolvers));
    }

    // Test to make sure the applications are valid
    for (uint ii=0; ii<m_applicationNames.size()-1; ++ii)
    {
      dataRepository::SynchronizedGroup& applicationA = solverApplications.GetGroup(m_applicationNames[ii]);
      dataRepository::SynchronizedGroup& applicationB = solverApplications.GetGroup(m_applicationNames[ii+1]);

      ViewWrapper<real64>::rtype endTime = applicationA.getData<real64>(keys::endTime);
      ViewWrapper<real64>::rtype beginTime = applicationB.getData<real64>(keys::beginTime);
      if (fabs(*(beginTime) - *(endTime)) > 1e-6)
      {
        std::cout << "Error in solver application times: " << m_applicationNames[ii] << std::endl;
        throw std::invalid_argument("Solver application times must be contiguous!");
      }
    }
  }
}


void ProblemManager::InitializeObjects()
{
  DomainPartition& domain  = getDomainPartition();

  // Initialize solvers
  dataRepository::SynchronizedGroup& solvers = GetGroup<dataRepository::SynchronizedGroup>(keys::solvers);
  for (std::vector<string>::iterator ss=m_activeSolvers.begin(); ss!=m_activeSolvers.end(); ++ss)
  {
    SolverBase& currentSolver = solvers.GetGroup<SolverBase>( *ss );
    currentSolver.Initialize( domain );
  }
}


void ProblemManager::RunSimulation()
{
  DomainPartition& domain  = getDomainPartition();
  dataRepository::SynchronizedGroup& solvers = GetGroup<dataRepository::SynchronizedGroup>(keys::solvers);
  dataRepository::SynchronizedGroup& solverApplications = solvers.GetGroup<dataRepository::SynchronizedGroup>(keys::solverApplications);

  double time = 0.0;
  int cycle = 0;
  real64 dt = 0.0;

  for( std::vector<string>::iterator app=m_applicationNames.begin(); app!=m_applicationNames.end(); ++app)
  {
    dataRepository::SynchronizedGroup& currentApplication = solverApplications.GetGroup(*app);
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

      for (std::vector<string>::iterator ss=m_applicationSolvers[*app].begin(); ss!=m_applicationSolvers[*app].end(); ++ss)
      {
        SolverBase& currentSolver = solvers.GetGroup<SolverBase>( *ss );
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

REGISTER_CATALOG_ENTRY( SynchronizedGroup, ProblemManager, std::string const &, SynchronizedGroup * const )

} /* namespace geosx */
