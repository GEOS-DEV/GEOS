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
  }
}

using namespace dataRepository;

ProblemManager::ProblemManager( const std::string& name,
                                SynchronizedGroup * const parent ) :
  SynchronizedGroup( name, parent )
{}

ProblemManager::~ProblemManager()
{}

void ProblemManager::Registration( dataRepository::SynchronizedGroup * const )
{
  RegisterGroup<DomainPartition>(keys::domain);
  RegisterGroup<SynchronizedGroup>("solvers");

  SynchronizedGroup& commandLine = RegisterGroup<SynchronizedGroup >(keys::commandLine);
  commandLine.RegisterViewWrapper<std::string>(keys::inputFileName);
  // commandLine.RegisterViewWrapper<std::string>(keys::restartFileName);
  commandLine.RegisterViewWrapper<bool>(keys::beginFromRestart);
  commandLine.RegisterViewWrapper<int32>(keys::xPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::yPartitionsOverride);
  commandLine.RegisterViewWrapper<int32>(keys::zPartitionsOverride);
  commandLine.RegisterViewWrapper<bool>(keys::overridePartitionNumbers);
}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{
  dataRepository::SynchronizedGroup& commandLine = GetGroup<SynchronizedGroup>(keys::commandLine);
  
  // ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);
  // ViewWrapper<std::string>::rtype  restartFileName = commandLine.getData<std::string>(keys::restartFileName);
  ViewWrapper<bool>::rtype         beginFromRestart = commandLine.getData<bool>(keys::beginFromRestart);
  ViewWrapper<int32>::rtype        xPartitionsOverride = commandLine.getData<int32>(keys::xPartitionsOverride);
  ViewWrapper<int32>::rtype        yPartitionsOverride = commandLine.getData<int32>(keys::yPartitionsOverride);
  ViewWrapper<int32>::rtype        zPartitionsOverride = commandLine.getData<int32>(keys::zPartitionsOverride);
  ViewWrapper<bool>::rtype         overridePartitionNumbers = commandLine.getData<bool>(keys::overridePartitionNumbers);


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
        *(xPartitionsOverride) = std::stoi(optarg);
        *(overridePartitionNumbers) = true;
      }
      else if( stringutilities::streq( std::string("ypar"), long_options[option_index].name ) )
      {
        *(yPartitionsOverride) = std::stoi(optarg);
        *(overridePartitionNumbers) = true;
      }
      else if( stringutilities::streq( std::string("zpar"), long_options[option_index].name ) )
      {
        *(zPartitionsOverride) = std::stoi(optarg);
        *(overridePartitionNumbers) = true;
      }
    }
    break;
    case 'a':   // Leave Empty: Included for totalview - does nothing
      break;

    case 'i':   // Record input file
    {
      // inputFileName = optarg;
    }
    break;

    case 'r':   // From restart
    {
      *(beginFromRestart) = true;
      // restartFileName = optarg;
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
  // Initialize python and numpy
  std::cout << "Initializing python...";
  Py_Initialize() ;
  import_array();
  std::cout << "  done!" << std::endl;
  // Add a test here to make sure a supported version of python is available
}


void ProblemManager::ClosePythonInterpreter()
{
  // Add any other cleanup here
  Py_Finalize();
}


void ProblemManager::ParseInputFile()
{
  dataRepository::SynchronizedGroup& commandLine = GetGroup<SynchronizedGroup>(keys::commandLine);
  
  // ViewWrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);

  FiniteElementSpace feSpace( keys::FE_Space, this);

  dataRepository::SynchronizedGroup& solvers = GetGroup<dataRepository::SynchronizedGroup>(keys::solvers);
  DomainPartition& domain  = getDomainPartition();

  std::string newName("new solver");
  std::string newName2("new solver2");

  std::unique_ptr<SolverBase> solver = SolverBase::CatalogInterface::Factory("SolidMechanics_LagrangianFEM", newName, &domain );
  auto& solver1 = solvers.RegisterGroup( newName, std::move(solver) );
  solver = SolverBase::CatalogInterface::Factory( "NewComponent", newName2, &domain );
  auto& solver2 = solvers.RegisterGroup( newName2, std::move(solver) );

  solver1.getData<string>(std::string("name"));

  solver1.Registration( &domain );
  solver2.Registration( &domain );

  double time = 0.0;
  double dt = 5.0e-5;
  for( int i=0 ; i<10 ; ++i )
  {
    solver1.TimeStep( time, dt, i, domain );
    time += dt;
    std::cout<<std::endl;
  }


  // Preprocess the xml file using python
  InitializePythonInterpreter();
  PyObject *pModule = PyImport_ImportModule("pygeos");
  if (pModule == NULL)
  {
    PyErr_Print();
    throw std::invalid_argument("Could not find the pygeos module in PYTHONPATH!");
  }
  PyObject *pPreprocessorFunction = PyObject_GetAttrString(pModule, "PreprocessGEOSXML");
  PyObject *pPreprocessorInputStr = Py_BuildValue("(s)", "input.xml");
  PyObject *pPreprocessorResult = PyObject_CallObject(pPreprocessorFunction, pPreprocessorInputStr);
  std::string processedInputFile = PyString_AsString(pPreprocessorResult);


  // Load preprocessed xml file and check for errors
  xmlResult = xmlDocument.load_file(processedInputFile.c_str());
  if (!xmlResult)
  {
    std::cout << "XML parsed with errors!" << std::endl;
    std::cout << "Error description: " << xmlResult.description() << std::endl;
    std::cout << "Error offset: " << xmlResult.offset << std::endl;
  }
  xmlProblemNode = xmlDocument.child("Problem");
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
