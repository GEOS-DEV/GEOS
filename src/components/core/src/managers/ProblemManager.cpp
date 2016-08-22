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

<<<<<<< HEAD
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

using namespace TICPP;
=======
>>>>>>> feature/settgast/addLegacyCode
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

<<<<<<< HEAD
  WrapperCollection& commandLine = RegisterChildWrapperCollection<WrapperCollection >(keys::commandLine);
  commandLine.RegisterWrapper<std::string>(keys::inputFileName);
  // commandLine.RegisterWrapper<std::string>(keys::restartFileName);
  commandLine.RegisterWrapper<bool>(keys::beginFromRestart);
  commandLine.RegisterWrapper<int32>(keys::xPartitionsOverride);
  commandLine.RegisterWrapper<int32>(keys::yPartitionsOverride);
  commandLine.RegisterWrapper<int32>(keys::zPartitionsOverride);
  commandLine.RegisterWrapper<bool>(keys::overridePartitionNumbers);  
=======
  RegisterViewWrapper< std::unordered_map<string,string> >("simulationParameterMap");
>>>>>>> feature/settgast/addLegacyCode
}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{
  dataRepository::WrapperCollection& commandLine = GetChildWrapperCollection<dataRepository::WrapperCollection>(keys::commandLine);
  
  // Wrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);
  // Wrapper<std::string>::rtype  restartFileName = commandLine.getData<std::string>(keys::restartFileName);
  Wrapper<bool>::rtype         beginFromRestart = commandLine.getData<bool>(keys::beginFromRestart);
  Wrapper<int32>::rtype        xPartitionsOverride = commandLine.getData<int32>(keys::xPartitionsOverride);
  Wrapper<int32>::rtype        yPartitionsOverride = commandLine.getData<int32>(keys::yPartitionsOverride);
  Wrapper<int32>::rtype        zPartitionsOverride = commandLine.getData<int32>(keys::zPartitionsOverride);
  Wrapper<bool>::rtype         overridePartitionNumbers = commandLine.getData<bool>(keys::overridePartitionNumbers);


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
<<<<<<< HEAD
=======
      else if( stringutilities::streq( std::string("include"), long_options[option_index].name ) )
      {
        commandLineIncludedFileList.push_back(optarg);
      }
      else if( stringutilities::streq( std::string("write_XML"), long_options[option_index].name ) )
      {
        m_doWriteXML = true;
        RegisterViewWrapper<string>("xmlOutputFileName").reference() = optarg;
      }

>>>>>>> feature/settgast/addLegacyCode
    }
    break;
    case 'a':   // Leave Empty: Included for totalview - does nothing
      break;


<<<<<<< HEAD
    case 'i':   // Record input file
    {
      // inputFileName = optarg;
=======

    //  if (!fileRootString.empty())
    //    m_FileManager.SetRoot(fileRootString.c_str());
    //  if (!inputFileString.empty())
    //    m_FileManager.SetInputFilename(inputFileString.c_str());
    //  if (!meshFileString.empty())
    //    m_FileManager.SetGeometryFilename(meshFileString.c_str());
    //  if (!demeshFileString.empty())
    //    m_FileManager.SetDiscreteElementGeometryFilename(demeshFileString.c_str());
    //  if (!edemeshFileString.empty())
    //    m_FileManager.SetEllipsoidalDiscreteElementGeometryFilename(edemeshFileString.c_str());
    //#ifdef SRC_EXTERNAL
    //  if (!fpmeshFileString.empty())
    //    m_FileManager.SetFaultPatchElementGeometryFilename(fpmeshFileString.c_str());
    //#endif

    case 'f':   // Record file root
    {
      RegisterViewWrapper<string>("fileRootString").reference() = optarg;
    }
    break;

    case 'i':   // Record input file
    {
      RegisterViewWrapper<string>("inputFileString").reference() = optarg;
    }
    break;

    case 'm':   // Record mesh file
    {
      RegisterViewWrapper<string>("meshFileString").reference() = optarg;
    }
    break;
    case 'd':   // Record discrete element mesh file
    {
      RegisterViewWrapper<string>("demeshFileString").reference() = optarg;
    }
    break;
    case 'e':   // Record ellipsoidal discrete element mesh file
    {
      RegisterViewWrapper<string>("edemeshFileString").reference() = optarg;
    }
    break;

    case 's':   // Record seismicity fault patch mesh file
    {
      RegisterViewWrapper<string>("fpmeshFileString").reference() = optarg;
>>>>>>> feature/settgast/addLegacyCode
    }
    break;

    case 'r':   // From restart
    {
<<<<<<< HEAD
      *(beginFromRestart) = true;
      // restartFileName = optarg;
=======
//      m_beginFromRestart = true;
      RegisterViewWrapper<int32>("beginFromRestart").reference() = true;
//      m_beginFromRestartFileName = optarg;
      RegisterViewWrapper<string>("beginFromRestartFileName").reference() = optarg;

>>>>>>> feature/settgast/addLegacyCode
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


<<<<<<< HEAD
void ProblemManager::InitializePythonInterpreter()
{
  // Initialize python and numpy
  std::cout << "Initializing python...";
  Py_Initialize() ;
  import_array();
  std::cout << "  done!" << std::endl;
=======
  // this option concerns flags for the default values of variables
//  switch(defaultVariableReportLevel)
//  {
//  case 0:
//    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::silent);
//    break;
//  case 1:
//    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::recordDefaults);
//    break;
//  case 2:
//    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::reportDefaults);
//    break;
//  case 3:
//    HierarchicalDataNode::SetDefaultReportLevel(HierarchicalDataNode::disableDefaults);
//    break;
//  }
>>>>>>> feature/settgast/addLegacyCode

  // Add a test here to make sure a supported version of python is available
}

<<<<<<< HEAD
void ProblemManager::ClosePythonInterpreter()
{
  // Add any other cleanup here
  Py_Finalize();
=======
//  if(setPartitions_flag)
//    m_partition.setPartitions(xPartitions,yPartitions,zPartitions );
//
//////////////////////////////////////////////////////////////////////////////
>>>>>>> feature/settgast/addLegacyCode
}


void ProblemManager::ParseInputFile()
{
  dataRepository::WrapperCollection& commandLine = GetChildWrapperCollection<dataRepository::WrapperCollection>(keys::commandLine);
  
  // Wrapper<std::string>::rtype  inputFileName = commandLine.getData<std::string>(keys::inputFileName);

<<<<<<< HEAD
  // Hard-wired parse:
  FiniteElementSpace feSpace( keys::FE_Space , this);
=======
  FiniteElementSpace feSpace( keys::FE_Space, this);
>>>>>>> feature/settgast/addLegacyCode

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
