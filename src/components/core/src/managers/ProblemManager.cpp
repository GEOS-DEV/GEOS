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

namespace geosx
{

using namespace dataRepository;

ProblemManager::ProblemManager( const std::string& name,
                                ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

ProblemManager::~ProblemManager()
{}

void ProblemManager::Registration( dataRepository::ManagedGroup * const )
{
  RegisterGroup<DomainPartition>(keys::domain);
  RegisterGroup<ManagedGroup>("solvers");

  RegisterViewWrapper< std::unordered_map<string,string> >("simulationParameterMap");
}

void ProblemManager::ParseCommandLineInput( int const& argc, char* const argv[])
{
  std::unordered_map<string,string>& simulationParameterMap = getReference< std::unordered_map<string,string> >(keys::simulationParameterMap);

  // Default values
  std::string fileRootString = "";
//  std::string meshFileString, demeshFileString, edemeshFileString, inputFileString, fpmeshFileString;
  int displaySolvers_flag = 0;
  int displayFields_flag = 0;
  int displayUnits_flag = 0;
  int displayHistory_flag = 0;
  int reportParameter_flag = 0;

  int defaultVariableReportLevel = 0;

  bool setPartitions_flag = false;
  unsigned int xPartitions = 1;
  unsigned int yPartitions = 1;
  unsigned int zPartitions = 1;

  bool m_doWriteXML = false;

  string_array commandLineIncludedFileList;

  // Get command line input
  while (true)
  {
    static struct option long_options[] =
    {
      { "help", no_argument, 0, 'h' },
      { "version", no_argument, 0, 'v' },
      { "solvers", no_argument, &displaySolvers_flag, 1 },
      { "fields", no_argument, &displayFields_flag, 1 },
      { "units", no_argument, &displayUnits_flag, 1 },
      { "history", no_argument, &displayHistory_flag, 1 },
      { "record_defaults", no_argument, &defaultVariableReportLevel, 1 },
      { "report_defaults", no_argument, &defaultVariableReportLevel, 2 },
      { "disable_defaults", no_argument, &defaultVariableReportLevel, 3 },
      { "report_parameters",no_argument,&reportParameter_flag,1},
      { "xpar", required_argument, 0, 0 },
      { "ypar", required_argument, 0, 0 },
      { "zpar", required_argument, 0, 0 },
      { "include", required_argument, 0, 0 },
      { "write_XML", required_argument, 0, 0 },
      { 0, 0, 0, 0 } };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long_only(argc, argv, "ahvf:p:m:d:i:e:r:", long_options, &option_index);

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
        xPartitions = std::stoi(optarg);
        setPartitions_flag = true;
      }
      else if( stringutilities::streq( std::string("ypar"), long_options[option_index].name ) )
      {
        yPartitions = std::stoi(optarg);
        setPartitions_flag = true;
      }
      else if( stringutilities::streq( std::string("zpar"), long_options[option_index].name ) )
      {
        zPartitions = std::stoi(optarg);
        setPartitions_flag = true;
      }
      else if( stringutilities::streq( std::string("include"), long_options[option_index].name ) )
      {
        commandLineIncludedFileList.push_back(optarg);
      }
      else if( stringutilities::streq( std::string("write_XML"), long_options[option_index].name ) )
      {
        m_doWriteXML = true;
        RegisterViewWrapper<string>("xmlOutputFileName").reference() = optarg;
      }

    }
    break;
    case 'a':   // Leave Empty: Included for totalview - does nothing
      break;



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
    }
    break;

    case 'r':   // From restart
    {
//      m_beginFromRestart = true;
      RegisterViewWrapper<int32>("beginFromRestart").reference() = true;
//      m_beginFromRestartFileName = optarg;
      RegisterViewWrapper<string>("beginFromRestartFileName").reference() = optarg;

    }
    break;

    case 'p':   // Record model parameter key=value
    {
      std::string keyValStr = optarg;
      string_array keyVal = stringutilities::Tokenize(keyValStr, "=");
      if (keyVal.size() == 2)
      {
        simulationParameterMap[keyVal[0]] = keyVal[1];
      }
      else
      {
        SLIC_ERROR("Error reading command line input: Parameter " + keyValStr);
      }
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

  if (displaySolvers_flag)
  {
//    DisplaySolvers();
    exit(0);
  }


  if (displayUnits_flag)
  {
//    DisplayUnits();
    exit(0);
  }


  if (displayHistory_flag)
  {
//    DisplayVersionHistory(m_rank);
    exit(0);
  }

//  if (displayFields_flag)
//  {
//    m_displayFields = true; // nb we want to initialize solvers, BC's initial conditions etc before reporting fields.
//    m_displaySplash = false;
//  }

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

  // this option sets the flag to report the current value of parameters
//  if(reportParameter_flag)
//  {
//    m_echoParameters = true;
//  }

//  if(setPartitions_flag)
//    m_partition.setPartitions(xPartitions,yPartitions,zPartitions );
//
//////////////////////////////////////////////////////////////////////////////
}

void ProblemManager::ParseInputFile()
{

  FiniteElementSpace feSpace( keys::FE_Space, this);

  dataRepository::ManagedGroup& solvers = GetGroup<dataRepository::ManagedGroup>(keys::solvers);
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
