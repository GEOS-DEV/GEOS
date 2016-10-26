/*
 * EventManager.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: sherman
 */

#include "EventManager.hpp"

#include "DocumentationNode.hpp"
#include "pugixml/src/pugixml.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

EventManager::EventManager( std::string const & name,
                            ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{
}

EventManager::~EventManager()
{
}



void EventManager::CreateSolverApplication(pugi::xml_node const & applicationNode)
{
  // Register a new solver application (Note: these must be identified by a unique name)
  std::string applicationName = applicationNode.attribute("name").value();
  dataRepository::ManagedGroup& newApplication = RegisterGroup<ManagedGroup>(applicationName);

  cxx_utilities::DocumentationNode * const docNode = getDocumentationNode();
  docNode->setSchemaType("Node");
  docNode->setName("Application");

  cxx_utilities::DocumentationNode * const appDocNode = newApplication.getDocumentationNode();

  appDocNode->AllocateChildNode( keys::beginTime,
                                 keys::beginTime,
                                 -1,
                                 "real64",
                                 "double",
                                 "application start time",
                                 "application start time",
                                 "0.0",
                                 "",
                                 1,
                                 0 );

  appDocNode->AllocateChildNode( keys::endTime,
                                 keys::endTime,
                                 -1,
                                 "real64",
                                 "double",
                                 "application endTime",
                                 "application endTime",
                                 "1.0e9",
                                 "",
                                 1,
                                 0 );

  appDocNode->AllocateChildNode( keys::dt,
                                 keys::dt,
                                 -1,
                                 "real64",
                                 "double",
                                 "application dt",
                                 "application dt",
                                 "-1.0",
                                 "",
                                 1,
                                 0 );

  appDocNode->AllocateChildNode( keys::solvers,
                                 keys::solvers,
                                 -1,
                                 "string",
                                 "string",
                                 "application solvers",
                                 "application solvers",
                                 "",
                                 "",
                                 1,
                                 0 );


  
  // This should be done automatically...
  newApplication.RegisterViewWrapper<real64>(keys::beginTime);
  newApplication.RegisterViewWrapper<real64>(keys::endTime);
  newApplication.RegisterViewWrapper<real64>(keys::dt);
  newApplication.RegisterViewWrapper<string_array>(keys::solvers);

  // Read application values from the xml
  *(newApplication.getData<real64>(keys::beginTime)) = applicationNode.attribute("beginTime").as_double(0.0);
  *(newApplication.getData<real64>(keys::endTime)) = applicationNode.attribute("endTime").as_double(0.0);
  *(newApplication.getData<real64>(keys::dt)) = applicationNode.attribute("dt").as_double(-1.0);

  // Store the solver list in this application
  std::vector<std::string> newApplicationSolvers;
  applicationNode.attribute("solvers").load_string_array(newApplicationSolvers, "");
  newApplication.resize(newApplicationSolvers.size());
  ViewWrapper<string_array>::rtype solvers = newApplication.getData<string_array>(keys::solvers);
  for (uint jj=0; jj<newApplicationSolvers.size(); ++jj)
  {
    solvers[static_cast<int>(jj)] = newApplicationSolvers[jj];
  }
}



void EventManager::ReadXML( pugi::xml_node const & problemNode )
{
  pugi::xml_node topLevelNode = problemNode.child("SolverApplications");
  if (topLevelNode == NULL)
  {
    throw std::invalid_argument("SolverApplications block not present in input xml file!");
  }
  else
  {
    // Allow other event types here?
    for (pugi::xml_node applicationNode=topLevelNode.first_child(); applicationNode; applicationNode=applicationNode.next_sibling())
    {
      CreateSolverApplication(applicationNode);
    }
  }
}


void EventManager::CheckEventTiming()
{
  cxx_utilities::DocumentationNode * const docNode = getDocumentationNode();

  for (std::map<std::string,DocumentationNode>::iterator eit=docNode->m_child.begin(); eit!=docNode->m_child.end(); ++eit)
  {
    dataRepository::ManagedGroup& applicationA = GetGroup(eit->first);
    ViewWrapper<real64>::rtype endTime = applicationA.getData<real64>(keys::endTime);

    if (++eit != docNode->m_child.end())
    {
      dataRepository::ManagedGroup& applicationB = GetGroup(eit->first);
      ViewWrapper<real64>::rtype beginTime = applicationB.getData<real64>(keys::beginTime);

      if (fabs(*(beginTime) - *(endTime)) > 1e-6)
      {
        std::cout << "Error in solver application times: " << eit->first << std::endl;
        throw std::invalid_argument("Solver application times must be contiguous!");
      }

      --eit;
    }
  }
}



} /* namespace geosx */
