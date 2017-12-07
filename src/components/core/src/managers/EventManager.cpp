/*
 * EventManager.cpp
 *
 *  Created on: Oct 5, 2016
 *      Author: sherman
 */

#include "EventManager.hpp"
#include "Event.hpp"

#include "DocumentationNode.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

EventManager::EventManager( std::string const & name,
                            ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{
  this->RegisterViewWrapper< map<double,std::unique_ptr<Event> > >(keys::Events);
}

EventManager::~EventManager()
{}

void EventManager::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  // Set the name to SolverApplications for now
  docNode->setName("Events");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Contains the set of solver applications");
}


void EventManager::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Event: " << childKey << ", " << childName << std::endl;
  this->RegisterGroup<SolverApplication>( childName );
}


void EventManager::CheckEventTiming()
{
  cxx_utilities::DocumentationNode * const docNode = getDocumentationNode();

  for (std::map<std::string,DocumentationNode>::iterator eit=docNode->m_child.begin() ; eit!=docNode->m_child.end() ; ++eit)
  {
    dataRepository::ManagedGroup * applicationA = GetGroup(eit->first);
    ViewWrapper<real64>::rtype endTime = applicationA->getData<real64>(keys::endTime);

    if (++eit != docNode->m_child.end())
    {
      dataRepository::ManagedGroup * applicationB = GetGroup(eit->first);
      ViewWrapper<real64>::rtype beginTime = applicationB->getData<real64>(keys::beginTime);

      if (fabs(*(beginTime) - *(endTime)) > 1e-6)
      {
        std::cout << "Error in solver application times: " << eit->first << std::endl;
        throw std::invalid_argument("Solver application times must be contiguous!");
      }

      --eit;
    }
  }
}



SolverApplication::SolverApplication( std::string const & name,
                                      ManagedGroup * const parent ):
  ManagedGroup( name, parent)
{}

SolverApplication::~SolverApplication()
{}

void SolverApplication::FillDocumentationNode()
{
  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  // docNode->setName(this->CatalogName());
  docNode->setName("Application");
  docNode->setSchemaType("Node");
  docNode->setShortDescription("Describes the timing of the solver application");

  docNode->AllocateChildNode( keys::beginTime,
                              keys::beginTime,
                              -1,
                              "real64",
                              "real64",
                              "application start time",
                              "application start time",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::endTime,
                              keys::endTime,
                              -1,
                              "real64",
                              "real64",
                              "application endTime",
                              "application endTime",
                              "1.0e9",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::dt,
                              keys::dt,
                              -1,
                              "real64",
                              "real64",
                              "application dt",
                              "application dt",
                              "-1.0",
                              "",
                              0,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::solvers,
                              keys::solvers,
                              -1,
                              "string_array",
                              "string_array",
                              "application solvers",
                              "application solvers",
                              "",
                              "",
                              0,
                              1,
                              0 );

}


} /* namespace geosx */
