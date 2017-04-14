/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"
#include "Managers/ElementManager.hpp"
#include "Managers/NodeManager.hpp"

namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ) :
  ManagedGroup( name, parent )
{}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::BuildDataStructure( dataRepository::ManagedGroup * const /*domain*/ )
{
  this->RegisterViewWrapper<real64>(keys::maxDt);
  this->RegisterViewWrapper<real64>(keys::courant);
}

void SolverBase::FillDocumentationNode( dataRepository::ManagedGroup * const  )
{


  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed groups, this could be done automatically
  docNode->setSchemaType("Node");

  docNode->AllocateChildNode( keys::courant,
                              keys::courant,
                              -1,
                              "real64",
                              "real64",
                              "courant Number",
                              "courant Number",
                              "0.7",
                              "",
                              1,
                              1,
                              0 );

  docNode->AllocateChildNode( keys::maxDt,
                              keys::maxDt,
                              -1,
                              "real64",
                              "real64",
                              "Maximum Stable Timestep",
                              "Maximum Stable Timestep",
                              "0.0",
                              "",
                              0,
                              1,
                              0 );



}


void SolverBase::Initialize( dataRepository::ManagedGroup& /*domain*/ )
{
  *(this->getData<real64>(keys::courant)) = std::numeric_limits<real64>::max();
}

} /* namespace ANST */
