/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"


namespace geosx
{

using namespace dataRepository;

SolverBase::SolverBase( std::string const & name,
                        ManagedGroup * const parent ):
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
{}

void SolverBase::FillDocumentationNode()
{


  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();
  docNode->setName(this->CatalogName());    // If this method lived in Managed
                                            // groups, this could be done
                                            // automatically
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



void SolverBase::TimeStep( real64 const& time_n,
                           real64 const& dt,
                           const int cycleNumber,
                           ManagedGroup * domain )
{
}



void SolverBase::CreateChild( string const & childKey, string const & childName )
{
  std::cout << "Adding Solver: " << childKey << ", " << childName << std::endl;
  this->RegisterGroup( childName, CatalogInterface::Factory( childKey, childName, this ) );
}


//void SolverBase::Initialize( dataRepository::ManagedGroup& /*domain*/ )
//{
//  *(this->getData<real64>(keys::courant)) =
// std::numeric_limits<real64>::max();
//}

} /* namespace ANST */
