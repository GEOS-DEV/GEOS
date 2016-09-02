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
                        SynchronizedGroup * const parent ) :
  SynchronizedGroup( name, parent )
{}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void SolverBase::ReadXML( pugi::xml_node solverNode )
{
  *(this->getData<real64>(keys::courant)) = solverNode.attribute("courant").as_double(0.5);
}

void SolverBase::Registration( dataRepository::SynchronizedGroup * const /*domain*/ )
{
  this->RegisterViewWrapper<real64>(keys::maxDt);
  this->RegisterViewWrapper<real64>(keys::courant);
}

void SolverBase::Initialize( dataRepository::SynchronizedGroup& /*domain*/ )
{
  *(this->getData<real64>(keys::courant)) = std::numeric_limits<real64>::max();
}

} /* namespace ANST */
