/*
 * SolverBase.cpp
 *
 *  Created on: Dec 2, 2014
 *      Author: rrsettgast
 */

#include "SolverBase.hpp"

namespace geosx
{

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
  double courant = solverNode.attribute("courant").as_double(0.5);
}

} /* namespace ANST */
