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
                        WrapperCollection * const parent ) :
  WrapperCollection( name, parent )
{}

SolverBase::~SolverBase()
{}

SolverBase::CatalogInterface::CatalogType& SolverBase::GetCatalogue()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

} /* namespace ANST */
