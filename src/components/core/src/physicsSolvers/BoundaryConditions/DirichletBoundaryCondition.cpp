/*
 * DirichletBoundaryCondition.cpp
 *
 *  Created on: Jun 2, 2017
 *      Author: rrsettgast
 */

#include "DirichletBoundaryCondition.hpp"

namespace geosx
{
using namespace dataRepository;

DirichletBoundaryCondition::DirichletBoundaryCondition( string const & name, ManagedGroup *const parent ):
  BoundaryConditionBase(name,parent)
{
  // TODO Auto-generated constructor stub

}

DirichletBoundaryCondition::~DirichletBoundaryCondition()
{
  // TODO Auto-generated destructor stub
}



REGISTER_CATALOG_ENTRY( BoundaryConditionBase, DirichletBoundaryCondition, string const &, ManagedGroup * const )

} /* namespace geosx */
