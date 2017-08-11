/*
 * DirichletBoundaryCondition.hpp
 *
 *  Created on: Jun 2, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_

#include "BoundaryConditionBase.hpp"

namespace geosx
{

class DirichletBoundaryCondition : public BoundaryConditionBase
{
public:
  DirichletBoundaryCondition( string const & name, dataRepository::ManagedGroup *const parent );
  DirichletBoundaryCondition() = delete;
  virtual ~DirichletBoundaryCondition();

  static string CatalogName() { return "Dirichlet"; }

  template< typename T >
  static void ApplyBounaryConditionDefaultMethod( BoundaryConditionBase const & bc,
                                                  real64 const time,
                                                  array<T> & field );
};


template< typename T >
void DirichletBoundaryCondition::ApplyBounaryConditionDefaultMethod(
                                 BoundaryConditionBase const & bc,
                                 real64 const time,
                                 array<T> & field )
{
  for( auto & fieldValue : field )
  {
    fieldValue = bc.GetValue(time);
  }
}



} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_ */
