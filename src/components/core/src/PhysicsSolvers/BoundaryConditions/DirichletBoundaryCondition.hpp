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

  static string CatalogName() { return "DirichletBoundaryCondition"; }

};

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_BOUNDARYCONDITIONS_DIRICHLETBOUNDARYCONDITION_HPP_ */
