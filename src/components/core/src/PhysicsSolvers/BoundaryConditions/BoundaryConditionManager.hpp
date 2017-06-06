/*
 * BoundaryConditionManager.hpp
 *
 *  Created on: May 26, 2017
 *      Author: rrsettgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_
#define SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_

namespace geosx
{

class BoundaryConditionManager : public dataRepository::ManagedGroup
{
public:
  BoundaryConditionManager();
  virtual ~BoundaryConditionManager();

  void ReadXML( pugi::xml_node const & problemNode );


  template< typename T, typename Solver, typename BCFunctionPtr, typename ... ARGS>
  void ApplyBoundaryCondition( Solver * const solverPtr,
                               BCFunctionPtr boundaryConditionFunctionPtr,
                               string const & bcName,
                               realT const time,
                               realT const dt );
};


template< typename T, typename Solver, typename BCFunctionPtr, typename ... ARGS>
void BoundaryConditionManager::ApplyBoundaryCondition( Solver* solverPtr,
                                                       BCFunctionPtr boundaryConditionFunctionPtr,
                                                       string const & bcName,
                                                       realT const time,
                                                       realT const dt )
{


}



} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_BOUNDARYCONDITIONS_BOUNDARYCONDITIONMANAGER_HPP_ */
