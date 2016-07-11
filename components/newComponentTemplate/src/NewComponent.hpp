/*
 * NewComponent.hpp
 *
 *  Created on: Jun 8, 2016
 *      Author: settgast
 */

#ifndef COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#define COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_
#include "PhysicsSolvers/SolverBase.hpp"


namespace geosx
{

class WrapperCollection;


class NewComponent : public SolverBase
{
public:
  NewComponent( std::string const & name);
  virtual ~NewComponent();

  static std::string CatalogueName() { return "NewComponent"; }

  virtual void RegisterDataObjects( WrapperCollection& domain ) ;


  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         WrapperCollection& domain ) ;

private:
  NewComponent() = delete;
  NewComponent(const NewComponent&) = delete;
  NewComponent(const NewComponent&&) = delete;
  NewComponent& operator=(const NewComponent&) = delete;
  NewComponent& operator=(const NewComponent&&) = delete;
};

} /* namespace geosx */

#endif /* COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_ */
