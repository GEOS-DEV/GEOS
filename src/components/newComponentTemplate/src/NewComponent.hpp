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
namespace dataRepository
{
class ManagedGroup;
}

class NewComponent : public SolverBase
{
public:
  NewComponent( std::string const & name,
                ManagedGroup * const parent);
  virtual ~NewComponent();

  static std::string CatalogName() { return "NewComponent"; }

  virtual void ReadXML_PostProcess( ) override;

  virtual void InitializePreSubGroups( ManagedGroup * const group ) override;

  virtual void InitializePostSubGroups( ManagedGroup * const group ) override;


  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         dataRepository::ManagedGroup * domain ) override;

private:
  NewComponent() = delete;
  NewComponent(const NewComponent&) = delete;
  NewComponent(const NewComponent&&) = delete;
  NewComponent& operator=(const NewComponent&) = delete;
  NewComponent& operator=(const NewComponent&&) = delete;
};

} /* namespace geosx */

#endif /* COMPONENTS_NEWCOMPONENTTEMPLATE_SRC_NEWCOMPONENT_HPP_ */
