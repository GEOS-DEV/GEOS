/*
 * NewtonianMechanics.hpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#ifndef SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_
#define SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_

#include "SolverBase.hpp"


namespace geosx
{
namespace dataRepository
{
class SynchronizedGroup;
}

class SolidMechanics_LagrangianFEM : public SolverBase
{
public:
  SolidMechanics_LagrangianFEM( const std::string& name,
                                SynchronizedGroup * const parent );


  virtual ~SolidMechanics_LagrangianFEM();

  static std::string CatalogName() { return "SolidMechanics_LagrangianFEM"; }

  void TimeStepExplicit( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         dataRepository::SynchronizedGroup& domain );

  virtual void ReadXML( pugi::xml_node solverNode );

  virtual void Registration( dataRepository::SynchronizedGroup * const domain );

  virtual void Initialize( dataRepository::SynchronizedGroup& domain );

  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         dataRepository::SynchronizedGroup& domain ) override;


private:
  SolidMechanics_LagrangianFEM();

};

namespace Integration
{

template<typename T>
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ dy,
               T * const __restrict__ y,
               real64 const dx,
               uint64 const length )
{
  for( uint64 a=0 ; a<length ; ++a )
  {
    dy[a] = dydx[a] * dx;
    y[a] += dy[a];
  }
}

template<typename T>
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ y,
               real64 const dx,
               uint64 const length )
{
  for( uint64 a=0 ; a<length ; ++a )
  {
    y[a] += dydx[a] * dx;
  }
}

} // namespace integration

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
