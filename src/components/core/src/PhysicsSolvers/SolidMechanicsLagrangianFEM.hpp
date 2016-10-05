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
class ManagedGroup;
}

class SolidMechanics_LagrangianFEM : public SolverBase
{
public:
  SolidMechanics_LagrangianFEM( const std::string& name,
                                ManagedGroup * const parent );


  virtual ~SolidMechanics_LagrangianFEM();

  static string CatalogName() { return "SolidMechanics_LagrangianFEM"; }

  virtual void FillDocumentationNode( dataRepository::ManagedGroup * const group ) override;

  virtual void ReadXML( pugi::xml_node const & solverNode ) override;

  virtual void BuildDataStructure( dataRepository::ManagedGroup * const domain ) override;

  virtual void Initialize( dataRepository::ManagedGroup& domain ) override;

  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         dataRepository::ManagedGroup& domain ) override;

  void TimeStepExplicit( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         dataRepository::ManagedGroup& domain );



private:
  SolidMechanics_LagrangianFEM();

};

namespace Integration
{

template<typename T>
#if CONTAINERARRAY_RETURN_PTR == 1
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ dy,
               T * const __restrict__ y,
               real64 const dx,
               localIndex const length )
#else
void OnePoint( T dydx,
               T dy,
               T y,
               real64 const dx,
               localIndex const length )
#endif
{
  for( auto a=0 ; a<length ; ++a )
  {
    dy[a] = dydx[a] * dx;
    y[a] += dy[a];
  }
}

template<typename T>
#if CONTAINERARRAY_RETURN_PTR == 1
void OnePoint( T const * const __restrict__ dydx,
               T * const __restrict__ y,
               real64 const dx,
               localIndex const length )
#else
void OnePoint( T const &  dydx,
               T & y,
               real64 const dx,
               localIndex const length )
#endif
{
  for( auto a=0 ; a<length ; ++a )
  {
    y[a] += dydx[a] * dx;
  }
}

} // namespace integration

} /* namespace geosx */

#endif /* SOLID_MECHANICS_LAGRANGIAN_FEM_HPP_ */
