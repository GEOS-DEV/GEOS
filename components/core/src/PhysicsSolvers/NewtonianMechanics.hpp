/*
 * NewtonianMechanics.hpp
 *
 *  Created on: Dec 4, 2014
 *      Author: rrsettgast
 */

#ifndef NEWTONIANMECHANICS_HPP_
#define NEWTONIANMECHANICS_HPP_

#include "SolverBase.hpp"


namespace geosx
{
class DataObjectManager;

class NewtonianMechanics : public SolverBase
{
public:
  NewtonianMechanics( const std::string& name );


  virtual ~NewtonianMechanics();

  static std::string CatalogueName() { return "NewtonianMechanics"; }

  virtual void RegisterDataObjects( DataObjectManager& domain ) override;


  virtual void TimeStep( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         DataObjectManager& domain ) override;

  void TimeStepExplicit( real64 const& time_n,
                         real64 const& dt,
                         int32 const cycleNumber,
                         DataObjectManager& domain );





private:
  NewtonianMechanics();

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

}

} /* namespace ANST */

#endif /* NEWTONIANMECHANICS_HPP_ */
