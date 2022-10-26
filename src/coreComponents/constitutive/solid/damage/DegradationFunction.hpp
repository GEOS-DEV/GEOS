#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_

#include <vector>

namespace geosx
{
namespace constitutive
{

class DegradationFunction
{
  public:
  //Lorentz type degradation (cohesive fracture)
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getValue( real64 const d, real64 const eps, real64 const ell, real64 const gc, real64 const psic, localIndex const degOption ) 
  {
    real64 gd;
    if(degOption == 0)//quadratic degradation
    {
      gd = pow( 1-d, 2 );
    }
    if(degOption == 1)//quasi-quadratic (cohesive) degradation
    {
      //these operation could be moved to damage.hpp, and that would reduce the number of params passed here, however, that would make
      //us do these extra operations even in the quadratic degradation case.
      real64 const m = 3*gc/(8*ell*psic);
      real64 const p = 1;
      gd = pow( 1-d, 2 ) /( pow( 1-d, 2 ) + m * d * (1 + p*d) );   
    }
    return (1-eps)*gd + eps;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getDerivative( real64 const d, real64 const eps, real64 const ell, real64 const gc, real64 const psic, localIndex const degOption ) 
  {
    real64 gdprime;
    if(degOption == 0)//quadratic degradation
    {
      gdprime = -2.0*(1 - d);
    }
    if(degOption == 1)//quasi-quadratic (cohesive) degradation
    {
      //these operation could be moved to damage.hpp, and that would reduce the number of params passed here, however, that would make
      //us do these extra operations even in the quadratic degradation case.
      real64 const m = 3*gc/(8*ell*psic);
      real64 const p = 1;
      gdprime = -m*(1 - d)*(1 + (2*p + 1)*d) / pow( pow( 1-d, 2 ) + m*d*(1+p*d), 2 );
    }
    return (1-eps)*gdprime;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getSecondDerivative( real64 const d, real64 const eps, real64 const ell, real64 const gc, real64 const psic, localIndex const degOption ) 
  {
    real64 gdprimeprime;
    if(degOption == 0)//quadratic degradation
    {
      gdprimeprime = 2.0;
    }
    if(degOption == 1)//quasi-quadratic (cohesive) degradation
    {
      //these operation could be moved to damage.hpp, and that would reduce the number of params passed here, however, that would make
      //us do these extra operations even in the quadratic degradation case.
      real64 const m = 3*gc/(8*ell*psic);
      real64 const p = 1;
      gdprimeprime = -2*m*( pow( d, 3 )*(2*m*p*p + m*p + 2*p + 1) + pow( d, 2 )*(-3*m*p*p -3*p) + d*(-3*m*p - 3) + (-m+p+2) )/pow( pow( 1-d, 2 ) + m*d*(1+p*d), 3 );
    }
    return (1-eps)*gdprimeprime;
  }

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_ */
