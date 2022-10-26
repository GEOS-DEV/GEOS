#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_

#include <vector>

namespace geosx
{
namespace constitutive
{

class PressureFunction
{
  public:
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getValue( real64 const d, localIndex const indOption) 
  {
    real64 md;
    if(indOption == 0)//linear indicator
    {
      md = d;
    }
    if(indOption == 1)//cosine indicator
    {
      md = cos(M_PI*d);
    }
    return md;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getDerivative( real64 const d, localIndex const indOption ) 
  {
    real64 mdprime;
    if(indOption == 0)//linear indicator
    {
      mdprime = 1.0;
    }
    if(indOption == 1)//cosine indicator
    {
      mdprime = sin(M_PI*d);
    }
    return mdprime;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getSecondDerivative( real64 const d, localIndex const indOption ) 
  {
    real64 mdprimeprime;
    if(indOption == 0)//linear indicator
    {
      mdprimeprime = 0.0;
    }
    if(indOption == 1)//cosine indicator
    {
      mdprimeprime = -cos(M_PI*d);
    }
    return mdprimeprime;
  }

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_ */
