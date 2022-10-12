#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_

#include <vector>

namespace geosx
{
namespace constitutive
{

//Template class declaration
template< char const *INDICATOR_TYPE >
class PressureFunction
{};

//Standard linear indicator function
template<>
class PressureFunction<"Linear">
{
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getValue( real64 const d) const
  {
    return d;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getDerivative( real64 const d ) const
  {
    return 1.0;
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getSecondDerivative( real64 const d ) const
  {
    return 0.0;
  }

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_PRESSURE_FUNCTION_HPP_ */
