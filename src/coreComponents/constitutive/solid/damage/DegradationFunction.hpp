#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_

#include <vector>

namespace geosx
{
namespace constitutive
{

//Declare enumeration
enum degradationTypes { QUADRATIC = "Quadratic", QUASIQUADRATIC = "QuasiQuadratic"};  

//Template class declaration
template< degradationTypes D >
class DegradationFunction
{};
// template< char const *DEGRADATION_TYPE >
// class DegradationFunction
// {};

//Standard quadratic degradation functions
template<>
class DegradationFunction<QUADRATIC>
{
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getValue( real64 const d, real64 const eps, std::vector<real64> const paramValues ) const
  {
//    if( m_extDrivingForceFlag )
//    {
      //pf = fmax( fmin( 1.0, m_damage( k, q )), 0.0 );
//    }
    GEOSX_UNUSED_VAR( paramValues );
    return ((1 - eps)*(1 - d)*(1 - d) + eps);
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getDerivative( real64 const d, real64 const eps, std::vector<real64> const paramValues ) const
  {
    GEOSX_UNUSED_VAR( paramValues );
    return -2.0*(1-eps)*(1 - d);
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getSecondDerivative( real64 const d, real64 const eps, std::vector<real64> const paramValues ) const
  {
    GEOSX_UNUSED_VAR( paramValues );
    GEOSX_UNUSED_VAR( d );
    return 2.0*(1-eps);
  }

  ///accessor for list of material parameters needed for this degradation function  
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static std::vector<string> getParameterList()
  {
    return parameterList;
  }  
  
  private:
  std::vector<string> parameterList = {}; //no parameters needed for quadratic degradation function  

};

template<>
class DegradationFunction<QUASIQUADRATIC>
{
  //Lorentz type degradation (cohesive fracture)
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getValue( real64 const d, real64 const eps, std::vector<real64> const paramValues ) const
  {
    GEOSX_UNUSED_VAR( eps );
    real64 const ell = paramValues[0]; 
    real64 const gc = paramValues[1];
    real64 const psic = paramValues[2];
    real64 m = 3*gc/(8*ell*psic);
    real64 p = 1;
    return pow( 1 - m_damage( k, q ), 2 ) /( pow( 1 - m_damage( k, q ), 2 ) + m * m_damage( k, q ) * (1 + p*m_damage( k, q )) );
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static real64 getDerivative( real64 const d, real64 const eps, std::vector<real64> const paramValues ) const
  {
    GEOSX_UNUSED_VAR( eps );
    real64 const ell = paramValues[0]; 
    real64 const gc = paramValues[1];
    real64 const psic = paramValues[2];
    real64 m = 3*gc/(8*ell*psic);
    real64 p = 1;
    return -m*(1 - d)*(1 + (2*p + 1)*d) / pow( pow( 1-d, 2 ) + m*d*(1+p*d), 2 );
  }


  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  virtual real64 getSecondDerivative( real64 const d, std::vector<real64> const paramValues ) const override
  {
    GEOSX_UNUSED_VAR( eps );
    real64 const ell = paramValues[0]; 
    real64 const gc = paramValues[1];
    real64 const psic = paramValues[2];
    real64 m = 3*gc/(8*ell*psic);
    real64 p = 1;
    return -2*m*( pow( d, 3 )*(2*m*p*p + m*p + 2*p + 1) + pow( d, 2 )*(-3*m*p*p -3*p) + d*(-3*m*p - 3) + (-m+p+2) )/pow( pow( 1-d, 2 ) + m*d*(1+p*d), 3 );
  }

  ///accessor for list of material parameters needed for this degradation function  
  GEOSX_FORCE_INLINE
  GEOSX_HOST_DEVICE
  static std::vector<string> getParameterList()
  {
    return parameterList;
  }  
  
  private:
  std::vector<string> parameterList = {"lengthScale","criticalFractureEnergy","criticalStrainEnergy"}; //parameters used when computing g(d)  

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGE_DEGRADATION_FUNCTION_HPP_ */
