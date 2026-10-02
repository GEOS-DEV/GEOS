/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYDROSTATICCOORDINATE_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_HYDROSTATICCOORDINATE_HPP_

#include "common/DataTypes.hpp"
#include <cmath>
#include <limits>

namespace geos
{
/** A coordinate basis for hydrostatic initialization, never a change of physical gravity. */
struct HydrostaticCoordinate
{
  real64 up[3] = { 0.0, 0.0, 1.0 };
  real64 gravityComponent = 0.0;
  bool gravityAligned = false;

  HydrostaticCoordinate( real64 const (&gravity)[3], bool const aligned ):
    gravityComponent( gravity[2] ), gravityAligned( aligned )
  {
    if( aligned )
    {
      real64 const magnitude = std::hypot( gravity[0], gravity[1], gravity[2] );
      gravityComponent = -magnitude;
      if( magnitude > 0.0 )
      {
        for( int i = 0; i < 3; ++i ) up[i] = -gravity[i] / magnitude;
      }
    }
  }

  /** Half-open partition, with an explicit finite-precision tie band. */
  GEOS_HOST_DEVICE
  bool isBelow( real64 const value, real64 const contact, real64 const projectionError = 0.0 ) const
  {
    if( !gravityAligned ) return value < contact;
    real64 const a = value < 0.0 ? -value : value;
    real64 const b = contact < 0.0 ? -contact : contact;
    real64 const scalarScale = a > b ? a : b;
    real64 const tolerance = projectionError + 2.0 * std::numeric_limits< real64 >::epsilon() * scalarScale;
    return value < contact - tolerance;
  }

  /** Bound the represented point/normal projection's rounding scale. */
  template< typename POINT >
  GEOS_HOST_DEVICE
  real64 projectionErrorBound( POINT const & point ) const
  {
    real64 scale = 0.0;
    for( int i = 0; i < 3; ++i )
    {
      real64 const product = up[i] * point[i];
      scale += product < 0.0 ? -product : product;
    }
    return 8.0 * std::numeric_limits< real64 >::epsilon() * scale;
  }

  template< typename POINT >
  GEOS_HOST_DEVICE
  real64 elevation( POINT const & point ) const
  {
    if( !gravityAligned ) return up[0] * point[0] + up[1] * point[1] + up[2] * point[2];
    // Error-free product residuals and TwoSum avoid losing small potential
    // distances when large world-coordinate terms nearly cancel. This does
    // not invent precision missing from the represented geometry itself.
    real64 product[3], productError[3];
    for( int i = 0; i < 3; ++i )
    {
      product[i] = up[i] * point[i];
      productError[i] = ::fma( up[i], point[i], -product[i] );
    }
    real64 const first = product[0] + product[1];
    real64 const firstVirtual = first - product[0];
    real64 const firstError = ( product[0] - ( first - firstVirtual ) ) + ( product[1] - firstVirtual );
    real64 const total = first + product[2];
    real64 const secondVirtual = total - first;
    real64 const secondError = ( first - ( total - secondVirtual ) ) + ( product[2] - secondVirtual );
    return total + ( firstError + secondError + productError[0] + productError[1] + productError[2] );
  }

};
/** Collision-free internal table key for distinct mesh/region/subregion owners. */
inline string scopedHydrostaticTableName( string const & prefix, string const & path )
{
  char const digits[] = "0123456789abcdef";
  string name = prefix + "__";
  for( unsigned char c : path )
  {
    name += digits[c >> 4]; name += digits[c & 15];
  }
  return name;
}
}
#endif
