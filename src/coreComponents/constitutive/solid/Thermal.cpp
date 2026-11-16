/* SPDX-License-Identifier: LGPL-2.1-only */
#include "Thermal.hpp"
#include "CeramicDamage.hpp"
#include "SurfaceInformedPolymer.hpp"
#include "StrainHardeningPolymer.hpp"
#include "Graphite.hpp"
#include "Hyperelastic.hpp"
#include "HyperelasticMMS.hpp"
#include "Geomechanics.hpp"
#include "Chiumenti.hpp"
#include "VonMisesJ.hpp"
#include "PerfectlyPlastic.hpp"
#include "ElasticTransverseIsotropic.hpp"

namespace geos
{
namespace constitutive
{

#define REGISTER_THERMAL_SOLID( BASE ) \
  using Thermal##BASE = Thermal< BASE >; \
  using LinearThermal##BASE = Thermal< BASE, LinearThermalLaw >; \
  REGISTER_CATALOG_ENTRY( ConstitutiveBase, Thermal##BASE, string const &, dataRepository::Group * const ) \
  REGISTER_CATALOG_ENTRY( ConstitutiveBase, LinearThermal##BASE, string const &, dataRepository::Group * const )

REGISTER_THERMAL_SOLID( ElasticIsotropic )
REGISTER_THERMAL_SOLID( CeramicDamage )
REGISTER_THERMAL_SOLID( SurfaceInformedPolymer )
REGISTER_THERMAL_SOLID( StrainHardeningPolymer )
REGISTER_THERMAL_SOLID( Graphite )
REGISTER_THERMAL_SOLID( Hyperelastic )
REGISTER_THERMAL_SOLID( HyperelasticMMS )
REGISTER_THERMAL_SOLID( Geomechanics )
REGISTER_THERMAL_SOLID( Chiumenti )
REGISTER_THERMAL_SOLID( VonMisesJ )
REGISTER_THERMAL_SOLID( PerfectlyPlastic )
REGISTER_THERMAL_SOLID( ElasticTransverseIsotropic )

#undef REGISTER_THERMAL_SOLID
} // namespace constitutive
} // namespace geos
