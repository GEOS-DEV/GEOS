#include "FrictionDriverRunTest.hpp"
#include "constitutive/contact/CoulombFriction.hpp"
#include "constitutive/contact/FrictionlessContact.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include <type_traits>


namespace geos
{

template
void
FrictionDriver::runTest( constitutive::CoulombFriction &, SolidMechanicsAugmentedLagrangianContact&, const arrayView2d< real64 > & );

template
void
FrictionDriver::runTest( constitutive::FrictionlessContact &, SolidMechanicsAugmentedLagrangianContact&,const arrayView2d< real64 > & );

template
void
FrictionDriver::runTest( constitutive::RateAndStateFriction< std::integral_constant< bool, true > > &, SolidMechanicsAugmentedLagrangianContact&,const arrayView2d< real64 > & );

template
void
FrictionDriver::runTest( constitutive::RateAndStateFriction< std::integral_constant< bool, false > > &, SolidMechanicsAugmentedLagrangianContact&,const arrayView2d< real64 > & );

}
