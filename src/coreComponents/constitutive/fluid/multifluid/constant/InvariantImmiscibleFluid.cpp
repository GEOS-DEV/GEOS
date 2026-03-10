#include "InvariantImmiscibleFluid.hpp"
#include "common/format/StringUtilities.hpp"
#include "constitutive/ConstitutiveBase.hpp"

namespace geos
{
namespace constitutive
{

InvariantImmiscibleFluid::InvariantImmiscibleFluid( string const & name, Group * parent )
  : MultiFluidBase( name, parent )
{
  // Override input flags for mandatory options

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "List of fluid components (e.g. CH4, H2O, C5H12)" );

  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "List of fluid phases (e.g. gas, water, oil)" );

  registerWrapper( viewKeyStruct::componentMolarWeightString(), &m_componentMolarWeight )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Component molar weights [kg/mol]" );

  // Densities: constant phase densities
  registerWrapper( "densities", &m_densities )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Constant phase densities" );

  // Viscosities: constant phase viscosities
  registerWrapper( "viscosities", &m_viscosities )
    .setInputFlag( dataRepository::InputFlags::REQUIRED )
    .setDescription( "Constant phase viscosities" );
}

integer InvariantImmiscibleFluid::getWaterPhaseIndex() const
{
  // find index of water phase in user ordering
  for( size_t i=0; i< m_phaseNames.size(); ++i )
  {
    if( stringutilities::toLower( m_phaseNames[i] ) == "water" )
      return static_cast< integer >(i);
  }
  return -1;
}

void InvariantImmiscibleFluid::postInputInitialization()
{
  // check base inputs
  MultiFluidBase::postInputInitialization();

  integer numPhase = numFluidPhases();
  // check densities and viscosities size
  GEOS_THROW_IF_NE_MSG( m_densities.size(), numPhase,
                        GEOS_FMT( "%s: 'Densities' must have %d values", getFullName(), numPhase ), InputError );
  GEOS_THROW_IF_NE_MSG( m_viscosities.size(), numPhase,
                        GEOS_FMT( "%s: 'Viscosities' must have %d values", getFullName(), numPhase ), InputError );

  integer numComponents = numFluidComponents();
  // check tacit assumption of one component per phase
  GEOS_THROW_IF_NE_MSG( numComponents, numPhase,
                        GEOS_FMT( "%d number of components must be equato to %d number of phases", getFullName(), numPhase ), InputError );
}

InvariantImmiscibleFluid::KernelWrapper InvariantImmiscibleFluid::createKernelWrapper() const
{
  return KernelWrapper(
    m_densities.toViewConst(),
    m_viscosities.toViewConst(),
    m_componentMolarWeight.toViewConst(),
    m_useMass,
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseFraction).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseDensity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseMassDensity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseViscosity).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseEnthalpy).toView(),
    const_cast< MultiFluidBase::PhaseProp & >(m_phaseInternalEnergy).toView(),
    const_cast< MultiFluidBase::PhaseComp & >(m_phaseCompFraction).toView(),
    const_cast< MultiFluidBase::FluidProp & >(m_totalDensity).toView()
    );
}

void InvariantImmiscibleFluid::checkTablesParameters( real64 pressure, real64 temperature ) const
{
  // No PVT tables for this model, so nothing to check
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
}

// Register this fluid model
REGISTER_CATALOG_ENTRY( ConstitutiveBase, InvariantImmiscibleFluid, string const &, geos::dataRepository::Group * const )

} // namespace constitutive
} // namespace geos
