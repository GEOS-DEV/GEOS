/* SPDX-License-Identifier: LGPL-2.1-only */
#include "SolidMechanicsMPM.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/ThermalLaw.hpp"
#include "physicsSolvers/solidMechanics/kernels/ExplicitMPM.hpp"
#include "common/MpiWrapper.hpp"
#include <vector>

namespace geos
{
namespace
{
using namespace dataRepository;
namespace thermo = constitutive::thermomechanics;

/** Save every particle-indexed constitutive array, including old/new stress,
 * plastic strain, yield strength, and adapter history. clone() is NOT a deep
 * snapshot in this repository. Indexed host serialization makes an independent
 * byte snapshot and restores only this particle, without overwriting neighbors.
 */
class ConstitutiveParticleTrial
{
public:
  ConstitutiveParticleTrial( constitutive::ContinuumBase & model, localIndex const p )
  {
    m_indices.resize( 1 ); m_indices[0] = p;
    model.forWrappers( [&]( WrapperBase & wrapper )
    {
      if( wrapper.numArrayDims() == 0 || wrapper.sizedFromParent() != 1 ) return;
      wrapper.move( hostMemorySpace, true );
      parallelDeviceEvents events;
      buffer_unit_type * cursor = nullptr;
      localIndex const size = wrapper.packByIndex< false >( cursor, m_indices.toViewConst(), false, false, events );
      Entry entry{ &wrapper, std::vector< buffer_unit_type >( size ) };
      cursor = entry.bytes.data();
      wrapper.packByIndex< true >( cursor, m_indices.toViewConst(), false, false, events );
      m_saved.push_back( std::move( entry ) );
    } );
  }
  void restore()
  {
    for( auto & entry : m_saved )
    {
      parallelDeviceEvents events;
      buffer_unit_type const * cursor = entry.bytes.data();
      entry.wrapper->unpackByIndex( cursor, m_indices.toViewConst(), false, false, events );
    }
  }
private:
  struct Entry { WrapperBase * wrapper; std::vector< buffer_unit_type > bytes; };
  array1d< localIndex > m_indices;
  std::vector< Entry > m_saved;
};

thermo::Vector unrotateStress( arraySlice1d< real64 const > const stress,
                          arraySlice2d< real64 const > const rotation )
{
  real64 rt[3][3] = {}, result[6] = {};
  LvArray::tensorOps::transpose< 3, 3 >( rt, rotation );
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( result, rt, stress );
  thermo::Vector v{};
  for( integer i = 0; i < 6; ++i ) v[i] = result[i];
  return v;
}
}

/** Explicit global MPM dynamics, with implicit/local thermoplastic integration.
 * This energy contract is for constant-stiffness, corotational small-elastic-
 * strain solids. A finite-strain hyperelastic, damage, EOS, or variable-modulus
 * model must provide a different stored-energy contract, not inherit this one.
 */
void SolidMechanicsMPM::updateThermomechanicalStress( real64 const dt, ParticleManager & particles )
{
  GEOS_MARK_FUNCTION;
  parallelDeviceSync();
  namespace thermo = constitutive::thermomechanics;
  auto & settings = getGroup< MPMThermalSettings >( "Thermal" );
  GEOS_ERROR_IF( !std::isfinite( dt ) || !( dt > 0.0 ), "Thermomechanical update requires a positive timestep." );
  real64 mechanicalWork = 0.0, plasticHeat = 0.0, reservoirWork = 0.0, maxResidual = 0.0;
  particles.forParticleSubRegions( [&]( ParticleSubRegion & sub )
  {
    for( char const * name : { "particleRank", "particleDeleteFlag", "particleMass", "particleReferenceVolume",
      "particleTemperature", "particleTemperatureRate", "particleHeatCapacity", "particleStress", "particleDeformationGradient",
      "particleRotation", "particleOldRotation", "particleVelocityGradient", "particleInternalEnergy",
      "particleRetainedPlasticEnergy", "particlePlasticWork", "particleThermalIsothermal", "particleThermalPreviousTemperature",
      "particleThermalReservoirEnergy", "particleThermalStrainIncrement", "particleThermalStretch" } )
      sub.getWrapperBase( name ).move( hostMemorySpace, true );
    sub.setActiveParticleIndices();
    if( sub.activeParticleIndices().size() == 0 ) return;
    string const & name = sub.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    auto & model = getConstitutiveModel< constitutive::ContinuumBase >( sub, name );
    auto * material = dynamic_cast< constitutive::ThermalMaterial * >( &model );
    GEOS_ERROR_IF( material == nullptr || !material->supportsThermomechanicalEnergy(),
                   "Material '" << name << "' needs a thermomechanical stored-energy adapter." );
    model.forWrappers( []( WrapperBase & wrapper )
    { if( wrapper.numArrayDims() > 0 ) wrapper.move( hostMemorySpace, true ); } );
    arrayView1d< real64 > const temperature = sub.getParticleTemperature();
    arrayView1d< real64 > const rate = sub.getParticleTemperatureRate();
    arrayView1d< real64 > const capacity = sub.getField< fields::mpm::particleHeatCapacity >();
    arrayView1d< real64 > const energy = sub.getField< fields::mpm::particleInternalEnergy >();
    arrayView1d< real64 > const retained = sub.getField< fields::mpm::particleRetainedPlasticEnergy >();
    arrayView1d< real64 > const workHistory = sub.getField< fields::mpm::particlePlasticWork >();
    arrayView1d< real64 > const previous = sub.getField< fields::mpm::particleThermalPreviousTemperature >();
    arrayView1d< real64 > const reservoir = sub.getField< fields::mpm::particleThermalReservoirEnergy >();
    arrayView1d< real64 > const strainThermal = sub.getField< fields::mpm::particleThermalStrainIncrement >();
    arrayView1d< real64 > const stretch = sub.getField< fields::mpm::particleThermalStretch >();
    arrayView1d< real64 const > const mass = sub.getField< fields::mpm::particleMass >();
    arrayView1d< real64 const > const referenceVolume = sub.getField< fields::mpm::particleReferenceVolume >();
    arrayView1d< integer const > const isothermal = sub.getField< fields::mpm::particleThermalIsothermal >();
    arrayView2d< real64 > const stress = sub.getField< fields::mpm::particleStress >();
    arrayView3d< real64 const > const f = sub.getField< fields::mpm::particleDeformationGradient >();
    arrayView3d< real64 const > const r = sub.getField< fields::mpm::particleRotation >();
    arrayView3d< real64 const > const oldR = sub.getField< fields::mpm::particleOldRotation >();
    arrayView3d< real64 const > const l = sub.getField< fields::mpm::particleVelocityGradient >();
    arrayView1d< real64 const > alpha;
    if( material->solverThermalExpansion() )
    {
      model.getWrapperBase( "thermalExpansionCoefficient" ).move( hostMemorySpace, true );
      alpha = model.getReference< array1d< real64 > >( "thermalExpansionCoefficient" ).toViewConst();
    }
    auto const active = sub.activeParticleIndices();
    for( localIndex pp = 0; pp < active.size(); ++pp )
    {
      localIndex const p = active[pp];
      real64 const predictor = temperature[p], initialTemperature = previous[p];
      real64 const initialEnergy = energy[p], retainedOld = retained[p];
      real64 const expansion = alpha.size() > 0 ? alpha[p] : 0.0;
      real64 const rho0 = mass[p]/referenceVolume[p];
      auto const compliance = material->thermomechanicalCompliance( p );
      auto const oldStress = unrotateStress( stress[p], oldR[p] );
      real64 d[6] = { dt*l[p][0][0], dt*l[p][1][1], dt*l[p][2][2],
        0.5*dt*( l[p][1][2]+l[p][2][1] ), 0.5*dt*( l[p][0][2]+l[p][2][0] ),
        0.5*dt*( l[p][0][1]+l[p][1][0] ) };
      real64 rt[3][3] = {}, rotated[6] = {};
      LvArray::tensorOps::transpose< 3, 3 >( rt, oldR[p] );
      LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( rotated, rt, d );
      thermo::Vector strain{};
      for( integer i = 0; i < 6; ++i ) strain[i] = rotated[i]*( i < 3 ? 1.0 : 2.0 );
      SortedArray< localIndex > single;
      single.insert( p );
      ConstitutiveParticleTrial snapshot( model, p );
      thermo::EnergyState finalState;
      thermo::WorkIncrement work{};
      real64 residual = 0.0;
      auto trial = [&]( real64 const t )
      {
        snapshot.restore();
        temperature[p] = t;
        strainThermal[p] = expansion*( t-initialTemperature );
        material->prepareThermomechanicalTrial( p, t );
        // Other material adapters may consume temperature as a constitutive
        // field. Do not run the global dependency updater on every local trial.
        if( model.hasWrapper( "temperature" ) )
          model.getReference< array1d< real64 > >( "temperature" )[p] = t;
        constitutive::ConstitutivePassThruMPM< constitutive::ContinuumBase >::execute( model, [&]( auto & solid )
        {
          auto wrapper = solid.createKernelUpdates();
          solidMechanicsMPMKernels::ParticleStateUpdateKernel::launch< serialPolicy >(
            single.toViewConst(), 1, wrapper, dt, f, r, oldR, l, stress, strainThermal );
        } );
        auto const newStress = unrotateStress( stress[p], r[p] );
        work = thermo::work( compliance, oldStress, newStress, strain, rho0, strainThermal[p] );
        real64 const retainedNew = retainedOld+( 1.0-material->plasticHeatFraction() )*work.plastic;
        finalState = thermo::state( compliance, newStress, rho0, expansion, retainedNew );
        residual = material->internalEnergy( t, finalState )-initialEnergy-work.total;
        if( !std::isfinite( t ) || t < 0.0 || !( material->thermalSpecificHeat( t ) > 0.0 ) ||
            !( material->stateHeatCapacity( t, finalState ) > 0.0 ) )
          throw std::domain_error( "Inadmissible thermomechanical temperature/heat capacity" );
        return residual;
      };
      real64 finalTemperature = predictor;
      real64 bath = 0.0;
      try
      {
        if( isothermal[p] != 0 )
        {
          trial( predictor );
          bath = residual;
        }
        else if( settings.mechanicalCoupling == "Lagged" )
        {
          trial( predictor );
          finalTemperature = material->temperatureFromInternalEnergy( initialEnergy+work.total, predictor, finalState );
        }
        else
        {
          real64 const tolerance = settings.couplingTolerance*std::max( 1.0,
            std::abs( material->thermalSpecificHeat( predictor )*std::max( 1.0, predictor ) ) );
          finalTemperature = thermo::solveTemperature( predictor, tolerance, settings.maxCouplingIterations, trial );
          // The last bracket evaluation is not necessarily the root. Replay
          // exactly that root; all history writes from previous trials vanish.
          trial( finalTemperature );
          maxResidual = std::max( maxResidual, std::abs( residual ) );
        }
      }
      catch( std::exception const & error )
      {
        snapshot.restore();
        temperature[p] = predictor;
        GEOS_ERROR( "Thermomechanical update failed for material '" << name << "', particle local index " << p << ": " << error.what() );
      }
      real64 const conductivity = material->thermalConductivity( finalTemperature );
      GEOS_ERROR_IF( !std::isfinite( conductivity ) || conductivity < 0.0 ||
                     !std::isfinite( finalTemperature ) || finalTemperature < 0.0,
                     "Accepted thermomechanical temperature is outside the thermal material range." );
      // The accepted trial has already saved the constitutive state once. Energy
      // includes stress work exactly once; plastic heat only partitions that work.
      temperature[p] = finalTemperature;
      retained[p] = finalState.retained;
      workHistory[p] += work.plastic;
      energy[p] = material->internalEnergy( finalTemperature, finalState );
      reservoir[p] += bath;
      rate[p] = ( finalTemperature-initialTemperature )/dt;
      capacity[p] = material->stateHeatCapacity( finalTemperature, finalState );
      real64 const consumedTemperature = settings.mechanicalCoupling == "Lagged" ? predictor : finalTemperature;
      strainThermal[p] = expansion*( consumedTemperature-initialTemperature );
      stretch[p] *= std::exp( strainThermal[p] );
      previous[p] = consumedTemperature;
      mechanicalWork += mass[p]*work.total;
      plasticHeat += mass[p]*material->plasticHeatFraction()*work.plastic;
      reservoirWork += mass[p]*bath;
      if( model.hasWrapper( "temperature" ) ) model.getReference< array1d< real64 > >( "temperature" )[p] = finalTemperature;
      if( model.hasWrapper( "temperatureRate" ) ) model.getReference< array1d< real64 > >( "temperatureRate" )[p] = rate[p];
    }
    // Mark touched host arrays before downstream device mechanics consume them.
    for( char const * key : { "particleStress", "particleTemperature", "particleTemperatureRate", "particleHeatCapacity",
                             "particleInternalEnergy", "particleThermalStrainIncrement", "particleThermalStretch" } )
      sub.getWrapperBase( key ).move( parallelDeviceMemorySpace, true );
    model.forWrappers( []( WrapperBase & wrapper )
    { if( wrapper.numArrayDims() > 0 ) wrapper.move( parallelDeviceMemorySpace, true ); } );
  } );
  real64 const globalWork = MpiWrapper::sum( mechanicalWork );
  real64 const globalHeat = MpiWrapper::sum( plasticHeat );
  real64 const globalBath = MpiWrapper::sum( reservoirWork );
  settings.cumulativeMechanicalWork += globalWork;
  settings.cumulativePlasticHeat += globalHeat;
  settings.cumulativeReservoirEnergy += globalBath;
  real64 const globalResidual = MpiWrapper::max( maxResidual );
  GEOS_LOG_RANK_0_IF( settings.diagnostics == 1,
    "MPMThermomechanics mechanicalWork=" << globalWork << " plasticHeat=" << globalHeat
    << " reservoirEnergy=" << globalBath << " maxSpecificEnergyResidual=" << globalResidual );
}
} // namespace geos
