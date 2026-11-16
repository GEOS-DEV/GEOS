/* SPDX-License-Identifier: LGPL-2.1-only */
#include "SolidMechanicsMPM.hpp"
#include "MPMThermalNumerics.hpp"
#include "constitutive/solid/ThermalLaw.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/ParticleRegion.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include <map>
#include <set>
#include <sstream>

namespace geos
{
using namespace dataRepository;
namespace
{
char const * const capacityKey = "gridThermalCapacity";
char const * const momentKey = "gridThermalTemperatureMoment";
char const * const temperatureKey = "gridThermalTemperature";
char const * const fixedCapacityKey = "gridThermalFixedCapacity";
char const * const fixedMomentKey = "gridThermalFixedMoment";
char const * const fixedSquareKey = "gridThermalFixedSquare";
char const * const rateKey = "gridThermalPower";
char const * const boundKey = "gridThermalConductanceBound";
char const * const incrementKey = "gridThermalEnergyIncrement";
char const * const geometryKey = "gridThermalGeometry";
constexpr integer geometryComponents = 13;

struct ThermalContribution
{
  localIndex node;
  localIndex slot;
  real64 weight;
  mpmThermal::Vector gradient;
};

struct ThermalParticle
{
  constitutive::ThermalMaterial const * material;
  real64 * temperature;
  real64 * energy;
  real64 * heatCapacity;
  real64 * temperatureRate;
  real64 * reservoir;
  real64 * strain;
  real64 * stretch;
  real64 * previousTemperature;
  constitutive::thermomechanics::EnergyState mechanicalState;
  integer * isothermal;
  real64 initialTemperature;
  real64 alpha;
  real64 mass;
  real64 volume;
  real64 specificHeat = 0.0;
  real64 conductivity = 0.0;
  bool bspline = false;
  mpmThermal::Vector center;
  mpmThermal::Vector normal;
  mpmThermal::Vector surface;
  std::vector< ThermalContribution > support;
  std::vector< integer > membership;
  std::vector< real64 > referenceMass;
};

struct ThermalInterface
{
  localIndex a;
  localIndex b;
  real64 conductance;
  bool perfect;
  bool bonded;
};

void moveParticleThermalData( ParticleSubRegion & sub )
{
  stdVector< string > const names = {
    "particleRank", "particleDeleteFlag", "particleID", "particleGroup", "particleColor",
    "particleCenter", "particleVolume", "particleTemperature", "particleTemperatureRate",
    "particleMass", "particleHeatCapacity", "particleSurfaceNormal", "particleSurfacePosition",
    "particleInternalEnergy", "particleRetainedPlasticEnergy", "particlePlasticWork", "particleThermalIsothermal",
    "particleReferenceVolume", "particleStress", "particleRotation", "particleThermalMaterialIndex", "particleThermalInitialized",
    "particleThermalPreviousTemperature", "particleThermalReservoirEnergy",
    "particleThermalStrainIncrement", "particleThermalStretch",
    "particleThermalConditionMembership", "particleThermalConditionReferenceMass" };
  for( auto const & name : names ) sub.getWrapperBase( name ).move( hostMemorySpace, true );
  // Rebuild the same sorted owned-particle list on host, preserving active ordinals.
  sub.setActiveParticleIndices();
}

bool finiteTemperature( real64 const t ) { return std::isfinite( t ) && t >= 0.0; }
}

bool SolidMechanicsMPM::thermalEnabled() const
{ return getGroup< MPMThermalSettings >( "Thermal" ).enabled == 1; }

void SolidMechanicsMPM::registerThermalParticleFields( ParticleSubRegion & sub )
{
  if( !thermalEnabled() ) return;
  using namespace fields::mpm;
  localIndex count = 0;
  getGroup< MPMThermalSettings >( "Thermal" ).forSubGroups< MPMThermalParticleCondition >(
    [&]( MPMThermalParticleCondition const & ) { ++count; } );
  sub.registerField< particleRetainedPlasticEnergy >( getName() );
  sub.registerField< particlePlasticWork >( getName() );
  sub.registerField< particleThermalIsothermal >( getName() );
  sub.registerField< particleThermalMaterialIndex >( getName() );
  sub.registerField< particleThermalInitialized >( getName() );
  sub.registerField< particleThermalPreviousTemperature >( getName() );
  sub.registerField< particleThermalReservoirEnergy >( getName() );
  sub.registerField< particleThermalStrainIncrement >( getName() );
  sub.registerField< particleThermalStretch >( getName() );
  sub.registerField< particleThermalConditionMembership >( getName() ).reference().resizeDimension< 1 >( count );
  sub.registerField< particleThermalConditionReferenceMass >( getName() ).reference().resizeDimension< 1 >( count );
}

void SolidMechanicsMPM::registerThermalGridFields( NodeManager & nodes )
{
  if( !thermalEnabled() ) return;
  for( char const * name : { capacityKey, momentKey, temperatureKey, fixedCapacityKey,
                            fixedMomentKey, fixedSquareKey, rateKey, boundKey, incrementKey } )
    nodes.registerWrapper< array2d< real64 > >( name ).setInputFlag( InputFlags::FALSE ).
      setRestartFlags( RestartFlags::NO_WRITE ).setPlotLevel( PlotLevel::NOPLOT );
  nodes.registerWrapper< array3d< real64 > >( geometryKey ).setInputFlag( InputFlags::FALSE ).
    setRestartFlags( RestartFlags::NO_WRITE ).setPlotLevel( PlotLevel::NOPLOT );
}

/** Host reference implementation. Explicit host/device transitions keep this
 * path usable with HIP/CUDA mechanics, without capturing polymorphic materials
 * or large contact graphs in a mechanical device kernel. Only owned particles
 * scatter. Additive data reduce ghost->owner; solved increments replace ghosts.
 */
void SolidMechanicsMPM::advanceThermalState( real64 const time, real64 const dt,
                                           DomainPartition & domain, ParticleManager & particles,
                                           NodeManager & nodes, MeshLevel & mesh )
{
  if( !thermalEnabled() ) return;
  GEOS_MARK_FUNCTION;
  GEOS_ERROR_IF( !( dt > 0.0 ) || !std::isfinite( dt ), "Thermal transport requires a finite positive time step." );
  parallelDeviceSync();
  MPMThermalSettings & settings = getGroup< MPMThermalSettings >( "Thermal" );
  // Old checkpoints carry a caloric energy plus an unrelated work ledger. Those
  // values cannot be reinterpreted as the new total-energy state without a
  // material-specific conversion. Fail explicitly instead of corrupting energy.
  if( settings.energyFormulationVersion == 0 )
  {
    integer previouslyInitialized = 0;
    particles.forParticleSubRegions( [&]( ParticleSubRegion & sub )
    {
      sub.getWrapperBase( fields::mpm::particleThermalInitialized::key() ).move( hostMemorySpace, true );
      arrayView1d< integer const > const flags = sub.getField< fields::mpm::particleThermalInitialized >();
      for( localIndex p = 0; p < flags.size(); ++p ) previouslyInitialized = std::max( previouslyInitialized, flags[p] );
    } );
    GEOS_ERROR_IF( MpiWrapper::max( previouslyInitialized ) != 0,
                   "This checkpoint predates unified thermomechanical energy. Start a fresh run; automatic energy conversion is unavailable." );
    settings.energyFormulationVersion = 2;
  }
  GEOS_ERROR_IF( settings.energyFormulationVersion != 2, "Unknown thermomechanical energy checkpoint version." );
  std::vector< MPMThermalParticleCondition const * > conditions;
  settings.forSubGroups< MPMThermalParticleCondition >(
    [&]( MPMThermalParticleCondition const & c ) { conditions.push_back( &c ); } );
  std::vector< string > materials;
  particles.forParticleSubRegions( [&]( ParticleSubRegion & sub )
  { materials.push_back( sub.getReference< string >( viewKeyStruct::solidMaterialNamesString() ) ); } );
  std::sort( materials.begin(), materials.end() );
  materials.erase( std::unique( materials.begin(), materials.end() ), materials.end() );
  string_array conditionNames;
  for( auto const * c : conditions ) conditionNames.push_back( c->getName() );
  auto checkNames = []( string_array & saved, auto const & current )
  {
    if( saved.empty() ) saved.assign( current.begin(), current.end() );
    GEOS_ERROR_IF( saved.size() != current.size() || !std::equal( saved.begin(), saved.end(), current.begin() ),
                   "Thermal condition/material names changed across restart; particle state indices would be ambiguous." );
  };
  checkNames( settings.conditionNames, conditionNames ); checkNames( settings.materialNames, materials );
  integer const materialCount = static_cast< integer >( materials.size() );
  GEOS_ERROR_IF( materialCount == 0, "Thermal MPM requires particle materials." );
  settings.forSubGroups< MPMThermalContactRule >( [&]( MPMThermalContactRule const & rule )
  {
    for( auto const & name : rule.materials )
      GEOS_ERROR_IF( std::find( materials.begin(), materials.end(), name ) == materials.end(),
                     "Thermal contact rule names unknown particle material '" << name << "'." );
    for( auto const group : rule.contactGroups )
      GEOS_ERROR_IF( group >= m_numContactGroups, "Thermal contact rule names an unknown contact group." );
  } );
  // Regions/material names are global input metadata. Check their ordering across
  // ranks before allocating field slots used by ordinary node communication.
  integer signature = 0;
  for( auto const & name : materials ) for( unsigned char const c : name ) signature = ( signature*31LL+c )%1000000007;
  GEOS_ERROR_IF( MpiWrapper::min( materialCount ) != MpiWrapper::max( materialCount ) ||
                 MpiWrapper::min( signature ) != MpiWrapper::max( signature ),
                 "Particle material metadata must have identical names on every rank." );
  integer fieldCount = m_numContactGroups * ( m_damageFieldPartitioning == 1 ? 2 : 1 );
  if( m_rigidBodyMode == 1 )
  {
    m_rigidBodyColors.move( hostMemorySpace, true );
    fieldCount = static_cast< integer >( m_rigidBodyColors.size() );
    GEOS_ERROR_IF( fieldCount == 0, "Rigid thermal transport requires initialized body colors." );
  }
  GEOS_ERROR_IF( fieldCount > std::numeric_limits< integer >::max()/materialCount,
                 "Too many thermal field/material slots." );
  localIndex const slots = fieldCount*materialCount;
  localIndex const nodeCount = nodes.size();
  for( char const * name : { capacityKey, momentKey, temperatureKey, fixedCapacityKey,
                            fixedMomentKey, fixedSquareKey, rateKey, boundKey, incrementKey } )
  {
    auto & storage = nodes.getReference< array2d< real64 > >( name );
    storage.move( hostMemorySpace, true ); storage.resize( nodeCount, slots ); storage.zero();
  }
  auto & geometryStorage = nodes.getReference< array3d< real64 > >( geometryKey );
  geometryStorage.move( hostMemorySpace, true ); geometryStorage.resize( nodeCount, slots, geometryComponents ); geometryStorage.zero();
  auto const capacity = nodes.getReference< array2d< real64 > >( capacityKey ).toView();
  auto const moment = nodes.getReference< array2d< real64 > >( momentKey ).toView();
  auto const temperature = nodes.getReference< array2d< real64 > >( temperatureKey ).toView();
  auto const fixedCapacity = nodes.getReference< array2d< real64 > >( fixedCapacityKey ).toView();
  auto const fixedMoment = nodes.getReference< array2d< real64 > >( fixedMomentKey ).toView();
  auto const fixedSquare = nodes.getReference< array2d< real64 > >( fixedSquareKey ).toView();
  auto const rate = nodes.getReference< array2d< real64 > >( rateKey ).toView();
  auto const bound = nodes.getReference< array2d< real64 > >( boundKey ).toView();
  auto const increment = nodes.getReference< array2d< real64 > >( incrementKey ).toView();
  auto const geometry = geometryStorage.toView();
  nodes.getWrapperBase( NodeManager::viewKeyStruct::referencePositionString() ).move( hostMemorySpace, true );
  auto const position = nodes.referencePosition();
  nodes.getWrapperBase( "ghostRank" ).move( hostMemorySpace, true );
  auto const ghostRank = static_cast< NodeManager const & >( nodes ).ghostRank();
  std::vector< ThermalParticle > point;
  real64 sourceEnergy = 0.0, reservoirEnergy = 0.0, initialEnergy = 0.0;
  localIndex subIndex = 0;
  particles.forParticleSubRegions( [&]( ParticleSubRegion & sub )
  {
    moveParticleThermalData( sub );
    string const materialName = sub.getReference< string >( viewKeyStruct::solidMaterialNamesString() );
    integer const material = static_cast< integer >( std::lower_bound( materials.begin(), materials.end(), materialName )-materials.begin() );
    auto & model = getConstitutiveModel< constitutive::ContinuumBase >( sub, materialName );
    auto const * thermal = dynamic_cast< constitutive::ThermalMaterial const * >( &model );
    auto const active = sub.activeParticleIndices();
    GEOS_ERROR_IF( thermal == nullptr && active.size() > 0,
                   "Thermal transport requires a Thermal/LinearThermal solid for particle material '" << materialName << "'." );
    GEOS_ERROR_IF( model.hasWrapper( "specificInternalEnergy" ) && active.size() > 0,
                   "An EOS owning specificInternalEnergy needs a caloric adapter before use with this thermal solver." );
    GEOS_ERROR_IF( thermal != nullptr && !thermal->supportsThermomechanicalEnergy() && active.size() > 0,
                   "Material '" << materialName << "' has no thermomechanical stored-energy contract. "
                   "Use Thermal/LinearThermal ElasticIsotropic or PerfectlyPlastic. "
                   "Anisotropic, damage, polymer, VonMisesJ, and finite-strain models require their own stored-energy adapters." );
    // Constitutive arrays are read by the host caloric closure. Move them before
    // obtaining any views; retain views, never temporary owning-array copies.
    model.forWrappers( []( WrapperBase & wrapper )
    { if( wrapper.numArrayDims() > 0 ) wrapper.move( hostMemorySpace, true ); } );
    auto const t = sub.getParticleTemperature();
    auto const tRate = sub.getParticleTemperatureRate();
    // getField() returns an owning array by reference. Keep views into the
    // persistent fields: ThermalParticle retains their element pointers after
    // this subregion callback returns, throughout the thermal substeps.
    arrayView1d< real64 const > const mass = sub.getField< fields::mpm::particleMass >();
    auto const volume = sub.getParticleVolume();
    arrayView1d< real64 > const heatCapacity = sub.getField< fields::mpm::particleHeatCapacity >();
    arrayView1d< real64 > const energy = sub.getField< fields::mpm::particleInternalEnergy >();
    arrayView1d< real64 const > const retained = sub.getField< fields::mpm::particleRetainedPlasticEnergy >();
    arrayView1d< real64 const > const referenceVolume = sub.getField< fields::mpm::particleReferenceVolume >();
    arrayView2d< real64 const > const stress = sub.getField< fields::mpm::particleStress >();
    arrayView3d< real64 const > const rotation = sub.getField< fields::mpm::particleRotation >();
    arrayView1d< integer > const isothermal = sub.getField< fields::mpm::particleThermalIsothermal >();
    arrayView1d< integer > const oldMaterial = sub.getField< fields::mpm::particleThermalMaterialIndex >();
    arrayView1d< integer > const initializedFlag = sub.getField< fields::mpm::particleThermalInitialized >();
    arrayView1d< real64 > const previousTemperature = sub.getField< fields::mpm::particleThermalPreviousTemperature >();
    arrayView1d< real64 > const reservoir = sub.getField< fields::mpm::particleThermalReservoirEnergy >();
    arrayView1d< real64 > const strain = sub.getField< fields::mpm::particleThermalStrainIncrement >();
    arrayView1d< real64 > const stretch = sub.getField< fields::mpm::particleThermalStretch >();
    arrayView2d< integer > const membership = sub.getField< fields::mpm::particleThermalConditionMembership >();
    arrayView2d< real64 > const referenceMass = sub.getField< fields::mpm::particleThermalConditionReferenceMass >();
    GEOS_ERROR_IF( membership.size( 1 ) != static_cast< localIndex >( conditions.size() ),
                   "Thermal condition layout changed across restart; retain the original condition names/order." );
    auto const ids = sub.getParticleID();
    auto const groups = sub.getParticleGroup();
    auto const colors = sub.getParticleColor();
    auto const center = sub.getParticleCenter();
    auto const normal = sub.getParticleSurfaceNormal();
    auto const surface = sub.getParticleSurfacePosition();
    auto const & region = dynamicCast< ParticleRegion const & >( sub.getParent().getParent() );
    auto & countsStorage = m_numEffectiveMappedNodes[subIndex]; countsStorage.move( hostMemorySpace, true );
    auto & mapStorage = m_effectiveMappedNodes[subIndex]; mapStorage.move( hostMemorySpace, true );
    auto & fieldStorage = m_effectiveMappedFields[subIndex]; fieldStorage.move( hostMemorySpace, true );
    auto & weightStorage = m_effectiveShapeFunctionValues[subIndex]; weightStorage.move( hostMemorySpace, true );
    auto & gradientStorage = m_effectiveShapeFunctionGradientValues[subIndex]; gradientStorage.move( hostMemorySpace, true );
    arrayView1d< real64 const > alpha;
    if( thermal != nullptr && thermal->solverThermalExpansion() )
    {
      model.getWrapperBase( "thermalExpansionCoefficient" ).move( hostMemorySpace, true );
      alpha = model.getReference< array1d< real64 > >( "thermalExpansionCoefficient" ).toViewConst();
    }
    for( localIndex pp = 0; pp < active.size(); ++pp )
    {
      localIndex const p = active[pp];
      GEOS_ERROR_IF( !( mass[p] > 0.0 ) || !( volume[p] > 0.0 ) || !finiteTemperature( t[p] ),
                     "Invalid particle mass, volume, or absolute temperature in thermal MPM." );
      constitutive::thermomechanics::Vector referenceStress{};
      real64 rt[3][3] = {}, unrotated[6] = {};
      LvArray::tensorOps::transpose< 3, 3 >( rt, rotation[p] );
      LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( unrotated, rt, stress[p] );
      for( integer i = 0; i < 6; ++i ) referenceStress[i] = unrotated[i];
      real64 const expansion = alpha.size() > 0 ? alpha[p] : 0.0;
      auto const storedState = constitutive::thermomechanics::state( thermal->thermomechanicalCompliance( p ),
        referenceStress, mass[p]/referenceVolume[p], expansion, retained[p] );
      bool const initialized = initializedFlag[p] != 0;
      if( !initialized )
      {
        energy[p] = thermal->internalEnergy( t[p], storedState ); previousTemperature[p] = t[p]; stretch[p] = 1.0;
        initializedFlag[p] = 1;
      }
      initialEnergy += mass[p]*energy[p];
      if( initialized && oldMaterial[p] != material )
      {
        // Specific thermal energy follows material-following particle fields.
        t[p] = thermal->temperatureFromInternalEnergy( energy[p], t[p], storedState );
      }
      else if( initialized )
      {
        // Account explicitly for temperature prescribed by an existing MPM event.
        real64 const imposedEnergy = thermal->internalEnergy( t[p], storedState );
        reservoir[p] += imposedEnergy-energy[p];
        reservoirEnergy += mass[p]*( imposedEnergy-energy[p] );
        energy[p] = imposedEnergy;
      }
      oldMaterial[p] = material; strain[p] = 0.0;
      ThermalParticle q;
      q.mechanicalState = storedState; q.isothermal = &isothermal[p]; isothermal[p] = 0;
      q.material = thermal; q.temperature = &t[p]; q.energy = &energy[p]; q.heatCapacity = &heatCapacity[p];
      q.temperatureRate = &tRate[p]; q.reservoir = &reservoir[p]; q.strain = &strain[p]; q.stretch = &stretch[p];
      q.previousTemperature = &previousTemperature[p];
      q.initialTemperature = previousTemperature[p]; q.alpha = alpha.size() > 0 ? alpha[p] : 0.0;
      GEOS_ERROR_IF( !std::isfinite( q.alpha ), "Thermal expansion coefficient must be finite." );
      q.mass = mass[p]; q.volume = volume[p];
      q.bspline = sub.getParticleType() == ParticleType::SinglePointBSpline;
      real64 const initialCapacity = thermal->stateHeatCapacity( t[p], storedState );
      GEOS_ERROR_IF( !finiteTemperature( t[p] ) || !std::isfinite( initialCapacity ) || !( initialCapacity > 0.0 ),
                     "Initial thermal state is outside the caloric material range." );
      for( integer k = 0; k < 3; ++k ) { q.center[k] = center[p][k]; q.normal[k] = normal[p][k]; q.surface[k] = surface[p][k]; }
      integer body = 0;
      if( m_rigidBodyMode == 1 )
      {
        auto it = std::find( m_rigidBodyColors.begin(), m_rigidBodyColors.end(), colors[p] );
        GEOS_ERROR_IF( it == m_rigidBodyColors.end(), "Unknown rigid-body color in thermal mapping." );
        body = static_cast< integer >( it-m_rigidBodyColors.begin() );
      }
      for( localIndex j = 0; j < countsStorage[pp]; ++j )
      {
        real64 const w = weightStorage[pp][j];
        GEOS_ERROR_IF( !std::isfinite( w ) || w < 0.0, "Thermal mapping requires nonnegative finite shape values." );
        integer const f = m_rigidBodyMode == 1 ? body : fieldStorage[pp][j];
        localIndex const node = mapStorage[pp][j];
        GEOS_ERROR_IF( node < 0 || node >= nodeCount || f < 0 || f >= fieldCount,
                       "Invalid active thermal particle-grid map." );
        ThermalContribution a{ node, f*materialCount+material, w, {{ 0.0, 0.0, 0.0 }} };
        for( integer k = 0; k < m_numDims; ++k ) a.gradient[k] = gradientStorage[pp][j][k];
        // Zero-value shape functions can still have a nonzero gradient. Keep
        // them when neighboring particles give the corresponding node capacity.
        if( !( w > 0.0 ) && !( mpmThermal::dot( a.gradient, a.gradient ) > 0.0 ) ) continue;
        q.support.push_back( a );
        capacity[node][a.slot] += w*q.mass*initialCapacity;
      }
      for( std::size_t c = 0; c < conditions.size(); ++c )
      {
        if( !( referenceMass[p][c] > 0.0 ) )
        {
          membership[p][c] = conditions[c]->matches( region.getName(), ids[p], groups[p], colors[p] ) ? 1 : 0;
          referenceMass[p][c] = q.mass;
        }
        q.membership.push_back( membership[p][c] ); q.referenceMass.push_back( referenceMass[p][c] );
      }
      point.push_back( std::move( q ) );
    }
    ++subIndex;
  } );
  syncGridFields( { capacityKey }, domain, nodes, mesh, MPI_SUM );
  // Drop negligible-capacity supports and differentiate the normalized retained
  // weights. This preserves sum(N)=1 and sum(grad N)=0, including particle/node
  // coincidences and mass-cutoff nodes, rather than silently discarding energy.
  for( auto & q : point )
  {
    auto const originalSupport = q.support;
    q.support.erase( std::remove_if( q.support.begin(), q.support.end(), [&]( ThermalContribution const & a )
    {
      real64 largest = 0.0;
      for( localIndex s = 0; s < slots; ++s ) largest = std::max( largest, capacity[a.node][s] );
      return !( capacity[a.node][a.slot] > 0.0 ) || capacity[a.node][a.slot] < settings.minimumCapacityFraction*largest;
    } ), q.support.end() );
    if( std::none_of( q.support.begin(), q.support.end(),
                     []( ThermalContribution const & a ) { return a.weight > 0.0; } ) && !originalSupport.empty() )
    {
      // A very light material must retain its energy and particle sources.
      // Collapse its unresolved support instead of deleting its thermal state.
      q.support.clear();
      q.support.push_back( *std::max_element( originalSupport.begin(), originalSupport.end(),
        []( ThermalContribution const & a, ThermalContribution const & b ) { return a.weight < b.weight; } ) );
    }
    real64 sum = 0.0; mpmThermal::Vector gradSum = {{ 0.0, 0.0, 0.0 }};
    for( auto const & a : q.support ) { sum += a.weight; for( integer k = 0; k < 3; ++k ) gradSum[k] += a.gradient[k]; }
    GEOS_ERROR_IF( !( sum > 0.0 ), "A thermal particle has no positive-capacity grid support." );
    for( auto & a : q.support )
    {
      for( integer k = 0; k < 3; ++k ) a.gradient[k] = ( a.gradient[k]-a.weight*gradSum[k]/sum )/sum;
      a.weight /= sum;
      auto const node = a.node, s = a.slot;
      geometry[node][s][0] += q.mass*a.weight;
      geometry[node][s][1] += q.volume*a.weight;
      if( q.bspline ) geometry[node][s][12] += q.mass*a.weight;
      bool const explicitSurface = mpmThermal::dot( q.normal, q.normal ) > 1.0e-24;
      for( integer k = 0; k < 3; ++k )
      {
        geometry[node][s][2+k] += q.volume*a.gradient[k] +
          m_explicitSurfaceNormalInfluence*q.volume*a.weight*q.normal[k]/m_hEl[k];
        geometry[node][s][5+k] += q.mass*a.weight*( q.center[k]-position[node][k] );
        if( explicitSurface ) geometry[node][s][8+k] += q.mass*a.weight*( q.center[k]+q.surface[k]-position[node][k] );
      }
      if( explicitSurface ) geometry[node][s][11] += q.mass*a.weight;
    }
  }
  syncGridFields( { geometryKey }, domain, nodes, mesh, MPI_SUM );
  // Mechanical fields have already been synchronized and normalized by gridTrialUpdate.
  for( char const * key : { viewKeyStruct::gridDamageString(), viewKeyStruct::gridMaxDamageString(),
                            viewKeyStruct::gridMaterialVolumeString(), viewKeyStruct::gridMappingNormalTensorString(),
                            viewKeyStruct::gridCenterOfMassString(), viewKeyStruct::gridSingleFieldStateFractionString() } )
    nodes.getWrapperBase( key ).move( hostMemorySpace, true );
  auto const damage = nodes.getReference< array2d< real64 > >( viewKeyStruct::gridDamageString() ).toViewConst();
  auto const maxDamage = nodes.getReference< array2d< real64 > >( viewKeyStruct::gridMaxDamageString() ).toViewConst();
  auto const materialVolume = nodes.getReference< array2d< real64 > >( viewKeyStruct::gridMaterialVolumeString() ).toViewConst();
  auto const normalTensor = nodes.getReference< array3d< real64 > >( viewKeyStruct::gridMappingNormalTensorString() ).toViewConst();
  auto const centerOfMass = nodes.getReference< array3d< real64 > >( viewKeyStruct::gridCenterOfMassString() ).toViewConst();
  auto const singleFraction = nodes.getReference< array2d< real64 > >( viewKeyStruct::gridSingleFieldStateFractionString() ).toViewConst();
  if( m_rigidBodyMode == 1 ) m_rigidBodyContactGroups.move( hostMemorySpace, true );
  std::vector< std::vector< ThermalInterface > > interfaces( nodeCount );
  mpmThermal::Vector h = {{ m_hEl[0], m_hEl[1], m_hEl[2] }};
  real64 const hMin = m_planeStrain == 1 ? std::min( h[0], h[1] ) : std::min( h[0], std::min( h[1], h[2] ) );
  auto & partition = dynamicCast< SpatialPartition & >( domain.getGroup( domain.groupKeys.partitionManager ) );
  auto const periodic = partition.getPeriodic();
  for( localIndex node = 0; node < nodeCount; ++node )
  {
    for( localIndex a = 0; a < slots; ++a ) if( geometry[node][a][0] > 0.0 )
    for( localIndex b = a+1; b < slots; ++b ) if( geometry[node][b][0] > 0.0 )
    {
      integer const fa = a/materialCount, fb = b/materialCount;
      integer const ga = m_rigidBodyMode == 1 ? m_rigidBodyContactGroups[fa] : fa%m_numContactGroups;
      integer const gb = m_rigidBodyMode == 1 ? m_rigidBodyContactGroups[fb] : fb%m_numContactGroups;
      auto const * rule = settings.contactRule( materials[a%materialCount], materials[b%materialCount], ga, gb );
      bool bonded = fa == fb && rule == nullptr;
      if( m_rigidBodyMode == 0 && fa != fb && ga == gb )
      {
        real64 const va = materialVolume[node][fa], vb = materialVolume[node][fb];
        real64 const fraction = va+vb > 0.0 ? ( va*singleFraction[node][fa]+vb*singleFraction[node][fb] )/( va+vb ) : 0.0;
        real64 const quality = computeSurfaceQualityFromMappingNormalTensor( m_planeStrain,
          normalTensor[node][fa], normalTensor[node][fb], va, vb );
        bonded = !evaluateSeparabilityCriterion( m_numContactGroups, m_maxSingleFieldStateFractionForSeparability,
          m_separabilityMinDamage, m_thinFeatureDFGThreshold, m_neighborRadius, m_surfaceQualityThreshold,
          fa, fb, damage[node][fa], damage[node][fb], maxDamage[node][fa], maxDamage[node][fb],
          centerOfMass[node][fa], centerOfMass[node][fb], quality, fraction );
      }
      if( bonded ) { interfaces[node].push_back( { a, b, 0.0, true, true } ); continue; }
      string const & mode = rule == nullptr ? settings.defaultContactMode : rule->mode;
      real64 const resistance = rule == nullptr ? settings.defaultContactResistance : rule->resistance;
      if( mode == "Insulated" ) continue;
      mpmThermal::Vector na, nb, normal, pa, pb, separation, pointOnPlane;
      for( integer k = 0; k < 3; ++k ) { na[k] = geometry[node][a][2+k]; nb[k] = geometry[node][b][2+k]; }
      if( m_planeStrain == 1 ) { na[2] = 0.0; nb[2] = 0.0; }
      if( !mpmThermal::normalize( na ) || !mpmThermal::normalize( nb ) || mpmThermal::dot( na, nb ) > -0.05 ) continue;
      for( integer k = 0; k < 3; ++k ) normal[k] = na[k]-nb[k];
      if( !mpmThermal::normalize( normal ) ) continue;
      for( integer k = 0; k < 3; ++k )
        separation[k] = geometry[node][b][5+k]/geometry[node][b][0]-geometry[node][a][5+k]/geometry[node][a][0];
      if( mpmThermal::dot( normal, separation ) < 0.0 ) for( auto & x : normal ) x = -x;
      bool const surfaceA = geometry[node][a][11] > 0.0, surfaceB = geometry[node][b][11] > 0.0;
      real64 inverseLengthSquared = 0.0;
      for( integer k = 0; k < m_numDims; ++k ) inverseLengthSquared += normal[k]*normal[k]/( h[k]*h[k] );
      real64 const normalLength = 1.0/std::sqrt( inverseLengthSquared );
      for( integer k = 0; k < 3; ++k )
      {
        pa[k] = surfaceA ? geometry[node][a][8+k]/geometry[node][a][11] : geometry[node][a][5+k]/geometry[node][a][0]+0.5*normalLength*normal[k];
        pb[k] = surfaceB ? geometry[node][b][8+k]/geometry[node][b][11] : geometry[node][b][5+k]/geometry[node][b][0]-0.5*normalLength*normal[k];
        separation[k] = pb[k]-pa[k]; pointOnPlane[k] = 0.5*( pa[k]+pb[k] );
      }
      if( mpmThermal::dot( normal, separation ) > settings.gapTolerance*hMin ) continue;
      mpmThermal::Vector lower, upper;
      real64 const cubicA = geometry[node][a][12]/geometry[node][a][0];
      real64 const cubicB = geometry[node][b][12]/geometry[node][b][0];
      real64 const supportScale = cubicA+cubicB > 0.0 ? 2.0 : 1.0;
      for( integer k = 0; k < 3; ++k )
      {
        lower[k] = periodic[k] ? -supportScale*h[k] : std::max( -supportScale*h[k], m_xGlobalMin[k]-position[node][k] );
        upper[k] = periodic[k] ? supportScale*h[k] : std::min( supportScale*h[k], m_xGlobalMax[k]-position[node][k] );
      }
      real64 const area = mpmThermal::planeNodalArea( normal, pointOnPlane, h, lower, upper, settings.areaQuadrature, cubicA, cubicB );
      if( !( area > 0.0 ) ) continue;
      bool const perfect = mode == "Perfect" || !( resistance > 0.0 );
      real64 const conductance = perfect ? 0.0 : area/resistance;
      GEOS_ERROR_IF( !std::isfinite( conductance ), "Thermal conductance overflow; use Perfect contact for zero resistance." );
      interfaces[node].push_back( { a, b, conductance, perfect, false } );
    }
  }
  real64 elapsed = 0.0;
  integer substeps = 0;
  settings.suggestedTimeStep = std::numeric_limits< real64 >::max();
  while( elapsed < dt )
  {
    GEOS_ERROR_IF( substeps >= settings.maxSubsteps, "Thermal maxSubsteps exceeded; reduce the mechanical time step." );
    real64 const currentTime = time+elapsed;
    for( char const * key : { capacityKey, momentKey, fixedCapacityKey, fixedMomentKey, fixedSquareKey, rateKey, boundKey, incrementKey } )
      nodes.getReference< array2d< real64 > >( key ).zero();
    std::vector< real64 > power( point.size(), 0.0 ), target( point.size(), 0.0 );
    std::vector< integer > thermostat( point.size(), -1 );
    real64 localSourceLimit = std::numeric_limits< real64 >::max();
    for( std::size_t p = 0; p < point.size(); ++p )
    {
      auto & q = point[p];
      for( std::size_t c = 0; c < conditions.size(); ++c )
      {
        auto const & condition = *conditions[c];
        if( q.membership[c] == 0 || !condition.active( currentTime ) ) continue;
        if( condition.isothermal )
        {
          real64 const t = condition.targetTemperature( currentTime );
          GEOS_ERROR_IF( thermostat[p] >= 0 && std::abs( target[p]-t ) > 1.0e-10*std::max( 1.0, std::abs( t ) ),
                         "Conflicting isothermal conditions select the same particle." );
          thermostat[p] = static_cast< integer >( c ); target[p] = t;
        }
        else power[p] += condition.specificPower*q.mass + condition.volumetricPower*q.volume +
                         condition.power*q.mass/q.referenceMass[c];
      }
      if( thermostat[p] >= 0 )
      {
        GEOS_ERROR_IF( !finiteTemperature( target[p] ), "An isothermal condition prescribed a negative/nonfinite absolute temperature." );
        real64 const u = q.material->internalEnergy( target[p], q.mechanicalState );
        real64 const heat = u-*q.energy;
        *q.reservoir += heat; reservoirEnergy += q.mass*heat;
        *q.temperature = target[p]; *q.energy = u;
      }
      q.specificHeat = q.material->stateHeatCapacity( *q.temperature, q.mechanicalState );
      q.conductivity = q.material->thermalConductivity( *q.temperature );
      GEOS_ERROR_IF( !finiteTemperature( *q.temperature ) || !std::isfinite( q.specificHeat ) || !( q.specificHeat > 0.0 ) ||
                     !std::isfinite( q.conductivity ) || q.conductivity < 0.0 || !std::isfinite( power[p] ),
                     "Thermal properties left their admissible temperature range." );
      if( power[p] < 0.0 && thermostat[p] < 0 )
      {
        GEOS_ERROR_IF( !( *q.temperature > 0.0 ), "A heat sink cannot remove energy below zero absolute temperature." );
        localSourceLimit = std::min( localSourceLimit, 0.5*q.mass*q.specificHeat*( *q.temperature )/( -power[p] ) );
      }
      for( auto const & a : q.support )
      {
        real64 const c = a.weight*q.mass*q.specificHeat;
        capacity[a.node][a.slot] += c; moment[a.node][a.slot] += c*( *q.temperature );
        if( thermostat[p] >= 0 )
        {
          fixedCapacity[a.node][a.slot] += c; fixedMoment[a.node][a.slot] += c*target[p];
          fixedSquare[a.node][a.slot] += c*target[p]*target[p];
        }
      }
    }
    syncGridFields( { capacityKey, momentKey, fixedCapacityKey, fixedMomentKey, fixedSquareKey }, domain, nodes, mesh, MPI_SUM );
    for( localIndex node = 0; node < nodeCount; ++node )
    {
      mpmThermal::UnionFind bonded( slots );
      for( auto const & pair : interfaces[node] ) if( pair.bonded ) bonded.join( pair.a, pair.b );
      std::vector< real64 > c( slots, 0.0 ), t( slots, 0.0 );
      for( localIndex s = 0; s < slots; ++s ) { auto const r = bonded.root( s ); c[r] += capacity[node][s]; t[r] += moment[node][s]; }
      for( localIndex s = 0; s < slots; ++s )
      {
        auto const r = bonded.root( s ); temperature[node][s] = c[r] > 0.0 ? t[r]/c[r] : 0.0;
      }
    }
    // K_ij = sum_p V_p k_p grad(N_i).grad(N_j). The absolute row-sum
    // bound is conservative for the actual retained CPDI/B-spline stencil.
    for( auto const & q : point )
    {
      std::vector< mpmThermal::Vector > gradients;
      std::vector< real64 > nodalTemperature;
      for( auto const & a : q.support )
      {
        gradients.push_back( a.gradient ); nodalTemperature.push_back( temperature[a.node][a.slot] );
      }
      auto const assembly = mpmThermal::assembleFourierParticle( q.volume, q.conductivity, gradients, nodalTemperature );
      for( std::size_t i = 0; i < q.support.size(); ++i )
      {
        auto const & a = q.support[i];
        rate[a.node][a.slot] += assembly.power[i];
        bound[a.node][a.slot] += assembly.absoluteRowBound[i];
      }
    }
    syncGridFields( { rateKey, boundKey }, domain, nodes, mesh, MPI_SUM );
    real64 localLimit = localSourceLimit;
    for( localIndex node = 0; node < nodeCount; ++node ) if( ghostRank[node] < 0 )
      for( localIndex s = 0; s < slots; ++s ) if( capacity[node][s] > 0.0 && bound[node][s] > 0.0 )
        localLimit = std::min( localLimit, settings.safetyFactor*2.0*capacity[node][s]/bound[node][s] );
    real64 const stable = MpiWrapper::min( localLimit );
    settings.suggestedTimeStep = std::min( settings.suggestedTimeStep, stable );
    real64 subDt = std::min( dt-elapsed, stable );
    real64 const nextEvent = settings.nextConditionTime( currentTime );
    if( nextEvent < time+dt ) subDt = std::min( subDt, nextEvent-currentTime );
    GEOS_ERROR_IF( !( subDt > 0.0 ) || !std::isfinite( subDt ) || !( elapsed+subDt > elapsed ),
                   "Thermal substep cannot advance time." );
    // Node owners alone solve the graph. Never SUM a derived solution on ghosts.
    for( localIndex node = 0; node < nodeCount; ++node ) if( ghostRank[node] < 0 )
    {
      std::vector< localIndex > active, lookup( slots, -1 );
      std::vector< mpmThermal::ThermalNode > state;
      for( localIndex s = 0; s < slots; ++s ) if( capacity[node][s] > 0.0 )
      {
        lookup[s] = static_cast< localIndex >( active.size() ); active.push_back( s );
        mpmThermal::ThermalNode n;
        n.capacity = capacity[node][s]; n.temperature = temperature[node][s]; n.bulkEnergy = subDt*rate[node][s];
        if( fixedCapacity[node][s] >= ( 1.0-1.0e-12 )*n.capacity )
        {
          real64 const t = fixedMoment[node][s]/fixedCapacity[node][s];
          real64 const variance = fixedSquare[node][s]/fixedCapacity[node][s]-t*t;
          // Distinct thermostatted particles can share a slot. In that case the
          // particle constraints, not an invalid common nodal value, own the BC.
          if( variance <= 1.0e-10*std::max( 1.0, t*t ) ) { n.prescribed = true; n.prescribedTemperature = t; }
        }
        state.push_back( n );
      }
      std::vector< mpmThermal::ThermalEdge > edges;
      for( auto const & pair : interfaces[node] ) if( lookup[pair.a] >= 0 && lookup[pair.b] >= 0 )
        edges.push_back( { static_cast< std::size_t >( lookup[pair.a] ), static_cast< std::size_t >( lookup[pair.b] ), pair.conductance, pair.perfect } );
      auto const energy = mpmThermal::exchangeGraph( subDt, state, edges );
      for( std::size_t i = 0; i < active.size(); ++i ) increment[node][active[i]] = energy[i];
    }
    replaceGridFieldsOwnerToGhost( { incrementKey }, domain, nodes, mesh );
    // Frozen-capacity G2P is exactly conservative: sum_p m_p*c_p*N_ip/C_i=1.
    // Invert each particle's own caloric policy after adding the specific energy.
    for( std::size_t p = 0; p < point.size(); ++p )
    {
      auto & q = point[p]; real64 du = subDt*power[p]/q.mass;
      for( auto const & a : q.support ) du += mpmThermal::gatherSpecificEnergy( q.specificHeat, a.weight,
                                               increment[a.node][a.slot], capacity[a.node][a.slot] );
      *q.energy += du; sourceEnergy += subDt*power[p];
      if( thermostat[p] >= 0 )
      {
        real64 const t = conditions[thermostat[p]]->targetTemperature( currentTime+subDt );
        GEOS_ERROR_IF( !finiteTemperature( t ), "An isothermal temperature history crossed zero absolute temperature." );
        real64 const u = q.material->internalEnergy( t, q.mechanicalState );
        real64 const bath = u-*q.energy; *q.reservoir += bath; reservoirEnergy += q.mass*bath; *q.energy = u;
      }
      *q.isothermal = thermostat[p] >= 0 ? 1 : 0;
      *q.temperature = q.material->temperatureFromInternalEnergy( *q.energy, *q.temperature, q.mechanicalState );
      GEOS_ERROR_IF( !finiteTemperature( *q.temperature ), "Thermal update produced an inadmissible absolute temperature; reduce the time step or source strength." );
    }
    elapsed = subDt >= dt-elapsed ? dt : elapsed+subDt;
    ++substeps;
  }
  real64 finalEnergy = 0.0;
  for( auto & q : point )
  {
    *q.heatCapacity = q.material->stateHeatCapacity( *q.temperature, q.mechanicalState );
    real64 const conductivity = q.material->thermalConductivity( *q.temperature );
    GEOS_ERROR_IF( !std::isfinite( *q.heatCapacity ) || !( *q.heatCapacity > 0.0 ) ||
                   !std::isfinite( conductivity ) || conductivity < 0.0,
                   "Final thermal state is outside the material's admissible property range." );
    *q.temperatureRate = ( *q.temperature-q.initialTemperature )/dt;
    *q.strain = q.alpha*( *q.temperature-q.initialTemperature );
    // Mechanics consumes the COMPLETE increment, including locally generated
    // plastic heat. It owns committing previousTemperature and thermal stretch.
    if( m_rigidBodyMode == 1 ) *q.previousTemperature = *q.temperature;
    GEOS_ERROR_IF( !std::isfinite( *q.stretch ) || !( *q.stretch > 0.0 ), "Invalid thermal stretch." );
    finalEnergy += q.mass*( *q.energy );
  }
  real64 const globalSource = MpiWrapper::sum( sourceEnergy );
  real64 const globalReservoir = MpiWrapper::sum( reservoirEnergy );
  settings.cumulativeSourceEnergy += globalSource;
  settings.cumulativeReservoirEnergy += globalReservoir;
  if( settings.diagnostics != 0 )
  {
    real64 const change = MpiWrapper::sum( finalEnergy-initialEnergy );
    GEOS_LOG_RANK_0( "MPMThermal time=" << time+dt << " substeps=" << substeps << " energyChange=" << change
                    << " sourceEnergy=" << globalSource << " reservoirEnergy=" << globalReservoir
                    << " balanceResidual=" << change-globalSource-globalReservoir );
  }
  // Only particle fields needed by the next mechanical device kernels return to
  // device. Thermal masks/energy and temporary grids may remain host-resident.
  particles.forParticleSubRegions( [&]( ParticleSubRegion & sub )
  {
    for( char const * name : { "particleTemperature", "particleTemperatureRate", "particleHeatCapacity",
                              "particleThermalStrainIncrement", "particleThermalStretch" } )
      sub.getWrapperBase( name ).move( parallelDeviceMemorySpace, true );
  } );
}

} // namespace geos
