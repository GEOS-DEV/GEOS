/* SPDX-License-Identifier: LGPL-2.1-only */
#include "HydrostaticStateValidation.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "common/MpiWrapper.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{
namespace
{
char const * const recordKey = "__hydrostaticAcceptedPhasePressure";
}

bool hasHydrostaticStateRecords( ElementSubRegionBase const & subRegion )
{
  return subRegion.hasWrapper( recordKey );
}

void recordHydrostaticState( ElementSubRegionBase & subRegion, localIndex element,
                             integer numPhases, HydrostaticPhaseVector const & phasePressure, real64 tolerance )
{
  if( !hasHydrostaticStateRecords( subRegion ) )
  {
    auto & wrapper = subRegion.registerWrapper< array2d< real64 > >( recordKey );
    wrapper.setInputFlag( dataRepository::InputFlags::FALSE ).setRestartFlags( dataRepository::RestartFlags::NO_WRITE )
      .setPlotLevel( dataRepository::PlotLevel::NOPLOT );
    wrapper.reference().resize( subRegion.size(), numPhases + 1 );
    wrapper.reference().setValues< serialPolicy >( 0.0 );
  }
  auto & records = subRegion.getReference< array2d< real64 > >( recordKey );
  records.move( hostMemorySpace, true );
  GEOS_THROW_IF( records.size( 1 ) != numPhases + 1, "Inconsistent hydrostatic validation phase count", InputError );
  records[element][0] = tolerance;
  for( integer ip = 0; ip < numPhases; ++ip ) records[element][ip+1] = phasePressure[ip];
}

void realizeSelectedHydrostaticState( ElementSubRegionBase & subRegion, std::function< void() > const & update )
{
  auto const records = subRegion.getReference< array2d< real64 > >( recordKey ).toViewConst();
  records.move( hostMemorySpace, false );
  array1d< localIndex > outside;
  for( localIndex k = 0; k < subRegion.size(); ++k ) if( records[k][0] <= 0.0 ) outside.emplace_back( k );
  if( outside.empty() ) { update(); return; }
  auto const indices = outside.toViewConst();
  struct Saved
  {
    dataRepository::WrapperBase * wrapper;
    std::vector< buffer_unit_type > bytes;
  };
  std::vector< Saved > saved;
  parallelDeviceEvents events;
  auto preserve = [&]( dataRepository::WrapperBase & wrapper )
  {
    wrapper.move( hostMemorySpace, false );
    buffer_unit_type * position = nullptr;
    localIndex const count = wrapper.packByIndex< false >( position, indices, false, false, events );
    saved.push_back( { &wrapper, std::vector< buffer_unit_type >( count ) } );
    position = saved.back().bytes.data();
    wrapper.packByIndex< true >( position, indices, false, false, events );
  };
  using namespace fields;
  // These are precisely the solver-owned fields written by updateFluidState.
  for( char const * key : { flow::globalCompFraction::key(), flow::dGlobalCompFraction_dGlobalCompDensity::key(),
                           flow::compAmount::key(), flow::phaseVolumeFraction::key(), flow::dPhaseVolumeFraction::key(),
                           flow::phaseMobility::key(), flow::dPhaseMobility::key() } )
    preserve( subRegion.getWrapperBase( key ) );
  subRegion.getConstitutiveModels().forSubGroups< constitutive::MultiFluidBase,
                                                  constitutive::RelativePermeabilityBase,
                                                  constitutive::CapillaryPressureBase >( [&]( auto & model )
  {
    model.forWrappers( [&]( dataRepository::WrapperBase & wrapper )
    {
      if( wrapper.sizedFromParent() != 0 && wrapper.isPackable( false ) ) preserve( wrapper );
    } );
  } );
  auto restore = [&]()
  {
    for( auto & entry : saved )
    {
      entry.wrapper->move( hostMemorySpace, true );
      buffer_unit_type const * position = entry.bytes.data();
      entry.wrapper->unpackByIndex( position, indices, false, false, events );
    }
  };
  try { update(); }
  catch( ... ) { restore(); throw; }
  restore();
}

void validateRealizedHydrostaticState( ElementSubRegionBase & subRegion,
                                      constitutive::CapillaryPressureBase const & capillary )
{
  integer invalid = 0;
  real64 maximumResidual = 0.0;
  if( hasHydrostaticStateRecords( subRegion ) )
  {
    auto const records = subRegion.getReference< array2d< real64 > >( recordKey ).toViewConst();
    arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
    arrayView2d< real64 const, compflow::USD_PHASE > const mobility = subRegion.getField< fields::flow::phaseMobility >();
    auto const capillaryPressure = capillary.phaseCapPressure();
    auto const ghosts = subRegion.ghostRank();
    records.move( hostMemorySpace, false ); pressure.move( hostMemorySpace, false );
    mobility.move( hostMemorySpace, false ); capillaryPressure.move( hostMemorySpace, false ); ghosts.move( hostMemorySpace, false );
    integer const numPhases = records.size( 1 ) - 1;
    for( localIndex k = 0; k < subRegion.size(); ++k )
    {
      if( records[k][0] <= 0.0 || ghosts[k] >= 0 ) continue;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        real64 const expected = records[k][ip+1];
        real64 const actual = pressure[k] - capillaryPressure[k][0][ip];
        real64 const residual = expected - actual;
        real64 const scale = std::max( 1.0, std::max( std::abs( expected ), std::abs( pressure[k] ) ) );
        real64 const allowance = records[k][0] + 128.0 * std::numeric_limits< real64 >::epsilon() * scale;
        real64 const error = mobility[k][ip] > 0.0 ? std::abs( residual ) : std::max( 0.0, residual );
        if( !std::isfinite( actual ) || !std::isfinite( mobility[k][ip] ) || mobility[k][ip] < 0.0 ) invalid = 1;
        if( error > allowance ) { invalid = 1; maximumResidual = std::max( maximumResidual, error ); }
      }
    }
    subRegion.deregisterWrapper( recordKey );
  }
  // Empty local selections still participate, avoiding rank-dependent checks.
  integer const globalInvalid = MpiWrapper::max( invalid );
  real64 const globalResidual = MpiWrapper::max( maximumResidual );
  GEOS_THROW_IF( globalInvalid != 0,
                 GEOS_FMT( "Realized hydrostatic EOS/capillary state violates mobile-phase pressure or endpoint complementarity after the normal component-density/EOS update (maximum pressure residual {} Pa). Check residual endpoints and relative permeability; allowLocalCompDensityChopping can introduce mobile phases, so set it to 0 when exact phase absence is required. No pressure tolerance was relaxed.", globalResidual ),
                 InputError, subRegion.getDataContext() );
}
}
}
