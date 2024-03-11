/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RegionMultiphaseStatistics.cpp
 */

#include "RegionMultiphaseStatistics.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
//#include "finiteVolume/FiniteVolumeManager.hpp"
//#include "finiteVolume/FluxApproximationBase.hpp"
//#include "mainInterface/ProblemManager.hpp"
//#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

// Container for region statistics
struct RegionMultiphaseStatistics::RegionStatistics
{
  /// Indices of minimum, maximum and total property values
  static constexpr integer MINIMUM = 0;
  static constexpr integer MAXIMUM = 1;
  static constexpr integer TOTAL = 2;
  static constexpr integer WEIGHTED_TOTAL = 3;
  static constexpr integer TYPE_END = 4;

  /// Indices for property type
  static constexpr integer PORE_VOLUME = 0;
  static constexpr integer PRESSURE = 1;
  static constexpr integer TEMPERATURE = 2;
  static constexpr integer PHASE_PORE_VOLUME = 3;
  static constexpr integer PHASE_MASS = 4;
  static constexpr integer PHASE_COMP_MASS = 5;
  static constexpr integer PROP_END = 6;

  static std::unique_ptr< array3d< real64 > > allocate( integer regionCount, integer phaseCount, integer componentCount )
  {
    std::unique_ptr< array3d< real64 > > properties = std::make_unique< array3d< real64 > >(
      TYPE_END,
      regionCount,
      getSize( phaseCount, componentCount ));
    LvArray::forValuesInSlice( properties->toSlice(), []( real64 & v ){ v = 0.0; } );

    static constexpr real64 minValue = -LvArray::NumericLimits< real64 >::max;
    static constexpr real64 maxValue =  LvArray::NumericLimits< real64 >::max;
    LvArray::forValuesInSlice( (*properties)[MINIMUM], []( real64 & v ){ v = maxValue; } );
    LvArray::forValuesInSlice( (*properties)[MAXIMUM], []( real64 & v ){ v = minValue; } );

    return std::move( properties );
  }

  static std::unique_ptr< array1d< localIndex > > allocateCounts( integer regionCount )
  {
    std::unique_ptr< array1d< localIndex > > elementCounts = std::make_unique< array1d< localIndex > >( regionCount );
    LvArray::forValuesInSlice( elementCounts->toSlice(), []( localIndex & v ){ v = 0; } );
    return std::move( elementCounts );
  }

  static real64 * getProperty( arrayView3d< real64 > props,
                               integer phaseCount,
                               integer componentCount,
                               integer region,
                               integer type,
                               integer property,
                               integer phase = 0,
                               integer component = 0 )
  {
    return &props( type, region, calculateIndex( phaseCount, componentCount, property, phase, component ));
  }

private:
  static integer getSize( integer phaseCount, integer componentCount )
  {
    // 1 pore volume, 1 pressure, 1 temperature, np phase volume, np phase mass, np*nc phase comp mass
    return 3 + 2*phaseCount + phaseCount*componentCount;
  }

  static integer calculateIndex( integer phaseCount,
                                 integer componentCount,
                                 integer property,
                                 integer phase,
                                 integer component )
  {
    if( property == PORE_VOLUME || property == PRESSURE || property == TEMPERATURE )
    {
      return property;
    }
    else if( property == PHASE_PORE_VOLUME )
    {
      return PHASE_PORE_VOLUME + phase;
    }
    else if( property == PHASE_MASS )
    {
      return PHASE_PORE_VOLUME + phaseCount;
    }
    else if( property == PHASE_COMP_MASS )
    {
      return PHASE_PORE_VOLUME + 2*phaseCount + phase*componentCount + component;
    }
    return 0;
  }
};

// The region statistics kernel
struct RegionMultiphaseStatistics::RegionStatisticsKernel
{
  template< typename POLICY >
  static void launch( localIndex const size,
                      integer const phaseCount,
                      integer const componentCount,
                      arrayView3d< real64 > const regionStatistics,
                      arrayView1d< localIndex > const & elementCounts,
                      arrayView1d< real64 const > const & regionIndices,
                      arrayView1d< real64 const > const & regionMarkers,
                      arrayView1d< integer const > const & elementGhostRank,
                      arrayView1d< real64 const > const & elementVolume,
                      arrayView1d< real64 const > const & pressure,
                      arrayView1d< real64 const > const & temperature,
                      arrayView1d< real64 const > const & referencePorosity,
                      arrayView2d< real64 const > const & porosity,
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDensity,
                      arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseComponentFraction,
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolumeFraction,
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseTrappedVolumeFraction,
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelativePermeability )
  {
    GEOS_UNUSED_VAR( regionStatistics );
    GEOS_UNUSED_VAR( regionMarkers );
    GEOS_UNUSED_VAR( elementVolume );
    GEOS_UNUSED_VAR( pressure );
    GEOS_UNUSED_VAR( temperature );
    GEOS_UNUSED_VAR( referencePorosity );
    GEOS_UNUSED_VAR( porosity );
    GEOS_UNUSED_VAR( phaseDensity );
    GEOS_UNUSED_VAR( phaseComponentFraction );
    GEOS_UNUSED_VAR( phaseVolumeFraction );
    GEOS_UNUSED_VAR( phaseTrappedVolumeFraction );
    GEOS_UNUSED_VAR( phaseRelativePermeability );

    integer const regionCount = regionIndices.size();

    auto populateMinMaxTotal = [regionStatistics,
                                phaseCount,
                                componentCount] GEOS_HOST_DEVICE ( real64 const value,
                                                                   real64 const weight,
                                                                   integer region,
                                                                   integer property,
                                                                   integer phase = 0,
                                                                   integer component = 0 )
    {
      real64 * targetValue = nullptr;
      targetValue = RegionStatistics::getProperty( regionStatistics, region, phaseCount, componentCount, RegionStatistics::MINIMUM, property, phase, component );
      RAJA::atomicMin< AtomicPolicy< POLICY > >( targetValue, value );
      targetValue = RegionStatistics::getProperty( regionStatistics, region, phaseCount, componentCount, RegionStatistics::MAXIMUM, property, phase, component );
      RAJA::atomicMax< AtomicPolicy< POLICY > >( targetValue, value );
      targetValue = RegionStatistics::getProperty( regionStatistics, region, phaseCount, componentCount, RegionStatistics::TOTAL, property, phase, component );
      RAJA::atomicAdd< AtomicPolicy< POLICY > >( targetValue, value );
      targetValue = RegionStatistics::getProperty( regionStatistics, region, phaseCount, componentCount, RegionStatistics::WEIGHTED_TOTAL, property, phase, component );
      RAJA::atomicAdd< AtomicPolicy< POLICY > >( targetValue, weight * value );
      return value;
    };

    forAll< POLICY >( size, [phaseCount,
                             componentCount,
                             regionCount,
                             regionStatistics,
                             regionIndices,
                             regionMarkers,
                             elementCounts,
                             elementGhostRank,
                             elementVolume,
                             pressure,
                             temperature,
                             referencePorosity,
                             porosity,
                             phaseDensity,
                             phaseComponentFraction,
                             phaseVolumeFraction,
                             phaseTrappedVolumeFraction,
                             phaseRelativePermeability,
                             &populateMinMaxTotal] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( elementGhostRank[ei] >= 0 )
      {
        return;
      }

      // Find the appropriate region
      integer region = 0;
      real64 const regionMarker = regionMarkers[ei];
      for(; region < regionCount; ++region )
      {
        if( LvArray::math::abs( regionIndices[region] - regionMarker ) < 0.1 )
        {
          break;
        }
      }

      if( regionCount <= region )
      {
        return;
      }

      RAJA::atomicAdd< AtomicPolicy< POLICY > >( &elementCounts[region], 1 );

      // Uncompacted pore volume which is also used as the weight
      real64 const weight = elementVolume[ei] * referencePorosity[ei];

      populateMinMaxTotal( weight, 1.0, region, RegionStatistics::PORE_VOLUME );

      populateMinMaxTotal( pressure[ei], weight, region, RegionStatistics::PRESSURE );
      populateMinMaxTotal( temperature[ei], weight, region, RegionStatistics::TEMPERATURE );

      for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
      {
        real64 const elementPhaseVolume = elementVolume[ei] * porosity[ei][0] * phaseVolumeFraction[ei][phaseIndex];
        real64 const elementPhaseDensity = phaseDensity[ei][0][phaseIndex];
        real64 const elementPhaseMass = elementPhaseVolume * elementPhaseDensity;

        populateMinMaxTotal( elementPhaseVolume, weight, region, RegionStatistics::PHASE_PORE_VOLUME, phaseIndex );
        populateMinMaxTotal( elementPhaseMass, weight, region, RegionStatistics::PHASE_MASS, phaseIndex );

        for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
        {
          real64 const elementComponentPhaseMass = elementPhaseMass * phaseComponentFraction[ei][0][phaseIndex][compIndex];
          populateMinMaxTotal( elementComponentPhaseMass, weight, region, RegionStatistics::PHASE_COMP_MASS, phaseIndex, compIndex );
        }
      }
    } );

    // Dummy loop to bring data back to the CPU
    forAll< serialPolicy >( 1, [regionStatistics, elementCounts] ( localIndex const )
    {
      GEOS_UNUSED_VAR( regionStatistics, elementCounts );
    } );
  }
};

RegionMultiphaseStatistics::RegionMultiphaseStatistics( const string & name,
                                                        Group * const parent ):
  Base( name, parent )
{
  registerWrapper( viewKeyStruct::regionNamesString(), &m_regionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The list of names of the regions" );

  registerWrapper( viewKeyStruct::regionIdentifiersString(), &m_regionIdentifiers ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The list of integer identifiers of the regions" );

  registerWrapper( viewKeyStruct::propertyNamesString(), &m_propertyNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The list of names of properties to report" );
}

void RegionMultiphaseStatistics::postProcessInput()
{
  Base::postProcessInput();

  // Make sure we have the same number of region names and region indices
  GEOS_ERROR_IF_NE_MSG( m_regionNames.size(), m_regionIdentifiers.size(),
                        GEOS_FMT( "The number of region names ({}) is different from the number of region identifiers ({}).",
                                  m_regionNames.size(), m_regionIdentifiers.size() ) );

  // Convert the property names to integers
  m_propertyNameTypes.clear();
  localIndex const propSize = m_propertyNames.size();
  if( propSize == 0 )
  {
    m_propertyNameTypes.emplace_back( PropertyNameType::Pressure );
    m_propertyNameTypes.emplace_back( PropertyNameType::PoreVolume );
    m_propertyNameTypes.emplace_back( PropertyNameType::Saturation );
  }
  else
  {
    m_propertyNameTypes.resize( propSize );
    for( integer i = 0; i < propSize; ++i )
    {
      m_propertyNameTypes[i] = EnumStrings< PropertyNameType >::fromString( m_propertyNames[i] );
    }
  }
}

void RegionMultiphaseStatistics::registerDataOnMesh( Group & meshBodies )
{
  Base::registerDataOnMesh( meshBodies );

  // The fields have to be registered in "registerDataOnMesh" (and not later)
  // otherwise they cannot be targeted by TimeHistory

  // For now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr )
  {
    return;
  }

  string const fieldName = GEOS_FMT( "{}{}", getName(), viewKeyStruct::fieldNameString() );

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions( regionNames, [&]( localIndex const,
                                                        ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( fieldName ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setDescription( fieldName ).
        setRegisteringObjects( getName() );
    } );
  } );
}

bool RegionMultiphaseStatistics::execute( real64 const time_n,
                                          real64 const dt,
                                          integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                          integer const GEOS_UNUSED_PARAM( eventCounter ),
                                          real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                          DomainPartition & domain )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          arrayView1d< string const > const & regionNames )
  {
    real64 const currentTime = time_n + dt;
    computeRegionStatistics( currentTime, mesh, regionNames );
  } );

  return false;
}

void RegionMultiphaseStatistics::computeRegionStatistics( real64 const time,
                                                          MeshLevel & mesh,
                                                          arrayView1d< string const > const & regionNames ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time );

  integer const numPhases = m_solver->numFluidPhases();
  integer const numComps = m_solver->numFluidComponents();

  string const fieldName = GEOS_FMT( "{}{}", getName(), viewKeyStruct::fieldNameString() );

  // Step 1: initialize the average/min/max quantities
  // All properties will be stored in one large array
  std::unique_ptr< array3d< real64 > > regionStatistics = RegionStatistics::allocate( m_regionNames.size(), numPhases, numComps );
  std::unique_ptr< array1d< localIndex > > elementCounts = RegionStatistics::allocateCounts( m_regionNames.size() );

  // Step 2: increment the average/min/max quantities for all the subRegions
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const elementGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const elementVolume = subRegion.getElementVolume();
    array1d< real64 > const regionMarkers = subRegion.getWrapper< array1d< real64 > >( fieldName ).reference();
    arrayView1d< real64 const > const pressure = subRegion.getField< fields::flow::pressure >();
    arrayView1d< real64 const > const temperature = subRegion.getField< fields::flow::temperature >();
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolumeFraction = subRegion.getField< fields::flow::phaseVolumeFraction >();

    Group const & constitutiveModels = subRegion.getGroup( ElementSubRegionBase::groupKeyStruct::constitutiveModelsString() );

    string const & solidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::solidNamesString() );
    CoupledSolidBase const & solid = constitutiveModels.getGroup< CoupledSolidBase >( solidName );
    arrayView1d< real64 const > const referencePorosity = solid.getReferencePorosity();
    arrayView2d< real64 const > const porosity = solid.getPorosity();

    string const & fluidName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::fluidNamesString() );
    MultiFluidBase const & fluid = constitutiveModels.getGroup< MultiFluidBase >( fluidName );
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDensity = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseComponentFraction = fluid.phaseCompFraction();

    string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseTrappedVolumeFraction = relperm.phaseTrappedVolFraction();
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseRelativePermeability = relperm.phaseRelPerm();

    RegionStatisticsKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                              numPhases,
                                                              numComps,
                                                              regionStatistics->toView(),
                                                              elementCounts->toView(),
                                                              m_regionIdentifiers.toViewConst(),
                                                              regionMarkers.toViewConst(),
                                                              elementGhostRank,
                                                              elementVolume,
                                                              pressure,
                                                              temperature,
                                                              referencePorosity,
                                                              porosity,
                                                              phaseDensity,
                                                              phaseComponentFraction,
                                                              phaseVolumeFraction,
                                                              phaseTrappedVolumeFraction,
                                                              phaseRelativePermeability );
  } );

  // Step 3: synchronize the results over the MPI ranks
  for( integer regionIndex = 0; regionIndex < m_regionIdentifiers.size(); ++regionIndex )
  {}

  for( integer i = 0; i < elementCounts->size(); ++i )
  {
    std::cout << "Region " << i + 1 << ": " << (*elementCounts)[i] << "\n";
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        RegionMultiphaseStatistics,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
