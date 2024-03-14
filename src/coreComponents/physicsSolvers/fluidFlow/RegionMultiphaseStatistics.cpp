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
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "mesh/MeshFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"

namespace geos
{

struct RegionMultiphaseStatistics::Statistics
{
  static constexpr integer MAX_NUM_PHASE = 3;
  static constexpr integer MAX_NUM_COMPS = 5;
  /// Indices for properties
  static constexpr integer STATIC_PORE_VOLUME = 0;
  static constexpr integer PORE_VOLUME = 1;

  static constexpr integer MINIMUM_PRESSURE = 2;
  static constexpr integer MAXIMUM_PRESSURE = 3;
  static constexpr integer AVERAGE_PRESSURE = 4;
  static constexpr integer MINIMUM_TEMPERATURE = 5;
  static constexpr integer MAXIMUM_TEMPERATURE = 6;
  static constexpr integer AVERAGE_TEMPERATURE = 7;
  // These are starting positions
  static constexpr integer PHASE_PORE_VOLUME = 8;
  static constexpr integer PHASE_DENSITY = PHASE_PORE_VOLUME + MAX_NUM_PHASE;
  static constexpr integer PHASE_VISCOSITY = PHASE_DENSITY + MAX_NUM_PHASE;
  static constexpr integer PHASE_MASS = PHASE_VISCOSITY + MAX_NUM_PHASE;
  static constexpr integer PHASE_COMP_MASS = PHASE_MASS + MAX_NUM_PHASE;
  // End marker
  static constexpr integer END = PHASE_COMP_MASS + MAX_NUM_COMPS*MAX_NUM_PHASE;

  GEOS_HOST_DEVICE Statistics( real64 const val )
  {
    for( integer i = 0; i < END; i++ )
    {
      m_data[i] = val;
    }
  }

  GEOS_HOST_DEVICE Statistics( Statistics const & rhs )
  {
    for( integer i = 0; i < END; i++ )
    {
      m_data[i] = rhs.m_data[i];
    }
  }

  GEOS_HOST_DEVICE Statistics const & operator=( Statistics const & rhs )
  {
    for( integer i = 0; i < END; i++ )
    {
      m_data[i] = rhs.m_data[i];
    }
    return *this;
  }

  Statistics( Statistics && ) = delete;
  ~Statistics() = default;

  GEOS_HOST_DEVICE bool operator!=( Statistics const & ) const
  {
    return false;
  }

  real64 m_data[END];
};

}

namespace RAJA
{
namespace operators
{
template< >
struct plus< geos::RegionMultiphaseStatistics::Statistics >
{
  using DataType = geos::RegionMultiphaseStatistics::Statistics;
  static constexpr int END = geos::RegionMultiphaseStatistics::Statistics::END;

  GEOS_HOST_DEVICE DataType operator()( const DataType & lhs, const DataType & rhs ) const
  {
    DataType res( lhs );
    for( int i = 0; i < END; ++i )
    {
      res.m_data[i] += rhs.m_data[i];
    }
    return res;
  }

  GEOS_HOST_DEVICE static DataType identity()
  {
    return DataType( 0.0 );
  }
};

template<>
struct minimum< geos::RegionMultiphaseStatistics::Statistics >
{
  using DataType = geos::RegionMultiphaseStatistics::Statistics;
  static constexpr int END = geos::RegionMultiphaseStatistics::Statistics::END;

  GEOS_HOST_DEVICE DataType operator()( const DataType & lhs, const DataType & rhs ) const
  {
    DataType res( lhs );
    for( int i = 0; i < END; ++i )
    {
      res.m_data[i] = LvArray::math::min( res.m_data[i], rhs.m_data[i] );
    }
    return res;
  }

  GEOS_HOST_DEVICE static DataType identity()
  {
    return DataType( LvArray::NumericLimits< geos::real64 >::max );
  }
};

}
}

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

// The region statistics kernel
struct RegionMultiphaseStatistics::RegionStatisticsKernel
{
  template< typename POLICY >
  static void launch( localIndex const size,
                      integer const phaseCount,
                      integer const componentCount,
                      integer const region,
                      arrayView2d< real64 > const regionStatistics,
                      arrayView1d< localIndex const > const & targetSet,
                      arrayView1d< real64 const > const & elementVolume,
                      arrayView1d< real64 const > const & pressure,
                      arrayView1d< real64 const > const & temperature,
                      arrayView1d< real64 const > const & referencePorosity,
                      arrayView2d< real64 const > const & porosity,
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDensity,
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseViscosity,
                      arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseComponentFraction,
                      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolumeFraction,
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseTrappedVolumeFraction,
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelativePermeability )
  {

    GEOS_UNUSED_VAR( region );

    RAJA::ReduceSum< ReducePolicy< POLICY >, Statistics > sumStats( RAJA::operators::plus< Statistics >::identity() );
    RAJA::ReduceMin< ReducePolicy< POLICY >, Statistics > minStats( RAJA::operators::minimum< Statistics >::identity() );

    forAll< POLICY >( size, [phaseCount,
                             componentCount,
                             targetSet,
                             elementVolume,
                             pressure,
                             temperature,
                             referencePorosity,
                             porosity,
                             phaseDensity,
                             phaseViscosity,
                             phaseComponentFraction,
                             phaseVolumeFraction,
                             phaseTrappedVolumeFraction,
                             phaseRelativePermeability,
                             sumStats,
                             minStats] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const ei = targetSet[i];

      RegionMultiphaseStatistics::Statistics stats( RAJA::operators::plus< Statistics >::identity() );

      // Uncompacted pore volume which is also used as the weight
      real64 const staticPoreVolume = elementVolume[ei] * referencePorosity[ei];
      real64 const poreVolume = elementVolume[ei] * porosity[ei][0];

      stats.m_data[Statistics::STATIC_PORE_VOLUME] = staticPoreVolume;
      stats.m_data[Statistics::PORE_VOLUME] = poreVolume;

      stats.m_data[Statistics::MINIMUM_PRESSURE] = pressure[ei];
      stats.m_data[Statistics::MAXIMUM_PRESSURE] = -pressure[ei];
      stats.m_data[Statistics::AVERAGE_PRESSURE] = staticPoreVolume * pressure[ei];

      stats.m_data[Statistics::MINIMUM_TEMPERATURE] = temperature[ei];
      stats.m_data[Statistics::MAXIMUM_TEMPERATURE] = -temperature[ei];
      stats.m_data[Statistics::AVERAGE_TEMPERATURE] = staticPoreVolume * temperature[ei];

      for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
      {
        real64 const elementPhaseVolume = poreVolume * phaseVolumeFraction[ei][phaseIndex];
        real64 const elementPhaseDensity = phaseDensity[ei][0][phaseIndex];
        real64 const elementPhaseViscosity = phaseViscosity[ei][0][phaseIndex];
        real64 const elementPhaseMass = elementPhaseVolume * elementPhaseDensity;

        stats.m_data[Statistics::PHASE_PORE_VOLUME + phaseIndex] = elementPhaseVolume;
        stats.m_data[Statistics::PHASE_MASS + phaseIndex] = elementPhaseMass;

        stats.m_data[Statistics::PHASE_DENSITY + phaseIndex] = elementPhaseVolume * elementPhaseDensity;
        stats.m_data[Statistics::PHASE_VISCOSITY + phaseIndex] = elementPhaseVolume * elementPhaseViscosity;

        integer const phaseMassIndex = Statistics::PHASE_COMP_MASS + phaseIndex*Statistics::MAX_NUM_COMPS;
        for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
        {
          real64 const elementComponentPhaseMass = elementPhaseMass * phaseComponentFraction[ei][0][phaseIndex][compIndex];
          stats.m_data[phaseMassIndex + compIndex] = elementComponentPhaseMass;
        }
      }

      // Atomic operations
      sumStats += stats;
      minStats.min( stats );
    } );

    auto const sumData = sumStats.get().m_data;
    auto const minData = minStats.get().m_data;

    arraySlice1d< real64 > regionData = regionStatistics[region];

    regionData[Statistics::MINIMUM_PRESSURE] = LvArray::math::min( regionData[Statistics::MINIMUM_PRESSURE], minData[Statistics::MINIMUM_PRESSURE] );
    regionData[Statistics::MAXIMUM_PRESSURE] = LvArray::math::max( regionData[Statistics::MAXIMUM_PRESSURE], -minData[Statistics::MAXIMUM_PRESSURE] );
    regionData[Statistics::MINIMUM_TEMPERATURE] = LvArray::math::min( regionData[Statistics::MINIMUM_TEMPERATURE], minData[Statistics::MINIMUM_TEMPERATURE] );
    regionData[Statistics::MAXIMUM_TEMPERATURE] = LvArray::math::max( regionData[Statistics::MAXIMUM_TEMPERATURE], -minData[Statistics::MAXIMUM_TEMPERATURE] );

    regionData[Statistics::STATIC_PORE_VOLUME] += sumData[Statistics::STATIC_PORE_VOLUME];
    regionData[Statistics::PORE_VOLUME] += sumData[Statistics::PORE_VOLUME];
    regionData[Statistics::AVERAGE_PRESSURE] += sumData[Statistics::AVERAGE_PRESSURE];
    regionData[Statistics::AVERAGE_TEMPERATURE] += sumData[Statistics::AVERAGE_TEMPERATURE];
    for( integer propIndex = Statistics::PHASE_PORE_VOLUME; propIndex < Statistics::END; ++propIndex )
    {
      regionData[propIndex] += sumData[propIndex];
    }
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
    m_propertyNameTypes.emplace_back( PropertyNameType::Temperature );
    m_propertyNameTypes.emplace_back( PropertyNameType::StaticPoreVolume );
    m_propertyNameTypes.emplace_back( PropertyNameType::PoreVolume );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhaseVolumeFraction );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhasePoreVolume );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhaseMass );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhaseDensity );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhaseViscosity );
    m_propertyNameTypes.emplace_back( PropertyNameType::PhaseComponentMass );
    m_propertyNameTypes.emplace_back( PropertyNameType::ComponentMass );
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

template< typename LAMBDA >
void RegionMultiphaseStatistics::forRegions( Group & meshBodies, LAMBDA && lambda )
{
  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const regionIndex,
                                                                  ElementSubRegionBase & subRegion )
    {
      string const regionName = regionNames[regionIndex];
      lambda( regionIndex, regionName, subRegion );
    } );
  } );
}

template< typename ARRAY >
std::ostream & RegionMultiphaseStatistics::writeArray( std::ostream & os, ARRAY const & array ) const
{
  integer size = array.size();
  if( 0 < size )
  {
    os << array[0];
    for( integer i = 1; i < array.size(); ++i )
    {
      os << "," << array[i];
    }
    os << "\n";
  }
  return os;
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

  forRegions( meshBodies, [&]( localIndex const,
                               string const & regionName,
                               ElementSubRegionBase & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( fieldName ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setDescription( fieldName ).
      setRegisteringObjects( getName() );

    for( string const name : m_regionNames )
    {
      string const wrapperName = GEOS_FMT( "{}_{}_{}_{}", getName(), name, regionName, subRegion.getName() );

      registerWrapper< array1d< localIndex > >( wrapperName ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setDescription( wrapperName ).
        setSizedFromParent( 0 );
    }
  } );

  // Create the csv file if requires
  if( 0 < getLogLevel() && 0 < m_writeCSV && MpiWrapper::commRank() == 0 )
  {
    initializeFile();
  }
}

void RegionMultiphaseStatistics::initializePostInitialConditionsPreSubGroups()
{
  string const fieldName = GEOS_FMT( "{}{}", getName(), viewKeyStruct::fieldNameString() );

  integer const regionCount = m_regionNames.size();
  array1d< integer > regionIdentifiers( regionCount );
  forAll< serialPolicy >( m_regionIdentifiers.size(), [&] ( localIndex const ei )
  {
    regionIdentifiers[ei] = LvArray::math::convert< integer >( m_regionIdentifiers[ei] + 0.5 );
  } );

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  forRegions( domain.getMeshBodies(), [&]( localIndex const,
                                           string const & solverRegionName,
                                           ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const regionMarkers = subRegion.getWrapper< array1d< real64 > >( fieldName )
                                                        .reference().toViewConst();
    arrayView1d< integer const > const elementGhostRank = subRegion.ghostRank();

    array1d< localIndex > regionElementCounts( regionCount );
    regionElementCounts.zero();

    array2d< localIndex > regionElements( regionCount, regionMarkers.size() );

    forAll< serialPolicy >( regionMarkers.size(), [&] ( localIndex const ei )
    {
      if( 0 <= elementGhostRank[ei] )
      {
        return;
      }
      integer const marker = LvArray::math::convert< integer >( regionMarkers[ei] + 0.5 );
      for( integer region = 0; region < regionCount; region++ )
      {
        if( regionIdentifiers[region] == marker )
        {
          regionElements( region, regionElementCounts[region]++ ) = ei;
          return;
        }
      }
    } );

    arrayView2d< localIndex const > const regionElementsView = regionElements.toViewConst();

    for( integer region = 0; region < regionCount; region++ )
    {
      string const wrapperName = GEOS_FMT( "{}_{}_{}_{}", getName(), m_regionNames[region], solverRegionName, subRegion.getName() );
      auto & regionWrapper = getWrapper< array1d< localIndex > >( wrapperName );
      regionWrapper.resize( regionElementCounts[region] );

      arrayView1d< localIndex > const regionArrayView = regionWrapper.reference().toView();

      forAll< parallelDevicePolicy<> >( regionElementCounts[region],
                                        [region,
                                         regionArrayView,
                                         regionElementsView]
                                        GEOS_HOST_DEVICE ( localIndex const ei )
      {
        regionArrayView[ei] = regionElementsView( region, ei );
      } );
    }
  } );
}

bool RegionMultiphaseStatistics::execute( real64 const time_n,
                                          real64 const dt,
                                          integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                          integer const GEOS_UNUSED_PARAM( eventCounter ),
                                          real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                          DomainPartition & domain )
{
  // If the logLevel is low and we don't need the csv then there's nothing to do
  if( getLogLevel() <= 0 && m_writeCSV == 0 )
  {
    return false;
  }

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

  integer const phaseCount = m_solver->numFluidPhases();
  integer const componentCount = m_solver->numFluidComponents();
  integer const regionCount = m_regionNames.size();

  string const fieldName = GEOS_FMT( "{}{}", getName(), viewKeyStruct::fieldNameString() );

  // Step 1: initialize the average/min/max quantities
  // All properties will be stored in one large array
  array2d< real64 > regionStatisticsData( regionCount, Statistics::END );
  regionStatisticsData.zero();
  arrayView2d< real64 > regionStatistics = regionStatisticsData.toView();
  for( integer region = 0; region < regionCount; ++region )
  {
    regionStatistics( region, Statistics::MINIMUM_PRESSURE ) =  LvArray::NumericLimits< real64 >::max;
    regionStatistics( region, Statistics::MAXIMUM_PRESSURE ) = -LvArray::NumericLimits< real64 >::max;
    regionStatistics( region, Statistics::MINIMUM_TEMPERATURE ) =  LvArray::NumericLimits< real64 >::max;
    regionStatistics( region, Statistics::MAXIMUM_TEMPERATURE ) = -LvArray::NumericLimits< real64 >::max;
  }

  // Step 2: increment the average/min/max quantities for all the subRegions
  mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const solverRegionIndex,
                                                                ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const elementVolume = subRegion.getElementVolume();
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
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseViscosity = fluid.phaseViscosity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseComponentFraction = fluid.phaseCompFraction();

    string const & relpermName = subRegion.getReference< string >( CompositionalMultiphaseBase::viewKeyStruct::relPermNamesString() );
    RelativePermeabilityBase const & relperm = constitutiveModels.getGroup< RelativePermeabilityBase >( relpermName );
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseTrappedVolumeFraction = relperm.phaseTrappedVolFraction();
    arrayView3d< real64 const, relperm::USD_RELPERM > const phaseRelativePermeability = relperm.phaseRelPerm();

    for( integer region = 0; region < regionCount; ++region )
    {
      string const wrapperName = GEOS_FMT( "{}_{}_{}_{}", getName(), m_regionNames[region], regionNames[solverRegionIndex], subRegion.getName() );
      arrayView1d< localIndex const > const targetSet = getWrapper< array1d< localIndex > >( wrapperName )
                                                          .reference()
                                                          .toViewConst();
      if( targetSet.size() == 0 )
      {
        continue;
      }

      RegionStatisticsKernel::launch< parallelDevicePolicy<> >( targetSet.size(),
                                                                phaseCount,
                                                                componentCount,
                                                                region,
                                                                regionStatistics,
                                                                targetSet,
                                                                elementVolume,
                                                                pressure,
                                                                temperature,
                                                                referencePorosity,
                                                                porosity,
                                                                phaseDensity,
                                                                phaseViscosity,
                                                                phaseComponentFraction,
                                                                phaseVolumeFraction,
                                                                phaseTrappedVolumeFraction,
                                                                phaseRelativePermeability );
    }

  } );

  // Step 3: Force data back to cpu if necessary
  forAll< serialPolicy >( 1, [regionStatistics]( localIndex const ){
    GEOS_UNUSED_VAR( regionStatistics );
  } );

  // Step 4: synchronize the results over the MPI ranks
  array1d< real64 > phaseVolumes( phaseCount );

  for( integer region =0; region<regionCount; region++ )
  {
    arraySlice1d< real64 > regData = regionStatistics[region];

    regData[Statistics::STATIC_PORE_VOLUME] = MpiWrapper::sum( regData[Statistics::STATIC_PORE_VOLUME] );
    regData[Statistics::PORE_VOLUME] = MpiWrapper::sum( regData[Statistics::PORE_VOLUME] );

    regData[Statistics::MINIMUM_PRESSURE] = MpiWrapper::min( regData[Statistics::MINIMUM_PRESSURE] );
    regData[Statistics::MAXIMUM_PRESSURE] = MpiWrapper::max( regData[Statistics::MAXIMUM_PRESSURE] );
    regData[Statistics::MINIMUM_TEMPERATURE] = MpiWrapper::min( regData[Statistics::MINIMUM_TEMPERATURE] );
    regData[Statistics::MAXIMUM_TEMPERATURE] = MpiWrapper::max( regData[Statistics::MAXIMUM_TEMPERATURE] );

    regData[Statistics::AVERAGE_PRESSURE] = MpiWrapper::sum( regData[Statistics::AVERAGE_PRESSURE] );
    regData[Statistics::AVERAGE_TEMPERATURE] = MpiWrapper::sum( regData[Statistics::AVERAGE_TEMPERATURE] );

    for( integer propIndex = Statistics::PHASE_PORE_VOLUME; propIndex < Statistics::END; ++propIndex )
    {
      regData[propIndex] = MpiWrapper::sum( regData[propIndex] );
    }

    // Actually calculate averages
    real64 const invStaticPoreVolume = safeInverse( regData[Statistics::STATIC_PORE_VOLUME] );
    regData[Statistics::AVERAGE_PRESSURE] *= invStaticPoreVolume;
    regData[Statistics::AVERAGE_TEMPERATURE] *= invStaticPoreVolume;

    // Calculate average density and viscosity
    for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
    {
      real64 const invPhaseVolume = safeInverse( regData[Statistics::PHASE_PORE_VOLUME + phaseIndex] );
      regData[Statistics::PHASE_DENSITY + phaseIndex] *= invPhaseVolume;
      regData[Statistics::PHASE_VISCOSITY + phaseIndex] *= invPhaseVolume;
    }
  }

  if( m_writeCSV > 0 && MpiWrapper::commRank() == 0 )
  {
    array1d< real64 > propColumns;
    populateColumns( time, propColumns, regionStatistics );
    std::ofstream outputFile( m_outputDir + "/" + viewKeyStruct::fileNameString() + ".csv", std::ios_base::app );
    writeArray( outputFile, propColumns );
    outputFile.close();
  }
}

void RegionMultiphaseStatistics::initializeFile() const
{
  string const fluidModelName = m_solver->referenceFluidModelName();

  DomainPartition const & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();
  MultiFluidBase const & fluidModel = constitutiveManager.getConstitutiveRelation< MultiFluidBase >( fluidModelName );
  auto fluidComponentNames = fluidModel.componentNames();
  auto fluidPhaseNames = fluidModel.phaseNames();

  integer const phaseCount = m_solver->numFluidPhases();
  integer const componentCount = m_solver->numFluidComponents();
  integer const regionCount = m_regionIdentifiers.size();
  auto regionNames = m_regionNames.toView();

  string_array propNames;
  string_array regNames;
  string_array propUnits;
  string_array phaseNames;
  string_array compNames;

  auto addColumn = [&]( string const propName,
                        string const propUnit,
                        integer const regionIndex=-1,
                        integer const phaseIndex=-1,
                        integer const compIndex=-1 )
  {
    propNames.emplace_back( propName );
    propUnits.emplace_back( propUnit );
    if( 0 <= regionIndex )
    {
      regNames.emplace_back( regionNames[regionIndex] );
    }
    else
    {
      regNames.emplace_back( "" );
    }
    if( 0 <= phaseIndex )
    {
      phaseNames.emplace_back( fluidPhaseNames[phaseIndex] );
    }
    else
    {
      phaseNames.emplace_back( "" );
    }
    if( 0 <= compIndex )
    {
      compNames.emplace_back( fluidComponentNames[compIndex] );
    }
    else
    {
      compNames.emplace_back( "" );
    }
  };

  integer const useMass = m_solver->getReference< integer >( CompositionalMultiphaseBase::viewKeyStruct::useMassFlagString() );
  string const massUnit = useMass ? "kg" : "mol";
  string const densityUnit = useMass ? "kg/m^3" : "mol/m^3";

  addColumn( "Time", "s" );
  for( integer region = 0; region < regionCount; ++region )
  {
    for( auto const & prop : m_propertyNameTypes )
    {
      string const propName = EnumStrings< PropertyNameType >::toString( prop );
      string propUnit = "";
      if( prop == PropertyNameType::Pressure )
      {
        propUnit = "Pa";
        addColumn( GEOS_FMT( "Min {}", propName ), propUnit, region );
        addColumn( GEOS_FMT( "Average {}", propName ), propUnit, region );
        addColumn( GEOS_FMT( "Max {}", propName ), propUnit, region );
      }

      if( prop == PropertyNameType::Temperature )
      {
        propUnit = "K";
        addColumn( GEOS_FMT( "Min {}", propName ), propUnit, region );
        addColumn( GEOS_FMT( "Average {}", propName ), propUnit, region );
        addColumn( GEOS_FMT( "Max {}", propName ), propUnit, region );
      }

      if( prop == PropertyNameType::StaticPoreVolume )
      {
        propUnit = "rm^3";
        addColumn( propName, propUnit, region );
      }

      if( prop == PropertyNameType::PoreVolume )
      {
        propUnit = "rm^3";
        addColumn( propName, propUnit, region );
      }

      if( prop == PropertyNameType::PhaseVolumeFraction )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          addColumn( propName, propUnit, region, phaseIndex );
        }
      }

      if( prop == PropertyNameType::PhasePoreVolume )
      {
        propUnit = "rm^3";
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          addColumn( propName, propUnit, region, phaseIndex );
        }
      }

      if( prop == PropertyNameType::PhaseMass )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          addColumn( propName, massUnit, region, phaseIndex );
        }
      }

      if( prop == PropertyNameType::PhaseDensity )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          addColumn( propName, densityUnit, region, phaseIndex );
        }
      }

      if( prop == PropertyNameType::PhaseViscosity )
      {
        propUnit = "Pa s";
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          addColumn( propName, propUnit, region, phaseIndex );
        }
      }

      if( prop == PropertyNameType::PhaseComponentMass )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
          {
            addColumn( propName, massUnit, region, phaseIndex, compIndex );
          }
        }
      }

      if( prop == PropertyNameType::ComponentMass )
      {
        for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
        {
          addColumn( propName, massUnit, region, -1, compIndex );
        }
      }
    }
  }

  std::ofstream outputFile( m_outputDir + "/" + viewKeyStruct::fileNameString() + ".csv" );

  writeArray( outputFile, propNames );
  writeArray( outputFile, propUnits );
  writeArray( outputFile, regNames );
  writeArray( outputFile, phaseNames );
  writeArray( outputFile, compNames );

  outputFile.close();
}

void RegionMultiphaseStatistics::populateColumns( real64 const time,
                                                  array1d< real64 > & columns,
                                                  arrayView2d< real64 const > const & regionStatistics ) const
{
  integer const phaseCount = m_solver->numFluidPhases();
  integer const componentCount = m_solver->numFluidComponents();
  integer const regionCount = m_regionIdentifiers.size();

  columns.clear();

  columns.emplace_back( time );
  for( integer region = 0; region < regionCount; ++region )
  {
    for( auto const & prop : m_propertyNameTypes )
    {
      if( prop == PropertyNameType::Pressure )
      {
        columns.emplace_back( regionStatistics( region, Statistics::MINIMUM_PRESSURE ) );
        columns.emplace_back( regionStatistics( region, Statistics::AVERAGE_PRESSURE ) );
        columns.emplace_back( regionStatistics( region, Statistics::MAXIMUM_PRESSURE ) );
      }

      if( prop == PropertyNameType::Temperature )
      {
        columns.emplace_back( regionStatistics( region, Statistics::MINIMUM_TEMPERATURE ) );
        columns.emplace_back( regionStatistics( region, Statistics::AVERAGE_TEMPERATURE ) );
        columns.emplace_back( regionStatistics( region, Statistics::MAXIMUM_TEMPERATURE ) );
      }

      if( prop == PropertyNameType::StaticPoreVolume )
      {
        columns.emplace_back( regionStatistics( region, Statistics::STATIC_PORE_VOLUME ) );
      }

      if( prop == PropertyNameType::PoreVolume )
      {
        columns.emplace_back( regionStatistics( region, Statistics::PORE_VOLUME ) );
      }

      if( prop == PropertyNameType::PhaseVolumeFraction )
      {
        real64 const invPoreVolume = safeInverse( regionStatistics( region, Statistics::PORE_VOLUME ));
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          columns.emplace_back( invPoreVolume * regionStatistics( region, Statistics::PHASE_PORE_VOLUME + phaseIndex ) );
        }
      }

      if( prop == PropertyNameType::PhasePoreVolume )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          columns.emplace_back( regionStatistics( region, Statistics::PHASE_PORE_VOLUME + phaseIndex ) );
        }
      }

      if( prop == PropertyNameType::PhaseMass )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          columns.emplace_back( regionStatistics( region, Statistics::PHASE_MASS + phaseIndex ) );
        }
      }

      if( prop == PropertyNameType::PhaseDensity )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          columns.emplace_back( regionStatistics( region, Statistics::PHASE_DENSITY + phaseIndex ) );
        }
      }

      if( prop == PropertyNameType::PhaseViscosity )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          columns.emplace_back( regionStatistics( region, Statistics::PHASE_VISCOSITY + phaseIndex ) );
        }
      }

      if( prop == PropertyNameType::PhaseComponentMass )
      {
        for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
        {
          integer const phaseMassIndex = Statistics::PHASE_COMP_MASS + phaseIndex*Statistics::MAX_NUM_COMPS;
          for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
          {
            columns.emplace_back( regionStatistics( region, phaseMassIndex + compIndex ) );
          }
        }
      }

      if( prop == PropertyNameType::ComponentMass )
      {
        for( integer compIndex = 0; compIndex < componentCount; ++compIndex )
        {
          real64 componentMass = 0.0;
          for( integer phaseIndex = 0; phaseIndex < phaseCount; ++phaseIndex )
          {
            integer const phaseMassIndex = Statistics::PHASE_COMP_MASS + phaseIndex*Statistics::MAX_NUM_COMPS;
            componentMass += regionStatistics( region, phaseMassIndex + compIndex );
          }
          columns.emplace_back( componentMass );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        RegionMultiphaseStatistics,
                        string const &,
                        dataRepository::Group * const )

} /* namespace geos */
