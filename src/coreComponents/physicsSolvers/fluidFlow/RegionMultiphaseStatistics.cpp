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
//#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
//#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
//#include "constitutive/solid/CoupledSolidBase.hpp"
//#include "finiteVolume/FiniteVolumeManager.hpp"
//#include "finiteVolume/FluxApproximationBase.hpp"
//#include "mainInterface/ProblemManager.hpp"
//#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "mesh/MeshFields.hpp"
//#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVM.hpp"
//#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseBaseKernels.hpp"
//#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFVMKernels.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;

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
  GEOS_ERROR_IF_NE(m_regionNames.size(), m_regionIdentifiers.size(), 
    GEOS_FMT( "The number of region names {} is different from the number of region identifiers {}.",
    m_regionNames.size(), m_regionIdentifiers.size() ) );
  
  // Convert the property names to integers
  m_propertyNameTypes.clear();
  localIndex const propSize = m_propertyNames.size();
  if( propSize == 0 )
  {
m_propertyNameTypes.emplace_back(PropertyNameType::Pressure);
m_propertyNameTypes.emplace_back(PropertyNameType::PoreVolume);
m_propertyNameTypes.emplace_back(PropertyNameType::Saturation);
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
        subRegion.registerWrapper<array1d<real64>>( fieldName ).
             setPlotLevel( PlotLevel::NOPLOT ).
             setRestartFlags( RestartFlags::NO_WRITE ).
             setDescription( fieldName ).
             setRegisteringObjects( getName() );
      });
  });
}

bool RegionMultiphaseStatistics::execute( real64 const time_n,
                                          real64 const dt,
                                          integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                          integer const GEOS_UNUSED_PARAM( eventCounter ),
                                          real64 const GEOS_UNUSED_PARAM( eventProgress ),
                                          DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( domain );
  std::cout << "*** Executing " << getName() << "\n";
  return false;
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        RegionMultiphaseStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geos */

