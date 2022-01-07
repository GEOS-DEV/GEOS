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
 * @file SinglePhaseProppantBase.cpp
 */


#include "SinglePhaseProppantBase.hpp"

#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/slurryFluidSelector.hpp"
#include "constitutive/fluid/SlurryFluidExtrinsicData.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/ProppantSolid.hpp"
#include "constitutive/solid/porosity/ProppantPorosity.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransportExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBaseKernels.hpp"

namespace geosx
{

using namespace constitutive;

using constitutive::SlurryFluidBase;

template< typename POROUSWRAPPER_TYPE >
void execute3( POROUSWRAPPER_TYPE porousWrapper,
               SurfaceElementSubRegion & subRegion,
               arrayView1d< real64 const > const & oldHydraulicAperture,
               arrayView1d< real64 const > const & newHydraulicAperture,
               arrayView1d< real64 const > const & proppantPackVolumeFraction )
{
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < porousWrapper.numGauss(); ++q )
    {
      porousWrapper.updateStateFromApertureAndProppantVolumeFraction( k, q,
                                                                      oldHydraulicAperture[k],
                                                                      newHydraulicAperture[k],
                                                                      proppantPackVolumeFraction[k] );
    }
  } );
}

SinglePhaseProppantBase::SinglePhaseProppantBase( const string & name,
                                                  Group * const parent ):
  SinglePhaseBase( name, parent )
{}

SinglePhaseProppantBase::~SinglePhaseProppantBase()
{}

void SinglePhaseProppantBase::validateFluidModels( DomainPartition const & domain ) const
{
  // Validate fluid models in regions
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel const & meshLevel = dynamicCast< MeshBody const * >( mesh.second )->getMeshLevel( 0 );
    validateModelMapping< SlurryFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
  }
}

SinglePhaseBase::FluidPropViews SinglePhaseProppantBase::getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const
{
  SlurryFluidBase const & slurryFluid = dynamicCast< SlurryFluidBase const & >( fluid );
  return { slurryFluid.density(),
           slurryFluid.dDensity_dPressure(),
           slurryFluid.viscosity(),
           slurryFluid.dViscosity_dPressure(),
           slurryFluid.getExtrinsicData< extrinsicMeshData::slurryfluid::density >().getDefaultValue(),
           slurryFluid.getExtrinsicData< extrinsicMeshData::slurryfluid::viscosity >().getDefaultValue() };
}

void SinglePhaseProppantBase::updateFluidModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();

  arrayView1d< real64 const > const dPres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

  arrayView1d< real64 const > const proppantConcentration =
    dataGroup.getExtrinsicData< extrinsicMeshData::proppant::proppantConcentration >();

  arrayView1d< real64 const > const dProppantConcentration =
    dataGroup.getExtrinsicData< extrinsicMeshData::proppant::deltaProppantConcentration >();

  arrayView2d< real64 const > const componentConcentration =
    dataGroup.getExtrinsicData< extrinsicMeshData::proppant::componentConcentration >();

  arrayView2d< real64 const > const cellBasedFlux =
    dataGroup.getExtrinsicData< extrinsicMeshData::proppant::cellBasedFlux >();

  arrayView1d< integer const > const isProppantBoundaryElement =
    dataGroup.getExtrinsicData< extrinsicMeshData::proppant::isProppantBoundary >();

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    SinglePhaseProppantBaseKernels::FluidUpdateKernel::launch( fluidWrapper,
                                                               pres,
                                                               dPres,
                                                               proppantConcentration,
                                                               dProppantConcentration,
                                                               componentConcentration,
                                                               cellBasedFlux,
                                                               isProppantBoundaryElement );
  } );
}


void SinglePhaseProppantBase::updatePorosityAndPermeability( SurfaceElementSubRegion & subRegion,
                                                             localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const proppantPackVolumeFraction =
    subRegion.getExtrinsicData< extrinsicMeshData::proppant::proppantPackVolumeFraction >();

  arrayView1d< real64 const > const newHydraulicAperture =
    subRegion.getExtrinsicData< extrinsicMeshData::flow::hydraulicAperture >();
  arrayView1d< real64 const > const oldHydraulicAperture =
    subRegion.getExtrinsicData< extrinsicMeshData::flow::aperture0 >();

  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( m_solidModelNames[targetIndex] );

  constitutive::ConstitutivePassThru< ProppantSolid< ProppantPorosity, ProppantPermeability > >::execute( porousSolid, [=, &subRegion] ( auto & castedProppantSolid )
  {
    typename TYPEOFREF( castedProppantSolid ) ::KernelWrapper porousWrapper = castedProppantSolid.createKernelUpdates();

    execute3( porousWrapper, subRegion, newHydraulicAperture, oldHydraulicAperture, proppantPackVolumeFraction );

  } );

}


}
