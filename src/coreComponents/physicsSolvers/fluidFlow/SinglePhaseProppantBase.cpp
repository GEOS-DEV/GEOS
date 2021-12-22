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
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/ProppantSolid.hpp"
#include "constitutive/solid/porosity/ProppantPorosity.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/proppantTransport/ProppantTransport.hpp"
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
           slurryFluid.getWrapper< array2d< real64 > >( SlurryFluidBase::viewKeyStruct::densityString() ).getDefaultValue(),
           slurryFluid.getWrapper< array2d< real64 > >( SlurryFluidBase::viewKeyStruct::viscosityString() ).getDefaultValue() };
}

void SinglePhaseProppantBase::updateFluidModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();

  arrayView1d< real64 const > const dPres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

  arrayView1d< real64 const > const proppantConcentration =
    dataGroup.getReference< array1d< real64 > >( ProppantTransport::viewKeyStruct::proppantConcentrationString() );

  arrayView1d< real64 const > const dProppantConcentration =
    dataGroup.getReference< array1d< real64 > >( ProppantTransport::viewKeyStruct::deltaProppantConcentrationString() );

  arrayView2d< real64 const > const componentConcentration =
    dataGroup.getReference< array2d< real64 > >( ProppantTransport::viewKeyStruct::componentConcentrationString() );

  arrayView2d< real64 const > const cellBasedFlux =
    dataGroup.getReference< array2d< real64 > >( ProppantTransport::viewKeyStruct::cellBasedFluxString() );

  arrayView1d< integer const > const isProppantBoundaryElement =
    dataGroup.getReference< array1d< integer > >( ProppantTransport::viewKeyStruct::isProppantBoundaryString() );

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
    subRegion.getReference< array1d< real64 > >( ProppantTransport::viewKeyStruct::proppantPackVolumeFractionString() );

  arrayView1d< real64 const > const newHydraulicAperture =
    subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::hydraulicAperture::key() );
  arrayView1d< real64 const > const oldHydraulicAperture =
    subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::aperture0::key() );

  CoupledSolidBase & porousSolid = subRegion.template getConstitutiveModel< CoupledSolidBase >( m_solidModelNames[targetIndex] );

  constitutive::ConstitutivePassThru< ProppantSolid< ProppantPorosity, ProppantPermeability > >::execute( porousSolid, [=, &subRegion] ( auto & castedProppantSolid )
  {
    typename TYPEOFREF( castedProppantSolid ) ::KernelWrapper porousWrapper = castedProppantSolid.createKernelUpdates();

    execute3( porousWrapper, subRegion, newHydraulicAperture, oldHydraulicAperture, proppantPackVolumeFraction );

  } );

}

void SinglePhaseProppantBase::resetViewsPrivate( ElementRegionManager const & elemManager )
{
  {
    using keys = SlurryFluidBase::viewKeyStruct;

    m_density.clear();
    m_density = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::densityString(),
                                                                             targetRegionNames(),
                                                                             fluidModelNames() );
    m_density.setName( getName() + "/accessors/" + keys::densityString() );

    m_dDens_dPres.clear();
    m_dDens_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dDens_dPresString(),
                                                                                 targetRegionNames(),
                                                                                 fluidModelNames() );
    m_dDens_dPres.setName( getName() + "/accessors/" + keys::dDens_dPresString() );

    m_viscosity.clear();
    m_viscosity = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::viscosityString(),
                                                                               targetRegionNames(),
                                                                               fluidModelNames() );
    m_viscosity.setName( getName() + "/accessors/" + keys::viscosityString() );

    m_dVisc_dPres.clear();
    m_dVisc_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dVisc_dPresString(),
                                                                                 targetRegionNames(),
                                                                                 fluidModelNames() );
    m_dVisc_dPres.setName( getName() + "/accessors/" + keys::dVisc_dPresString() );
  }

  {
    m_permeabilityMultiplier.clear();
    m_permeabilityMultiplier = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( extrinsicMeshData::permeability::permeabilityMultiplier::key(),
                                                                                            targetRegionNames(),
                                                                                            m_permeabilityModelNames );
    m_permeabilityMultiplier.setName( getName() + "/accessors/" + extrinsicMeshData::permeability::permeabilityMultiplier::key() );
  }
}

}
