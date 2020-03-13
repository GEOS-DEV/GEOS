/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseProppantBase.cpp
 */


#include "SinglePhaseProppantBase.hpp"

#include "ProppantTransport.hpp"

#include "constitutive/fluid/SlurryFluidBase.hpp"

namespace geosx
{

using constitutive::SlurryFluidBase;

SinglePhaseProppantBase::SinglePhaseProppantBase( const std::string & name,
                                                  Group * const parent ):
  SinglePhaseBase( name, parent )
{}

SinglePhaseProppantBase::~SinglePhaseProppantBase()
{
  // TODO Auto-generated destructor stub
}

void SinglePhaseProppantBase::UpdateState( Group * dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateMobility< SlurryFluidBase >( dataGroup );
}

void SinglePhaseProppantBase::UpdateFluidModel( Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  SlurryFluidBase * const fluid = GetConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidName );

  arrayView1d< real64 const > const & proppantConcentration = dataGroup->getReference< array1d< real64 > >(
    ProppantTransport::viewKeyStruct::proppantConcentrationString );

  arrayView1d< real64 const > const & dProppantConcentration = dataGroup->getReference< array1d< real64 > >(
    ProppantTransport::viewKeyStruct::deltaProppantConcentrationString );

  arrayView2d< real64 const > const & componentConcentration = dataGroup->getReference< array2d< real64 > >(
    ProppantTransport::viewKeyStruct::componentConcentrationString );

  arrayView1d< R1Tensor const > const & cellBasedFlux = dataGroup->getReference< array1d< R1Tensor > >( ProppantTransport::viewKeyStruct::cellBasedFluxString );

  arrayView1d< integer const > const & isProppantBoundaryElement  = dataGroup->getReference< array1d< integer > >(
    ProppantTransport::viewKeyStruct::isProppantBoundaryString );

  forall_in_range< RAJA::seq_exec >( 0, dataGroup->size(), [=] ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], proppantConcentration[a] + dProppantConcentration[a], componentConcentration[a], cellBasedFlux[a].L2_Norm(),
                        isProppantBoundaryElement[a], a, 0 );
  } );
}


void SinglePhaseProppantBase::ResetViewsPrivate( ElementRegionManager * const elemManager,
                                                 constitutive::ConstitutiveManager * const constitutiveManager )
{
  m_density =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SlurryFluidBase::viewKeyStruct::densityString,
                                                                                                constitutiveManager );
  m_dDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SlurryFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                                constitutiveManager );
  m_viscosity =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SlurryFluidBase::viewKeyStruct::viscosityString,
                                                                                                constitutiveManager );
  m_dVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SlurryFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                                constitutiveManager );

  m_poroMultiplier =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( ProppantTransport::viewKeyStruct::poroMultiplierString );

  m_transTMultiplier =
    elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( ProppantTransport::viewKeyStruct::transTMultiplierString );
}

}
