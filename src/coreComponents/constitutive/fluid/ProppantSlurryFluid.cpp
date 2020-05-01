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
 * @file ProppantSlurryFluid.cpp
 */

#include "ProppantSlurryFluid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

ProppantSlurryFluid::ProppantSlurryFluid( std::string const & name, Group * const parent ):
  SlurryFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::referenceProppantDensityString, &m_referenceProppantDensity )->
    setApplyDefaultValue( 1400.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference proppant density" );

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure )->
    setApplyDefaultValue( 1e5 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceDensityString, &m_referenceDensity )->
    setApplyDefaultValue( 1000.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity )->
    setApplyDefaultValue( 0.001 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString, &m_maxProppantConcentration )->
    setApplyDefaultValue( 0.6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum proppant concentration" );

}

ProppantSlurryFluid::~ProppantSlurryFluid() = default;

void ProppantSlurryFluid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  SlurryFluidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_density = m_referenceDensity;
  m_viscosity = m_referenceViscosity;

}

void
ProppantSlurryFluid::DeliverClone( string const & name,
                                   Group * const parent,
                                   std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< ProppantSlurryFluid >( name, parent );
  }
  SlurryFluidBase::DeliverClone( name, parent, clone );
  ProppantSlurryFluid * const newConstitutiveRelation = dynamic_cast< ProppantSlurryFluid * >(clone.get());

  newConstitutiveRelation->m_compressibility      = this->m_compressibility;
  newConstitutiveRelation->m_referenceProppantDensity        = this->m_referenceProppantDensity;
  newConstitutiveRelation->m_referencePressure    = this->m_referencePressure;
  newConstitutiveRelation->m_referenceDensity     = this->m_referenceDensity;
  newConstitutiveRelation->m_referenceViscosity   = this->m_referenceViscosity;
  newConstitutiveRelation->m_maxProppantConcentration   = this->m_maxProppantConcentration;

}

void ProppantSlurryFluid::PostProcessInput()
{
  SlurryFluidBase::PostProcessInput();

  GEOSX_ERROR_IF( m_compressibility < 0.0, "An invalid value of fluid compressibility ("
                  << m_compressibility << ") is specified" );

  GEOSX_ERROR_IF( m_referenceDensity <= 0.0, "An invalid value of fluid reference density (" << m_compressibility << ") is specified" );

  GEOSX_ERROR_IF( m_referenceViscosity <= 0.0, "An invalid value of fluid reference viscosity is specified" );

  GEOSX_ERROR_IF( m_maxProppantConcentration <= 0.0 || m_maxProppantConcentration > 1.0, "An invalid value of maximum proppant volume fraction is specified" );

}

void ProppantSlurryFluid::BatchUpdate( arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM(
                                         pressure ), arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM(
                                         proppantConcentration ), arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM(
                                         componentConcentration ), arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( shearRate ))
{
  GEOSX_ERROR( "BatchUpdate for ProppantSlurryFluid is not implemented" );
}

void ProppantSlurryFluid::PointUpdate( real64 const & pressure, real64 const & proppantConcentration,
                                       arraySlice1d< real64 const > const & componentConcentration, real64 const & GEOSX_UNUSED_PARAM(
                                         shearRate ), integer const & isProppantBoundary, localIndex const k, localIndex const q )

{

  ComputeFluidDensity( pressure, componentConcentration, m_componentDensity[k][q], m_dCompDens_dPres[k][q], m_dCompDens_dCompConc[k][q], m_fluidDensity[k][q],
                       m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q] );

  ComputeFluidViscosity( m_componentDensity[k][q], m_dCompDens_dPres[k][q], m_dCompDens_dCompConc[k][q], m_fluidDensity[k][q], m_dFluidDens_dPres[k][q],
                         m_dFluidDens_dCompConc[k][q], m_fluidViscosity[k][q], m_dFluidVisc_dPres[k][q], m_dFluidVisc_dCompConc[k][q] );

  Compute( proppantConcentration, m_fluidDensity[k][q], m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q], m_fluidViscosity[k][q],
           m_dFluidVisc_dPres[k][q], m_dFluidVisc_dCompConc[k][q], isProppantBoundary, m_density[k][q], m_dDens_dPres[k][q], m_dDens_dProppantConc[k][q],
           m_dDens_dCompConc[k][q], m_viscosity[k][q], m_dVisc_dPres[k][q], m_dVisc_dProppantConc[k][q], m_dVisc_dCompConc[k][q] );

}

void ProppantSlurryFluid::PointUpdateFluidProperty( real64 const & pressure, arraySlice1d< real64 const > const & componentConcentration, real64 const & GEOSX_UNUSED_PARAM(
                                                      shearRate ), localIndex const k, localIndex const q )
{

  ComputeFluidDensity( pressure, componentConcentration, m_componentDensity[k][q], m_dCompDens_dPres[k][q], m_dCompDens_dCompConc[k][q], m_fluidDensity[k][q],
                       m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q] );

  ComputeFluidViscosity( m_componentDensity[k][q], m_dCompDens_dPres[k][q], m_dCompDens_dCompConc[k][q], m_fluidDensity[k][q], m_dFluidDens_dPres[k][q],
                         m_dFluidDens_dCompConc[k][q], m_fluidViscosity[k][q], m_dFluidVisc_dPres[k][q], m_dFluidVisc_dCompConc[k][q] );

}

void ProppantSlurryFluid::PointUpdateComponentDensity( real64 const & pressure, arraySlice1d< real64 const > const & componentConcentration, localIndex const k,
                                                       localIndex const q )
{

  ComputeComponentDensity( pressure, componentConcentration, m_componentDensity[k][q], m_dCompDens_dPres[k][q], m_dCompDens_dCompConc[k][q] );

}



void ProppantSlurryFluid::ComputeFluidDensity( real64 const & pressure,
                                               arraySlice1d< real64 const > const & componentConcentration,
                                               arraySlice1d< real64 > const & componentDensity,
                                               arraySlice1d< real64 > const & dComponentDensity_dPressure,
                                               arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
                                               real64 & fluidDensity,
                                               real64 & dFluidDensity_dPressure,
                                               arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const
{

  localIndex const NC = numFluidComponents();

  real64 baseFluidDensity = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure));

  real64 dBaseFluidDensity_dPressure = m_compressibility * baseFluidDensity;

  fluidDensity = baseFluidDensity;

  dFluidDensity_dPressure = dBaseFluidDensity_dPressure;

  for( localIndex c = 0; c < NC; ++c )
    dFluidDensity_dComponentConcentration[c] = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    componentDensity[c] = componentConcentration[c] * m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

    dComponentDensity_dPressure[c] = m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
      dComponentDensity_dComponentConcentration[c][i] = 0.0;

    dComponentDensity_dComponentConcentration[c][c] = m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

    fluidDensity += componentDensity[c] - componentConcentration[c] * baseFluidDensity;

    dFluidDensity_dPressure += dComponentDensity_dPressure[c] - componentConcentration[c] * dBaseFluidDensity_dPressure;

    dFluidDensity_dComponentConcentration[c] += dComponentDensity_dComponentConcentration[c][c] - baseFluidDensity;

  }


  for( localIndex c = 0; c < NC; ++c )
    dFluidDensity_dComponentConcentration[c] = 0.0;

}

void ProppantSlurryFluid::ComputeComponentDensity( real64 const & pressure,
                                                   arraySlice1d< real64 const > const & componentConcentration,
                                                   arraySlice1d< real64 > const & componentDensity,
                                                   arraySlice1d< real64 > const & dComponentDensity_dPressure,
                                                   arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const
{

  localIndex const NC = numFluidComponents();

  for( localIndex c = 0; c < NC; ++c )
  {

    componentDensity[c] = componentConcentration[c] * m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

    dComponentDensity_dPressure[c] = m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
      dComponentDensity_dComponentConcentration[c][i] = 0.0;

    dComponentDensity_dComponentConcentration[c][c] = m_defaultDensity[c] * exp( m_defaultCompressibility[c] * (pressure - m_referencePressure));

  }

}


void ProppantSlurryFluid::ComputeFluidViscosity( arraySlice1d< real64 const > const & componentDensity,
                                                 arraySlice1d< real64 const > const & dComponentDensity_dPressure,
                                                 arraySlice2d< real64 const > const & GEOSX_UNUSED_PARAM( dComponentDensity_dComponentConcentration ),
                                                 real64 const & fluidDensity,
                                                 real64 const & dFluidDensity_dPressure,
                                                 arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                                                 real64 & fluidViscosity,
                                                 real64 & dFluidViscosity_dPressure,
                                                 arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const
{

  localIndex const NC = numFluidComponents();

  fluidViscosity = m_referenceViscosity;
  dFluidViscosity_dPressure = 0.0;

  for( localIndex c = 0; c < NC; ++c )
    dFluidViscosity_dComponentConcentration[c] = 0.0;

  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {

    fluidViscosity += componentDensity[c1] / fluidDensity * (m_defaultViscosity[c1] - m_referenceViscosity);

    dFluidViscosity_dPressure +=
      (dComponentDensity_dPressure[c1] / fluidDensity - componentDensity[c1] / fluidDensity / fluidDensity * dFluidDensity_dPressure) *
      (m_defaultViscosity[c1] - m_referenceViscosity);
    /*
       for(localIndex c2 = 0; c2 < NC; ++c2)
       {

        dFluidViscosity_dComponentConcentration[c2] += (dComponentDensity_dComponentConcentration[c1][c2] / fluidDensity
           - componentDensity[c1] / fluidDensity / fluidDensity * dFluidDensity_dComponentConcentration[c2]) *
           (m_defaultViscosity[c1] - m_referenceViscosity);

       }
     */
  }

}


void ProppantSlurryFluid::Compute( real64 const & proppantConcentration,
                                   real64 const & fluidDensity,
                                   real64 const & dFluidDensity_dPressure,
                                   arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                                   real64 const & fluidViscosity,
                                   real64 const & dFluidViscosity_dPressure,
                                   arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidViscosity_dComponentConcentration ),
                                   integer const & isProppantBoundary,
                                   real64 & density,
                                   real64 & dDensity_dPressure,
                                   real64 & dDensity_dProppantConcentration,
                                   arraySlice1d< real64 > const & dDensity_dComponentConcentration,
                                   real64 & viscosity,
                                   real64 & dViscosity_dPressure,
                                   real64 & dViscosity_dProppantConcentration,
                                   arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const

{

  localIndex const NC = numFluidComponents();


  real64 effectiveConcentration = proppantConcentration;
  if( effectiveConcentration > 0.95 * m_maxProppantConcentration || isProppantBoundary == 1 )
    effectiveConcentration = 0.0;

  density = (1.0 - effectiveConcentration) * fluidDensity + effectiveConcentration * m_referenceProppantDensity;

  dDensity_dPressure = (1.0 - effectiveConcentration) * dFluidDensity_dPressure;

  /*
     dDensity_dProppantConcentration = -fluidDensity + m_referenceProppantDensity;

     for(localIndex c = 0; c < NC; ++c)
     {

      dDensity_dComponentConcentration[c] = (1.0 - effectiveConcentration) * dFluidDensity_dComponentConcentration[c];

     }
   */

  dDensity_dProppantConcentration = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    dDensity_dComponentConcentration[c] = 0.0;

  }


  real64 coef = pow( 1.0 + 1.25 *  effectiveConcentration / (1.0 - effectiveConcentration / m_maxProppantConcentration), 2.0 );

  viscosity = fluidViscosity * coef;

  dViscosity_dPressure = dFluidViscosity_dPressure * coef;

  /*
     dViscosity_dProppantConcentration = fluidViscosity * 2.0 * (1.0 + 1.25 *  effectiveConcentration / (1.0 -
        effectiveConcentration / m_maxProppantConcentration)) * 1.25 * m_maxProppantConcentration *
        m_maxProppantConcentration / (m_maxProppantConcentration - effectiveConcentration) / (m_maxProppantConcentration
        - effectiveConcentration);

     for(localIndex c = 0; c < NC; ++c)
     dViscosity_dComponentConcentration[c] = dFluidViscosity_dComponentConcentration[c] * coef;
   */

  dViscosity_dProppantConcentration = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    dViscosity_dComponentConcentration[c] = 0.0;

  }

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantSlurryFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
