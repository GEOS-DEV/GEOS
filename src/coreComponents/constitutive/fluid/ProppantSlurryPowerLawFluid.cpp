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
{}

void ProppantSlurryFluid::PointUpdate( real64 const & pressure, real64 const & proppantConcentration,
                                       arraySlice1d< real64 const > const & componentConcentration, real64 const & shearRate, localIndex const k,
                                       localIndex const q )

{

  localIndex const NC = numFluidComponents();

  ComputeFluidDensity( NC, pressure, componentConcentration, m_fluidDensity[k][q], m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q] );

  Compute( NC, pressure, proppantConcentration, componentConcentration, shearRate, m_fluidDensity[k][q], m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q],
           m_density[k][q], m_dDens_dPres[k][q], m_dDens_dProppantConc[k][q], m_dDens_dCompConc[k][q], m_viscosity[k][q], m_dVisc_dPres[k][q],
           m_dVisc_dProppantConc[k][q], m_dVisc_dCompConc[k][q] );

}

void ProppantSlurryFluid::PointUpdateFluidDensity( real64 const & pressure, arraySlice1d< real64 const > const & componentConcentration, localIndex const k,
                                                   localIndex const q )
{

  localIndex const NC = numFluidComponents();

  ComputeFluidDensity( NC, pressure, componentConcentration, m_fluidDensity[k][q], m_dFluidDens_dPres[k][q], m_dFluidDens_dCompConc[k][q] );

}


void ProppantSlurryFluid::ComputeFluidDensity( localIndex const NC,
                                               real64 const & pressure,
                                               arraySlice1d< real64 const > const & componentConcentration,
                                               real64 & fluidDensity,
                                               real64 & dFluidDensity_dPressure,
                                               arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const
{

  real64 fluidConcentration = 1.0;

  // density

  real64 waterDensity = m_referenceDensity * exp( m_compressibility * (pressure - m_referencePressure));

  real64 dWaterDensity_dPres = m_compressibility * waterDensity;

  fluidDensity = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    fluidDensity += componentConcentration[c] * m_defaultDensity[c];
    dFluidDensity_dComponentConcentration[c] = m_defaultDensity[c];

    fluidConcentration -= componentConcentration[c];

  }

  fluidDensity += fluidConcentration * waterDensity;

  dFluidDensity_dPressure = fluidConcentration * dWaterDensity_dPres;


}


void ProppantSlurryFluid::Compute( localIndex const NC,
                                   real64 const & GEOSX_UNUSED_PARAM( pressure ),
                                   real64 const & proppantConcentration,
                                   arraySlice1d< real64 const > const & componentConcentration,
                                   real64 const & shearRate,
                                   real64 const & fluidDensity,
                                   real64 const & dFluidDensity_dPressure,
                                   arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                                   real64 & density,
                                   real64 & dDensity_dPressure,
                                   real64 & dDensity_dProppantConcentration,
                                   arraySlice1d< real64 > const & dDensity_dComponentConcentration,
                                   real64 & viscosity,
                                   real64 & dViscosity_dPressure,
                                   real64 & dViscosity_dProppantConcentration,
                                   arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const

{

  static real64 eps = 1e-6;

  real64 fluidConcentration = 1.0;


  for( localIndex c = 0; c < NC; ++c )
  {

    fluidConcentration -= componentConcentration[c];

  }

  density = (1.0 - proppantConcentration) * fluidDensity + proppantConcentration * m_referenceProppantDensity;

  dDensity_dPressure = (1.0 - proppantConcentration) * dFluidDensity_dPressure;

  dDensity_dProppantConcentration = -fluidDensity + m_referenceProppantDensity;

  for( localIndex c = 0; c < NC; ++c )
  {

    dDensity_dComponentConcentration[c] = (1.0 - proppantConcentration) * dFluidDensity_dComponentConcentration[c];

  }

  real64 nIndex = NC > 0 ? 0.0 : 1.0;
  real64 K = NC > 0 ? 0.0 : m_referenceViscosity;

  array1d< real64 > dNIndex_dC( NC );
  array1d< real64 > dK_dC( NC );

  if( fluidConcentration < 1.0 )
  {

    for( localIndex c = 0; c < NC; ++c )
    {

      nIndex +=  componentConcentration[c] * m_nIndices[c] / (1.0 - fluidConcentration);
      K +=  componentConcentration[c] * m_Ks[c] / (1.0 - fluidConcentration);
      dNIndex_dC[c] = m_nIndices[c] / (1.0 - fluidConcentration);
      dK_dC[c] = m_Ks[c] / (1.0 - fluidConcentration);

    }

    for( localIndex c = 0; c < NC; ++c )
    {

      dNIndex_dC[c] -= nIndex / (1.0 - fluidConcentration);
      dK_dC[c] -= K / (1.0 - fluidConcentration);

    }
  }

  real64 fluidViscosity = 0.0;

  array1d< real64 > dFluidViscosity_dC( NC );

  bool isNewtonian = (fabs( nIndex - 1.0 ) < eps || shearRate < eps || K < eps || nIndex < eps) ? 1 : 0;


  if( isNewtonian )
  {

    fluidViscosity = m_referenceViscosity;

    for( localIndex c = 0; c < NC; ++c )
    {

      dFluidViscosity_dC[c] = 0.0;

    }
  }
  else if( nIndex > 1.0 )
  {


    fluidViscosity = m_referenceViscosity * (1.0 + K / m_referenceViscosity * pow( shearRate, nIndex - 1.0 ));

    for( localIndex c = 0; c < NC; ++c )
    {

      dFluidViscosity_dC[c] = dK_dC[c] * pow( shearRate, nIndex - 1.0 ) + K * log( shearRate ) * pow( shearRate, nIndex - 1.0 ) * dNIndex_dC[c];

    }

  }
  else
  {

    fluidViscosity = m_referenceViscosity / (1.0 + m_referenceViscosity / K * pow( shearRate, 1.0 - nIndex ));

    for( localIndex c = 0; c < NC; ++c )
    {

      dFluidViscosity_dC[c] = fluidViscosity * m_referenceViscosity *
                              pow( shearRate, 1.0 - nIndex ) * (log( shearRate ) * dNIndex_dC[c] / K + dK_dC[c] / K / K);
    }

  }


  real64 coef = 0.0;

  if( isNewtonian )
  {

    viscosity = fluidViscosity * pow( 1.0 + 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration), 2.0 );


    dViscosity_dProppantConcentration = fluidViscosity * 2.0 *
                                        (1.0 + 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration)) * 1.25 *
                                        m_maxProppantConcentration *
                                        m_maxProppantConcentration / (m_maxProppantConcentration - proppantConcentration) /
                                        (m_maxProppantConcentration - proppantConcentration);

    for( localIndex c = 0; c < NC; ++c )
      dViscosity_dComponentConcentration[c] = 0.0;

  }
  else
  {

    coef = 0.75 * (exp( 1.5 * nIndex ) - 1.0) * exp( -(1.0 - nIndex) * shearRate / 1000.0 );

    viscosity = fluidViscosity * pow( 1.0 +  coef * 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration), 2.0 );

    dViscosity_dProppantConcentration = fluidViscosity * 2.0 *
                                        (1.0 + coef * 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration)) * 1.25 *
                                        m_maxProppantConcentration *
                                        m_maxProppantConcentration / (m_maxProppantConcentration - proppantConcentration) /
                                        (m_maxProppantConcentration - proppantConcentration);

    for( localIndex c = 0; c < NC; ++c )
    {

      dViscosity_dComponentConcentration[c] = dFluidViscosity_dC[c] *
                                              pow( 1.0 +  coef * 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration),
                                                   2.0 ) + fluidViscosity * 2.0 *
                                              (1.0 + coef * 1.25 *  proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration)) * 1.25 *
                                              proppantConcentration / (1.0 - proppantConcentration / m_maxProppantConcentration) *
                                              (0.75 * (exp( 1.5 * nIndex ) * 1.5 * dNIndex_dC[c]) * exp( -(1.0 - nIndex) * shearRate / 1000.0 ) +
                                               (exp( 1.5 * nIndex ) - 1.0) *
                                               exp( -(1.0 - nIndex) * shearRate / 1000.0 ) * shearRate / 1000.0 * dNIndex_dC[c]);

    }

  }

  dViscosity_dPressure = 0.0;

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantSlurryFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
