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
 * @file ParticleFluid.cpp
 */

#include "ParticleFluid.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;

namespace constitutive
{

ParticleFluid::ParticleSettlingModel ParticleFluid::stringToParticleSettlingModel( string const & str )
{
  if( str == "Stokes" )
  {
    return ParticleFluid::ParticleSettlingModel::Stokes;
  }
  else if( str == "Intermediate" )
  {
    return ParticleFluid::ParticleSettlingModel::Intermediate;
  }
  else
  {
    GEOSX_ERROR( "Unrecognized particle settling velocity model: " << str );
  }
  return ParticleFluid::ParticleSettlingModel::Stokes;
}

ParticleFluid::ParticleFluid( std::string const & name, Group * const parent ):
  ParticleFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::particleSettlingModelString, &m_particleSettlingModelString )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Particle settling velocity model" );

  registerWrapper( viewKeyStruct::proppantDensityString, &m_proppantDensity )->
    setApplyDefaultValue( 1400.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::fluidViscosityString, &m_fluidViscosity )->
    setApplyDefaultValue( 0.001 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fluid viscosity" );

  registerWrapper( viewKeyStruct::proppantDiameterString, &m_proppantDiameter )->
    setApplyDefaultValue( 200e-6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::hinderedSettlingCoefficientString, &m_hinderedSettlingCoefficient )->
    setApplyDefaultValue( 5.9 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Hindered settling coefficient" );

  registerWrapper( viewKeyStruct::collisionAlphaString, &m_collisionAlpha )->
    setApplyDefaultValue( 1.27 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Collision alpha coefficient" );

  registerWrapper( viewKeyStruct::slipConcentrationString, &m_slipConcentration )->
    setApplyDefaultValue( 0.1 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Slip concentration" );

  registerWrapper( viewKeyStruct::collisionBetaString, &m_collisionBeta )->
    setApplyDefaultValue( 1.5 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Collision beta coefficient" );

  registerWrapper( viewKeyStruct::bridgingFactorString, &m_bridgingFactor )->
    setApplyDefaultValue( 3.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Bridging factor" );

  registerWrapper( viewKeyStruct::sphericityString, &m_sphericity )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Sphericity" );

}

ParticleFluid::~ParticleFluid() = default;

void ParticleFluid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                              localIndex const numConstitutivePointsPerParentIndex )
{
  ParticleFluidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

}

void
ParticleFluid::DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< ParticleFluid >( name, parent );
  }
  ParticleFluidBase::DeliverClone( name, parent, clone );
  ParticleFluid * const newConstitutiveRelation = dynamic_cast< ParticleFluid * >(clone.get());

  newConstitutiveRelation->m_proppantDensity  = this->m_proppantDensity;
  newConstitutiveRelation->m_fluidViscosity   = this->m_fluidViscosity;
  newConstitutiveRelation->m_proppantDiameter = this->m_proppantDiameter;

  newConstitutiveRelation->m_hinderedSettlingCoefficient = this->m_hinderedSettlingCoefficient;
  newConstitutiveRelation->m_collisionAlpha = this->m_collisionAlpha;
  newConstitutiveRelation->m_slipConcentration = this->m_slipConcentration;
  newConstitutiveRelation->m_collisionBeta = this->m_collisionBeta;

  newConstitutiveRelation->m_bridgingFactor = this->m_bridgingFactor;
  newConstitutiveRelation->m_sphericity = this->m_sphericity;

  newConstitutiveRelation->m_particleSettlingModelString = this->m_particleSettlingModelString;
  newConstitutiveRelation->m_particleSettlingModel = this->m_particleSettlingModel;

}

void ParticleFluid::PostProcessInput()
{
  ParticleFluidBase::PostProcessInput();

  m_packPermeabilityCoef = pow( m_sphericity * m_proppantDiameter, 2.0 ) / 180.0;

  m_bridgingAperture = m_bridgingFactor * m_proppantDiameter;

  m_particleSettlingModel = stringToParticleSettlingModel( m_particleSettlingModelString );

}

void ParticleFluid::BatchUpdate( arrayView1d< real64 const > const & )
{}


void ParticleFluid::PointUpdate( localIndex const NC, real64 const & proppantConcentration, arraySlice1d< real64 const > const & componentConcentration,
                                 arraySlice1d< real64 const > const & nIndex, arraySlice1d< real64 const > const & KIndex, real64 const & fluidDensity,
                                 real64 const & dFluidDensity_dPressure, arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                                 localIndex const k )
{

  Compute( NC, proppantConcentration, componentConcentration, nIndex, KIndex, fluidDensity, dFluidDensity_dPressure, dFluidDensity_dComponentConcentration,
           m_settlingFactor[k], m_dSettlingFactor_dPressure[k], m_dSettlingFactor_dProppantConcentration[k], m_dSettlingFactor_dComponentConcentration[k],
           m_collisionFactor[k], m_dCollisionFactor_dProppantConcentration[k] );

}

void ParticleFluid::BatchUpdateMob( arrayView1d< real64 const > const &, arrayView1d< real64 const > const & )
{}

void ParticleFluid::PointUpdateMob( real64 const & concentration, real64 const & aperture, localIndex const k )
{
  ComputeMob( concentration, aperture, m_isProppantMobile[k], m_proppantPackPermeability[k] );

}

void ParticleFluid::Compute( localIndex const NC,
                             real64 const & proppantConcentration,
                             arraySlice1d< real64 const > const & componentConcentration,
                             arraySlice1d< real64 const > const & nIndices,
                             arraySlice1d< real64 const > const & Ks,
                             real64 const & fluidDensity,
                             real64 const & dFluidDensity_dPressure,
                             arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
                             real64 & settlingFactor,
                             real64 & dSettlingFactor_dPressure,
                             real64 & dSettlingFactor_dProppantConcentration,
                             arraySlice1d< real64 > const & dSettlingFactor_dComponentConcentration,
                             real64 & collisionFactor,
                             real64 & dCollisionFactor_dProppantConcentration ) const
{

  static real64 eps = 1e-6;

  real64 singleParticleSettlingVelocity = 0.0;
  real64 dSingleParticleSettlingVelocity_dPressure = 0.0;
  array1d< real64 > dSingleParticleSettlingVelocity_dComponentConcentration( NC );

  real64 fluidConcentration = 1.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    fluidConcentration -= componentConcentration[c];

  }

  real64 nIndex = NC > 0 ? 0.0 : 1.0;
  real64 K = 0.0;

  array1d< real64 > dNIndex_dC( NC );
  array1d< real64 > dK_dC( NC );

  if( fluidConcentration < 1.0 )
  {
    for( localIndex c = 0; c < NC; ++c )
    {

      nIndex +=  componentConcentration[c] * nIndices[c] / (1.0 - fluidConcentration);
      K +=  componentConcentration[c] * Ks[c] / (1.0 - fluidConcentration);
      dNIndex_dC[c] = nIndices[c] / (1.0 - fluidConcentration);
      dK_dC[c] = Ks[c] / (1.0 - fluidConcentration);

    }

    for( localIndex c = 0; c < NC; ++c )
    {

      dNIndex_dC[c] -= nIndex / (1.0 - fluidConcentration);
      dK_dC[c] -= K / (1.0 - fluidConcentration);

    }
  }

  bool isNewtonian = (fabs( nIndex - 1.0 ) < eps || nIndex < eps) ? 1 : 0;

  switch( m_particleSettlingModel )
  {
    case ParticleSettlingModel::Stokes:
    {
      if( isNewtonian )
      {

        singleParticleSettlingVelocity = 9.81 * (m_proppantDensity - fluidDensity ) * m_proppantDiameter * m_proppantDiameter / (18.0 * m_fluidViscosity);

        dSingleParticleSettlingVelocity_dPressure = -9.81 * dFluidDensity_dPressure * m_proppantDiameter * m_proppantDiameter / (18.0 * m_fluidViscosity);


        for( localIndex c = 0; c < NC; ++c )
        {

          dSingleParticleSettlingVelocity_dComponentConcentration[c] = -9.81 * dFluidDensity_dComponentConcentration[c] * m_proppantDiameter *
                                                                       m_proppantDiameter / (18.0 * m_fluidViscosity);

        }

      }
      else
      {

        singleParticleSettlingVelocity = pow( 9.81 * (m_proppantDensity - fluidDensity ) * m_proppantDiameter /18.0 * K, 1.0 / nIndex ) * m_proppantDiameter;

        dSingleParticleSettlingVelocity_dPressure = -singleParticleSettlingVelocity * (1.0 / nIndex) / (m_proppantDensity - fluidDensity ) *
                                                    dFluidDensity_dPressure;

        for( localIndex c = 0; c < NC; ++c )
        {

          dSingleParticleSettlingVelocity_dComponentConcentration[c] = singleParticleSettlingVelocity *
                                                                       (1.0 / nIndex *
                                                                        (-dFluidDensity_dComponentConcentration[c] / (m_proppantDensity - fluidDensity) +
                                                                         dK_dC[c] / K) -
                                                                        log( 9.81 * (m_proppantDensity - fluidDensity ) * m_proppantDiameter /18.0 * K ) /
                                                                        nIndex / nIndex * dNIndex_dC[c]);

        }

      }
    }
    break;
    case ParticleSettlingModel::Intermediate:
    {

      singleParticleSettlingVelocity = 20.34 *
                                       pow( m_proppantDensity - fluidDensity, 0.71 ) * pow( m_proppantDiameter, 1.14 ) / pow( fluidDensity, 0.29 ) / pow(
        m_fluidViscosity, 0.43 );

    }
    break;
    default:
      GEOSX_ERROR( "Particle settling model type not supported" );
  }


  settlingFactor = 0.0;

  dSettlingFactor_dProppantConcentration = 0.0;

  collisionFactor = 0.0;

  dCollisionFactor_dProppantConcentration = 0.0;

  if( proppantConcentration >= 0.0 && proppantConcentration < m_maxProppantConcentration )
  {

    // settlingFactor
    // Stokes + hindered settling model

    settlingFactor = singleParticleSettlingVelocity * exp( -m_hinderedSettlingCoefficient * proppantConcentration );

    dSettlingFactor_dPressure = dSingleParticleSettlingVelocity_dPressure * exp( -m_hinderedSettlingCoefficient * proppantConcentration );

    dSettlingFactor_dProppantConcentration = -m_hinderedSettlingCoefficient * settlingFactor;

    for( localIndex c = 0; c < NC; ++c )
    {

      dSettlingFactor_dComponentConcentration[c] = dSingleParticleSettlingVelocity_dComponentConcentration[c] * exp(
        -m_hinderedSettlingCoefficient * proppantConcentration );

    }

    // collisionFactor
    // Collision model (We need to check the other models)
    if( m_isCollisionalSlip )
    {
      real64 lambda = m_collisionAlpha - pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta );

      real64 dLambda_dC = -m_collisionBeta * pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta - 1.0 );

      if( proppantConcentration < m_slipConcentration )
        dLambda_dC = -dLambda_dC;

      collisionFactor = (lambda - 1.0) / (1.0 - proppantConcentration);

      dCollisionFactor_dProppantConcentration = (dLambda_dC  + collisionFactor) / (1.0 - proppantConcentration);
    }

  }
}

void ParticleFluid::ComputeMob( real64 const & concentration,
                                real64 const & aperture,
                                integer & isProppantMobile,
                                real64 & proppantPackPermeability ) const
{

  real64 minAperture;

  // We also need to consider other models/implementations to address proppant "bridging"

  if( concentration < 0.17 )
  {

    minAperture = (1.0 + (m_bridgingFactor - 1.0) / 0.17 * concentration) * m_proppantDiameter;

  }
  else
  {

    minAperture = m_bridgingAperture;

  }


  isProppantMobile = (aperture > minAperture) ? 1 : 0;


  if( !isProppantMobile || concentration >= m_maxProppantConcentration )
  {

    real64 poro = 1.0 - concentration;

    if( poro < 1.0 )
      proppantPackPermeability = m_packPermeabilityCoef * (poro * poro *poro)/((1.0 - poro) * (1.0 - poro));
    else
      proppantPackPermeability = 0.0;
  }
  else
  {

    //we may use an empirical correlation to relate pack peremeability to proppant concentration

    proppantPackPermeability = 0.0;

  }

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParticleFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
