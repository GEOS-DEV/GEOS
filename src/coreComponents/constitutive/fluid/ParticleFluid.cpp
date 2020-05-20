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
  else if( str == "Turbulence" )
  {
    return ParticleFluid::ParticleSettlingModel::Turbulence;
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

  newConstitutiveRelation->m_sphericity = this->m_sphericity;

  newConstitutiveRelation->m_particleSettlingModelString = this->m_particleSettlingModelString;
  newConstitutiveRelation->m_particleSettlingModel = this->m_particleSettlingModel;

}

void ParticleFluid::PostProcessInput()
{
  ParticleFluidBase::PostProcessInput();

  GEOSX_ERROR_IF( m_particleSettlingModelString.empty(),
                  "Invalid particle settling model in ParticleFluid " );

  GEOSX_ERROR_IF( m_proppantDensity < 500.0,
                  "Invalid proppantDensity in ParticleFluid, which must >= 500.0 " );

  GEOSX_ERROR_IF( m_proppantDiameter < 10e-6,
                  "Invalid proppantDiameter in ParticleFluid, which must >= 10e-6 " );

  GEOSX_ERROR_IF( m_hinderedSettlingCoefficient< 0.0 || m_hinderedSettlingCoefficient > 10.0,
                  "Invalid hinderedSettlingCoefficient in ParticleFluid, which must between 0 and 10 " );

  GEOSX_ERROR_IF( m_collisionAlpha < 1.0,
                  "Invalid collisionAlpha in ParticleFluid, which must >= 1 " );

  GEOSX_ERROR_IF( m_collisionBeta < 0.0,
                  "Invalid collisionBeta in ParticleFluid, which must >= 0" );

  GEOSX_ERROR_IF( m_slipConcentration > 0.3,
                  "Invalid slipConcentration in ParticleFluid, which must <= 0.3" );

  m_packPermeabilityCoef = pow( m_sphericity * m_proppantDiameter, 2.0 ) / 180.0;

  m_particleSettlingModel = stringToParticleSettlingModel( m_particleSettlingModelString );

}

void ParticleFluid::BatchUpdate( arrayView1d< real64 const > const & )
{
  GEOSX_ERROR( "BatchUpdate for ParticleFluid is not implemented" );
}


void ParticleFluid::PointUpdate( localIndex const GEOSX_UNUSED_PARAM( NC ),
                                 real64 const & GEOSX_UNUSED_PARAM( proppantConcentration ),
                                 arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( componentConcentration ),
                                 arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( nIndex ),
                                 arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( KIndex ),
                                 real64 const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                 real64 const & GEOSX_UNUSED_PARAM( dFluidDensity_dPressure ),
                                 arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                                 localIndex const GEOSX_UNUSED_PARAM( k ) )
{
  GEOSX_ERROR( "Function: PointUpdate is used for ParticlePowerLawFluid and is not implemented" );

}

void ParticleFluid::PointUpdate( localIndex const NC, real64 const & proppantConcentration, real64 const & fluidDensity, real64 const & dFluidDensity_dPressure,
                                 arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration, real64 const & fluidViscosity,
                                 real64 const & dFluidViscosity_dPressure, arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
                                 localIndex const k )
{

  Compute( NC, proppantConcentration, fluidDensity, dFluidDensity_dPressure, dFluidDensity_dComponentConcentration, fluidViscosity, dFluidViscosity_dPressure,
           dFluidViscosity_dComponentConcentration, m_settlingFactor[k], m_dSettlingFactor_dPressure[k], m_dSettlingFactor_dProppantConcentration[k],
           m_dSettlingFactor_dComponentConcentration[k], m_collisionFactor[k], m_dCollisionFactor_dProppantConcentration[k] );

}


void ParticleFluid::Compute( localIndex const NC,
                             real64 const & proppantConcentration,
                             real64 const & fluidDensity,
                             real64 const & GEOSX_UNUSED_PARAM( dFluidDensity_dPressure ),
                             arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDensity_dComponentConcentration ),
                             real64 const & fluidViscosity,
                             real64 const & GEOSX_UNUSED_PARAM( dFluidViscosity_dPressure ),
                             arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidViscosity_dComponentConcentration ),
                             real64 & settlingFactor,
                             real64 & dSettlingFactor_dPressure,
                             real64 & dSettlingFactor_dProppantConcentration,
                             arraySlice1d< real64 > const & dSettlingFactor_dComponentConcentration,
                             real64 & collisionFactor,
                             real64 & dCollisionFactor_dProppantConcentration ) const
{

  real64 singleParticleSettlingVelocity = 0.0;

  static real64 constCoef = 9.81 * m_proppantDiameter * m_proppantDiameter / 18.0;

  switch( m_particleSettlingModel )
  {
    case ParticleSettlingModel::Stokes:
    {

      singleParticleSettlingVelocity = constCoef * (m_proppantDensity - fluidDensity ) / fluidViscosity;
    }
    break;
    case ParticleSettlingModel::Intermediate:
    {

      singleParticleSettlingVelocity = 0.2 * pow( m_proppantDiameter, 1.18 ) * pow( 9.81 * (m_proppantDensity - fluidDensity) /fluidDensity, 0.72 ) * pow(
        fluidDensity / fluidViscosity, 0.45 );

    }
    break;
    case ParticleSettlingModel::Turbulence:
    {

      singleParticleSettlingVelocity = 1.74 * pow( m_proppantDiameter, 0.5 ) * pow( 9.81 * (m_proppantDensity - fluidDensity) /fluidDensity, 0.5 );

    }
    break;
    default:
      GEOSX_ERROR( "Particle settling model type not supported" );
  }


  settlingFactor = 0.0;

  dSettlingFactor_dPressure = 0.0;

  dSettlingFactor_dProppantConcentration = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {

    dSettlingFactor_dComponentConcentration[c] = 0.0;

  }

  collisionFactor = 0.0;

  dCollisionFactor_dProppantConcentration = 0.0;

  if( proppantConcentration >= 0.0 && proppantConcentration < m_maxProppantConcentration )
  {

    // settlingFactor

    settlingFactor = singleParticleSettlingVelocity * exp( -m_hinderedSettlingCoefficient * proppantConcentration );

    // collisionFactor
    // Collision model (We need to check the other models)
    if( m_isCollisionalSlip )
    {
      real64 lambda = m_collisionAlpha - pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta );

      real64 dLambda_dC = -m_collisionBeta * pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta - 1.0 );

      if( proppantConcentration < m_slipConcentration )
        dLambda_dC = -dLambda_dC;

      collisionFactor = (lambda - 1.0) / (1.0 - proppantConcentration);


    }

  }

}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParticleFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
