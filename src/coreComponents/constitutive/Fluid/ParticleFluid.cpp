/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

ParticleFluid::ParticleFluid( std::string const & name, Group * const parent ):
  ParticleFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::fluidDensityString, &m_fluidDensity, false )->
    setApplyDefaultValue(1000.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid density");

  registerWrapper( viewKeyStruct::proppantDensityString, &m_proppantDensity, false )->
    setApplyDefaultValue(1400.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Proppant density");

    registerWrapper( viewKeyStruct::fluidViscosityString, &m_fluidViscosity, false )->
    setApplyDefaultValue(0.001)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fluid viscosity");

  registerWrapper( viewKeyStruct::proppantDiameterString, &m_proppantDiameter, false )->
    setApplyDefaultValue(200e-6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Proppant diameter");

  registerWrapper( viewKeyStruct::hinderedSettlingCoefficientString, &m_hinderedSettlingCoefficient, false )->
    setApplyDefaultValue(5.9)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Hindered settling coefficient");

  registerWrapper( viewKeyStruct::collisionAlphaString, &m_collisionAlpha, false )->
    setApplyDefaultValue(1.27)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Collision alpha coefficient");

  registerWrapper( viewKeyStruct::slipConcentrationString, &m_slipConcentration, false )->
    setApplyDefaultValue(0.1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Slip concentration");      

  registerWrapper( viewKeyStruct::collisionBetaString, &m_collisionBeta, false )->
    setApplyDefaultValue(1.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Collision beta coefficient");

  registerWrapper( viewKeyStruct::bridgingFactorString, &m_bridgingFactor, false )->
    setApplyDefaultValue(3.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Bridging factor");

  registerWrapper( viewKeyStruct::sphericityString, &m_sphericity, false )->
    setApplyDefaultValue(1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Sphericity");        

}

ParticleFluid::~ParticleFluid() = default;

void ParticleFluid::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                             localIndex const numConstitutivePointsPerParentIndex )
{
  ParticleFluidBase::AllocateConstitutiveData(parent, numConstitutivePointsPerParentIndex);

}

void
ParticleFluid::DeliverClone( string const & name,
                                            Group * const parent,
                                            std::unique_ptr<ConstitutiveBase> & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique<ParticleFluid>( name, parent );
  }
  ParticleFluidBase::DeliverClone( name, parent, clone );
  ParticleFluid * const newConstitutiveRelation = dynamic_cast<ParticleFluid *>(clone.get());

  newConstitutiveRelation->m_fluidDensity     = this->m_fluidDensity;
  newConstitutiveRelation->m_proppantDensity  = this->m_proppantDensity;
  newConstitutiveRelation->m_fluidViscosity   = this->m_fluidViscosity;
  newConstitutiveRelation->m_proppantDiameter = this->m_proppantDiameter;
  newConstitutiveRelation->m_singleParticleSettlingVelocity = this->m_singleParticleSettlingVelocity;
  newConstitutiveRelation->m_hinderedSettlingCoefficient = this->m_hinderedSettlingCoefficient;
  newConstitutiveRelation->m_collisionAlpha = this->m_collisionAlpha;
  newConstitutiveRelation->m_slipConcentration = this->m_slipConcentration;
  newConstitutiveRelation->m_collisionBeta = this->m_collisionBeta;

  newConstitutiveRelation->m_bridgingFactor = this->m_bridgingFactor;
  newConstitutiveRelation->m_sphericity = this->m_sphericity;
  
}

void ParticleFluid::PostProcessInput()
{
  ParticleFluidBase::PostProcessInput();

  m_singleParticleSettlingVelocity = 9.81 * (m_proppantDensity - m_fluidDensity ) * m_proppantDiameter * m_proppantDiameter / (18.0 * m_fluidViscosity);

  m_packPermeabilityCoef = pow(m_sphericity * m_proppantDiameter, 2.0) / 180.0;

  m_bridgingAperture = m_bridgingFactor * m_proppantDiameter;
  
}

void ParticleFluid::PointUpdate( real64 const & concentration, localIndex const k)
{
  Compute( concentration, m_settlingFactor[k], m_dSettlingFactor_dConc[k], m_collisionFactor[k], m_dCollisionFactor_dConc[k]);

}

void ParticleFluid::PointUpdateMob( real64 const & concentration, real64 const & aperture, localIndex const k)
{
  ComputeMob( concentration, aperture, m_isProppantMobile[k], m_proppantPackPermeability[k]);

}  

void ParticleFluid::Compute( real64 const & concentration,
			     real64 & settlingFactor,
			     real64 & dSettlingFactor_dConc,
			     real64 & collisionFactor,
			     real64 & dCollisionFactor_dConc ) const
{

  settlingFactor = 0.0;

  dSettlingFactor_dConc = 0.0;

  collisionFactor = 0.0;

  dCollisionFactor_dConc = 0.0;

  if(concentration >= 0.0 && concentration < m_maxProppantConcentration)
    {

      // settlingFactor
      // Stokes + hindered settling model
      
      settlingFactor = m_singleParticleSettlingVelocity * exp(-m_hinderedSettlingCoefficient * concentration);

      dSettlingFactor_dConc = -m_hinderedSettlingCoefficient * settlingFactor;
  
      // collisionFactor
      // Collision model (We need to check the other models)
      if (m_isCollisionalSlip)
        {
          real64 lambda = m_collisionAlpha - pow(fabs(concentration - m_slipConcentration), m_collisionBeta);

          real64 dLambda_dC = -m_collisionBeta * pow(fabs(concentration - m_slipConcentration), m_collisionBeta - 1.0);

          if(concentration < m_slipConcentration)
            dLambda_dC = -dLambda_dC;

          collisionFactor = (1.0 - lambda * concentration) / (1.0 - concentration);

          dCollisionFactor_dConc = (collisionFactor - dLambda_dC * concentration - lambda) / (1.0 - concentration);
        }

    }
}

void ParticleFluid::ComputeMob( real64 const & concentration,
				real64 const & aperture,
				bool & isProppantMobile,
				real64 & proppantPackPermeability ) const
{

  real64 minAperture;

  // We also need to consider other models/implementations to address proppant "bridging"
  
  if(concentration < 0.17)
    {

      minAperture = (1.0 + (m_bridgingFactor - 1.0) / 0.17 * concentration) * m_proppantDiameter;

    }
  else
    {

      minAperture = m_bridgingAperture;

    }


  isProppantMobile = (aperture > minAperture) ? 1 : 0;


  if(!isProppantMobile || concentration >= m_maxProppantConcentration)
    {

      real64 poro = 1.0 - concentration;

      if(poro < 1.0)      
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
