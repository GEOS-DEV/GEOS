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
  * @file ParticleFluidBase.cpp
  */

#include "ParticleFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ParticleFluidBase::ParticleFluidBase( std::string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{

  registerWrapper( viewKeyStruct::settlingFactorString, &m_settlingFactor, false );
  registerWrapper( viewKeyStruct::dSettlingFactor_dPressureString, &m_dSettlingFactor_dPressure, false );
  registerWrapper( viewKeyStruct::dSettlingFactor_dProppantConcentrationString, &m_dSettlingFactor_dProppantConcentration, false );
    registerWrapper( viewKeyStruct::dSettlingFactor_dComponentConcentrationString, &m_dSettlingFactor_dComponentConcentration, false );  

  registerWrapper( viewKeyStruct::collisionFactorString, &m_collisionFactor, false );
  registerWrapper( viewKeyStruct::dCollisionFactor_dProppantConcentrationString, &m_dCollisionFactor_dProppantConcentration, false );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString, &m_maxProppantConcentration, false )->
    setApplyDefaultValue(0.6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Max proppant concentration");    

  registerWrapper( viewKeyStruct::isProppantMobileString, &m_isProppantMobile, false );

  registerWrapper( viewKeyStruct::isCollisionalSlipString, &m_isCollisionalSlip, false )->
      setApplyDefaultValue(0)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Whether the collisional component of the slip velocity is considered");

  registerWrapper( viewKeyStruct::proppantPackPermeabilityString, &m_proppantPackPermeability, false );  
  
}

ParticleFluidBase::~ParticleFluidBase() = default;

void ParticleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
}

void ParticleFluidBase::AllocateConstitutiveData( Group * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_settlingFactor.resize( parent->size());
  m_dSettlingFactor_dPressure.resize( parent->size());
  m_dSettlingFactor_dProppantConcentration.resize( parent->size());
  m_dSettlingFactor_dComponentConcentration.resize( parent->size(), MAX_NUM_COMPONENTS);  

  m_collisionFactor.resize( parent->size());
  m_dCollisionFactor_dProppantConcentration.resize( parent->size());

  m_isProppantMobile.resize( parent->size());  
  m_proppantPackPermeability.resize( parent->size());  
  
}


void
ParticleFluidBase::DeliverClone( string const & name,
                               Group * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOS_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  ParticleFluidBase * const newConstitutiveRelation = dynamic_cast<ParticleFluidBase *>(clone.get());

  newConstitutiveRelation->m_settlingFactor = m_settlingFactor;
  newConstitutiveRelation->m_dSettlingFactor_dPressure = m_dSettlingFactor_dPressure;
  newConstitutiveRelation->m_dSettlingFactor_dProppantConcentration = m_dSettlingFactor_dProppantConcentration;
  newConstitutiveRelation->m_dSettlingFactor_dComponentConcentration = m_dSettlingFactor_dComponentConcentration;  

  newConstitutiveRelation->m_collisionFactor = m_collisionFactor;
  newConstitutiveRelation->m_dCollisionFactor_dProppantConcentration = m_dCollisionFactor_dProppantConcentration;  

  newConstitutiveRelation->m_maxProppantConcentration = this->m_maxProppantConcentration;

  newConstitutiveRelation->m_isProppantMobile = this->m_isProppantMobile;

  newConstitutiveRelation->m_isCollisionalSlip = this->m_isCollisionalSlip;

  newConstitutiveRelation->m_proppantPackPermeability = this->m_proppantPackPermeability;      
  
}
} //namespace constitutive

} //namespace geosx
