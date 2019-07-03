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

ParticleFluidBase::ParticleFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent )
{

  RegisterViewWrapper( viewKeyStruct::settlingFactorString, &m_settlingFactor, false );
  RegisterViewWrapper( viewKeyStruct::dSettlingFactor_dConcString, &m_dSettlingFactor_dConc, false );

  RegisterViewWrapper( viewKeyStruct::collisionFactorString, &m_collisionFactor, false );
  RegisterViewWrapper( viewKeyStruct::dCollisionFactor_dConcString, &m_dCollisionFactor_dConc, false );

  RegisterViewWrapper( viewKeyStruct::maxProppantConcentrationString, &m_maxProppantConcentration, false )->
    setApplyDefaultValue(0.6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Max proppant concentration");    

  RegisterViewWrapper( viewKeyStruct::isProppantMobileString, &m_isProppantMobile, false );

  RegisterViewWrapper( viewKeyStruct::proppantPackPermeabilityString, &m_proppantPackPermeability, false );  
  
}

ParticleFluidBase::~ParticleFluidBase() = default;

void ParticleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
}

void ParticleFluidBase::AllocateConstitutiveData( ManagedGroup * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_settlingFactor.resize( parent->size());
  m_dSettlingFactor_dConc.resize( parent->size());

  m_collisionFactor.resize( parent->size());
  m_dCollisionFactor_dConc.resize( parent->size());

  m_isProppantMobile.resize( parent->size());  
  m_proppantPackPermeability.resize( parent->size());  
  
}


void
ParticleFluidBase::DeliverClone( string const & name,
                               ManagedGroup * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOS_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  ParticleFluidBase * const newConstitutiveRelation = dynamic_cast<ParticleFluidBase *>(clone.get());

  newConstitutiveRelation->m_settlingFactor = m_settlingFactor;
  newConstitutiveRelation->m_dSettlingFactor_dConc = m_dSettlingFactor_dConc;

  newConstitutiveRelation->m_collisionFactor = m_collisionFactor;
  newConstitutiveRelation->m_dCollisionFactor_dConc = m_dCollisionFactor_dConc;  

  newConstitutiveRelation->m_maxProppantConcentration = this->m_maxProppantConcentration;

  newConstitutiveRelation->m_isProppantMobile = this->m_isProppantMobile;

  newConstitutiveRelation->m_proppantPackPermeability = this->m_proppantPackPermeability;      
  
}
} //namespace constitutive

} //namespace geosx
