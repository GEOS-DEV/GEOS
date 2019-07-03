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
  * @file SlurryFluidBase.cpp
  */

#include "SlurryFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SlurryFluidBase::SlurryFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent )
{

  RegisterViewWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default value for density.");

  RegisterViewWrapper( viewKeyStruct::defaultViscosityString, &m_defaultViscosity, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default value for viscosity.");

  RegisterViewWrapper( viewKeyStruct::densityString, &m_density, false )->setPlotLevel( PlotLevel::LEVEL_0 );

  RegisterViewWrapper( viewKeyStruct::fluidDensityString, &m_fluidDensity, false )->setPlotLevel( PlotLevel::LEVEL_0 );

  RegisterViewWrapper( viewKeyStruct::dDens_dPresString, &m_dDens_dPres, false );
  RegisterViewWrapper( viewKeyStruct::dDens_dConcString, &m_dDens_dConc, false );
  RegisterViewWrapper( viewKeyStruct::dFluidDens_dPresString, &m_dFluidDens_dPres, false );  

  RegisterViewWrapper( viewKeyStruct::viscosityString, &m_viscosity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  RegisterViewWrapper( viewKeyStruct::dVisc_dPresString, &m_dVisc_dPres, false );
  RegisterViewWrapper( viewKeyStruct::dVisc_dConcString, &m_dVisc_dConc, false );
  
}

SlurryFluidBase::~SlurryFluidBase() = default;

void SlurryFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
  this->getWrapper< array2d<real64> >(viewKeyStruct::densityString)->setApplyDefaultValue(m_defaultDensity);
  this->getWrapper< array2d<real64> >(viewKeyStruct::viscosityString)->setApplyDefaultValue(m_defaultViscosity);

}

void SlurryFluidBase::AllocateConstitutiveData( ManagedGroup * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_fluidDensity.resize( parent->size(), numConstitutivePointsPerParentIndex );  
  m_dDens_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDens_dConc.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dFluidDens_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );  
  

  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dPres.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dVisc_dConc.resize( parent->size(), numConstitutivePointsPerParentIndex );
  
}


void
SlurryFluidBase::DeliverClone( string const & name,
                               ManagedGroup * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOS_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  SlurryFluidBase * const newConstitutiveRelation = dynamic_cast<SlurryFluidBase *>(clone.get());

  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_defaultViscosity = m_defaultViscosity;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_fluidDensity = m_fluidDensity;  
  newConstitutiveRelation->m_viscosity = m_viscosity;
  newConstitutiveRelation->m_dDens_dPres = m_dDens_dPres;
  newConstitutiveRelation->m_dDens_dConc = m_dDens_dConc;
  newConstitutiveRelation->m_dFluidDens_dPres = m_dFluidDens_dPres;  
  newConstitutiveRelation->m_dVisc_dPres = m_dVisc_dPres;
  newConstitutiveRelation->m_dVisc_dConc = m_dVisc_dConc;  
}
} //namespace constitutive

} //namespace geosx
