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
  * @file SingleFluidBase.cpp
  */

#include "SingleFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

SingleFluidBase::SingleFluidBase( std::string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{

  registerWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default value for density.");

  registerWrapper( viewKeyStruct::defaultViscosityString, &m_defaultViscosity, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Default value for viscosity.");

  registerWrapper( viewKeyStruct::densityString, &m_density, false )->setPlotLevel( PlotLevel::LEVEL_0 );

  registerWrapper( viewKeyStruct::dDens_dPresString, &m_dDensity_dPressure, false );

  registerWrapper( viewKeyStruct::viscosityString, &m_viscosity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dVisc_dPresString, &m_dViscosity_dPressure, false );
}

SingleFluidBase::~SingleFluidBase() = default;

void SingleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
  this->getWrapper< array2d<real64> >(viewKeyStruct::densityString)->setApplyDefaultValue(m_defaultDensity);
  this->getWrapper< array2d<real64> >(viewKeyStruct::viscosityString)->setApplyDefaultValue(m_defaultViscosity);

}

void SingleFluidBase::AllocateConstitutiveData( Group * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
}


void
SingleFluidBase::DeliverClone( string const & name,
                               Group * const parent,
                               std::unique_ptr<ConstitutiveBase> & clone ) const
{
  GEOS_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  SingleFluidBase * const newConstitutiveRelation = dynamic_cast<SingleFluidBase *>(clone.get());

  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_defaultViscosity = m_defaultViscosity;
  newConstitutiveRelation->m_density = m_density;
  newConstitutiveRelation->m_viscosity = m_viscosity;
  newConstitutiveRelation->m_dDensity_dPressure = m_dDensity_dPressure;
  newConstitutiveRelation->m_dViscosity_dPressure = m_dViscosity_dPressure;
}
} //namespace constitutive

} //namespace geosx
