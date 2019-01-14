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

SingleFluidBase::SingleFluidBase( std::string const & name, ManagedGroup * const parent )
  : ConstitutiveBase( name, parent )
{
  RegisterViewWrapper( viewKeyStruct::densityString, &m_density, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  RegisterViewWrapper( viewKeyStruct::dDens_dPresString, &m_dDensity_dPressure, false );

  RegisterViewWrapper( viewKeyStruct::viscosityString, &m_viscosity, false )->setPlotLevel( PlotLevel::LEVEL_0 );
  RegisterViewWrapper( viewKeyStruct::dVisc_dPresString, &m_dViscosity_dPressure, false );
}

SingleFluidBase::~SingleFluidBase() = default;

void SingleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
}

void SingleFluidBase::AllocateConstitutiveData( ManagedGroup * const parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
}


} //namespace constitutive

} //namespace geosx
