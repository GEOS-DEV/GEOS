/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file CoulombFriction.cpp
 */

#include "CoulombFriction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

CoulombFriction::CoulombFriction( string const & name, Group * const parent ):
  FrictionBase( name, parent )
{
  registerWrapper( viewKeyStruct::shearStiffnessString(), &m_shearStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the shear elastic stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::cohesionString(), &m_cohesion ).
    setApplyDefaultValue( 1e6 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setRegisteringObjects( this->getName() ).
    setDescription( "Cohesion for each cell" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).
    setApplyDefaultValue( 0.4 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setRegisteringObjects( this->getName() ).
    setDescription( "Friction coefficient for each cell" );

  registerWrapper( viewKeyStruct::elasticSlipString(), &m_elasticSlip ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Elastic Slip" );

  registerWrapper( viewKeyStruct::defaultCohesionString(), &m_defaultCohesion ).
    setApplyDefaultValue( 1e6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default cohesion value" );

  registerWrapper( viewKeyStruct::defaultFrictionCoefficientString(), &m_defaultFrictionCoefficient ).
    setApplyDefaultValue( 0.4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default friction coefficient value" );
}

void CoulombFriction::postInputInitialization()
{
  GEOS_THROW_IF( m_defaultFrictionCoefficient < 0.0,
                 getFullName() << ": The provided default friction coefficient is less than zero. Value: " << m_defaultFrictionCoefficient,
                 InputError, getDataContext() );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::cohesionString() ).
    setApplyDefaultValue( m_defaultCohesion );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::frictionCoefficientString() ).
    setApplyDefaultValue( m_defaultFrictionCoefficient );
}

void CoulombFriction::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_elasticSlip.resize( 0, 2 );

  FrictionBase::allocateConstitutiveData( parent, numPts );
}

CoulombFrictionUpdates CoulombFriction::createKernelUpdates() const
{
  return CoulombFrictionUpdates( m_displacementJumpThreshold,
                                 m_shearStiffness,
                                 m_cohesion.toViewConst(),
                                 m_frictionCoefficient.toViewConst(),
                                 m_elasticSlip );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CoulombFriction, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
