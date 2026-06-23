/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file Liquid.cpp
 */

#include "Liquid.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

Liquid::Liquid( string const & name, Group * const parent ):
  ContinuumBase( name, parent ),
  m_defaultBulkModulus( 1.0 ),
  m_bulkModulus(),
  m_defaultViscosity( 0.0 ),
  m_viscosity(),
  m_jacobian(),
  m_velocityGradient()
{
  // register default values
  registerWrapper( viewKeyStruct::defaultBulkModulusString(), &m_defaultBulkModulus ).
    setApplyDefaultValue( m_defaultBulkModulus ).
    setInputFlag( InputFlags::REQUIRED).
    setDescription( "Default bulk modulus" );

  registerWrapper( viewKeyStruct::defaultViscosityString(), &m_defaultViscosity ).
    setApplyDefaultValue( m_defaultViscosity ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default viscosity" );

  // register fields
  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setApplyDefaultValue( m_defaultBulkModulus ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Bulk modulus" );

  registerWrapper( viewKeyStruct::viscosityString(), &m_viscosity ).
    setApplyDefaultValue( m_defaultViscosity ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Viscosity" );

  registerWrapper( viewKeyStruct::jacobianString(), &m_jacobian ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Jacobian" );

  registerWrapper( viewKeyStruct::velocityGradientString(), &m_velocityGradient ).
    setApplyDefaultValue( 0.0).
    setInputFlag( InputFlags::FALSE ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Velocity gradient" );
}


Liquid::~Liquid()
{}


void Liquid::allocateConstitutiveData( dataRepository::Group & parent,
                                    localIndex const numConstitutivePointsPerParentIndex )
{
  ContinuumBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_bulkModulus.resize( 0 );
  m_viscosity.resize( 0 );
  m_jacobian.resize( 0, numConstitutivePointsPerParentIndex );
  m_velocityGradient.resize( 0, 3, 3 );
}


void Liquid::postInputInitialization()
{
  ContinuumBase::postInputInitialization();

  GEOS_THROW_IF( m_defaultBulkModulus <= 0.0, "Reference pressure must be greater than 0.0", InputError );
  GEOS_THROW_IF( m_defaultViscosity < 0.0, "Reference temperature must be greater or equal than 0.0", InputError );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, Liquid, std::string const &, Group * const )
}
} /* namespace geos */
