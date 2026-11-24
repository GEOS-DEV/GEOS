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
  // register wrappers

  registerWrapper( viewKeyStruct::shearStiffnessString(), &m_shearStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the shear elastic stiffness. Units of Pressure/length" );

  registerWrapper( viewKeyStruct::elasticSlipString(), &m_elasticSlip ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Elastic Slip" );

  registerWrapper( viewKeyStruct::defaultCohesionString(), &m_defaultCohesion ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default cohesion value" );

  registerWrapper( viewKeyStruct::defaultFrictionCoefficientString(), &m_defaultFrictionCoefficient ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default friction coefficient value" );

  // register fields
  registerField< fields::contact::cohesion >( &m_cohesion );
  registerField< fields::contact::frictionCoefficient >( &m_frictionCoefficient );

}

void CoulombFriction::postInputInitialization()
{
  GEOS_THROW_IF( m_defaultCohesion < 0.0,
                 GEOS_FMT( ": The provided default cohesion is less than zero. Value: {}",
                           m_defaultCohesion ),
                 InputError, getDataContext() );

  GEOS_THROW_IF( m_defaultFrictionCoefficient < 0.0,
                 GEOS_FMT( ": The provided default friction coefficient is less than zero. Value: {}",
                           m_defaultFrictionCoefficient ),
                 InputError, getDataContext() );

  getField< fields::contact::cohesion >().
    setApplyDefaultValue( m_defaultCohesion );

  getField< fields::contact::frictionCoefficient >().
    setApplyDefaultValue( m_defaultFrictionCoefficient );
}

void CoulombFriction::initializePostInitialConditionsPreSubGroups()
{
  FrictionBase::initializePostInitialConditionsPreSubGroups();

  localIndex negativeCohesionCount = 0;
  localIndex negativeFrictionCoefficientCount = 0;
  localIndex firstNegativeCohesionIndex = -1;
  localIndex firstNegativeFrictionCoefficientIndex = -1;

  for( localIndex k = 0; k < m_cohesion.size(); ++k )
  {
    if( m_cohesion[k] < 0.0 )
    {
      ++negativeCohesionCount;
      if( firstNegativeCohesionIndex < 0 )
      {
        firstNegativeCohesionIndex = k;
      }
    }

    if( m_frictionCoefficient[k] < 0.0 )
    {
      ++negativeFrictionCoefficientCount;
      if( firstNegativeFrictionCoefficientIndex < 0 )
      {
        firstNegativeFrictionCoefficientIndex = k;
      }
    }
  }

  GEOS_THROW_IF( negativeCohesionCount > 0 || negativeFrictionCoefficientCount > 0,
                 GEOS_FMT( "Negative Coulomb properties detected in per-cell data: "
                           "cohesion count = {} (first index = {}), "
                           "friction coefficient count = {} (first index = {})",
                           negativeCohesionCount,
                           firstNegativeCohesionIndex,
                           negativeFrictionCoefficientCount,
                           firstNegativeFrictionCoefficientIndex ),
                 InputError,
                 getDataContext() );
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
