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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

SolidBase::SolidBase( string const & name, Group * const parent ):
  ContinuumBase( name, parent ),
  m_thermalExpansionCoefficient()
{
  registerWrapper( viewKeyStruct::defaultThermalExpansionCoefficientString(), &m_defaultThermalExpansionCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Default Linear Thermal Expansion Coefficient of the Solid Rock Frame" );

  // register fields
  registerField< fields::solid::thermalExpansionCoefficient >( &m_thermalExpansionCoefficient );
}


void SolidBase::postInputInitialization()
{
  ContinuumBase::postInputInitialization();

  getField< fields::solid::thermalExpansionCoefficient >().
    setApplyDefaultValue( m_defaultThermalExpansionCoefficient );
}


void SolidBase::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  ContinuumBase::allocateConstitutiveData( parent, numPts );

  m_thermalExpansionCoefficient.resize( 0 );
}


void SolidBase::saveConvergedState() const
{
  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView3d< real64 const, solid::STRESS_USD > newStress = m_newStress;
  arrayView3d< real64, solid::STRESS_USD > oldStress = m_oldStress;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      LvArray::tensorOps::copy< 6 >( oldStress[k][q], newStress[k][q] );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SolidBase, string const &, Group * const )
} /* namespace constitutive */
} /* namespace geos */
