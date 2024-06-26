/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseThermalConductivity.cpp
 */

#include "SinglePhaseThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseThermalConductivity::SinglePhaseThermalConductivity( string const & name, Group * const parent ):
  SinglePhaseThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultThermalConductivityComponentsString(), &m_defaultThermalConductivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz diagonal components of the default thermal conductivity tensor [J/(s.m.K)]" );

  registerWrapper( viewKeyStruct::thermalConductivityGradientComponentsString(), &m_thermalConductivityGradientComponents ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( {0.0,0.0,0.0} ).
    setDescription( "xx, yy, and zz diagonal components of the thermal conductivity gradient tensor w.r.t. temperature [J/(s.m.K^2)]" );

  registerWrapper( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "The reference temperature at which the conductivity components are equal to the default values" );
}

std::unique_ptr< ConstitutiveBase >
SinglePhaseThermalConductivity::deliverClone( string const & name,
                                              Group * const parent ) const
{
  return SinglePhaseThermalConductivityBase::deliverClone( name, parent );
}

void SinglePhaseThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                               localIndex const numConstitutivePointsPerParentIndex )
{
  SinglePhaseThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SinglePhaseThermalConductivity::postProcessInput()
{
  GEOS_THROW_IF( m_defaultThermalConductivityComponents[0] <= 0 ||
                 m_defaultThermalConductivityComponents[1] <= 0 ||
                 m_defaultThermalConductivityComponents[2] <= 0,
                 GEOS_FMT( "{}: the components of the default thermal conductivity tensor must be strictly positive",
                           getFullName() ),
                 InputError );

}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SinglePhaseThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geos
