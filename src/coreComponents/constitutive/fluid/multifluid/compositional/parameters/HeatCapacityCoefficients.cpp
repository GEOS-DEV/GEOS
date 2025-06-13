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
 * @file HeatCapacityCoefficients.cpp
 */

#include "HeatCapacityCoefficients.hpp"
#include "ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "dataRepository/InputFlags.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

HeatCapacityCoefficients::HeatCapacityCoefficients( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters >
HeatCapacityCoefficients::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< HeatCapacityCoefficients >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< HeatCapacityCoefficients >( std::move( parameters ) );
}

void HeatCapacityCoefficients::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::enthalpyReferenceTemperatureString(), &m_referenceTemperature ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "The reference temperature for enthalpy calculation" );

  fluid->registerWrapper( viewKeyStruct::referenceEnthalpyString(), &m_referenceEnthalpy ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "The enthalpy of each component at the reference temperature" );

  fluid->registerWrapper( viewKeyStruct::componentHeatCapacityCoefficientsString(), &m_coefficients ).
    setInputFlag( dataRepository::InputFlags::REQUIRED ).
    setDescription( "The polynomial coefficients for the specific heat capacity of each component in each phase" );
}

void HeatCapacityCoefficients::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                            ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  integer const numPhases = fluid->numFluidPhases();
  integer const numComps = fluid->numFluidComponents();

  // If the reference temperatures are given, then there must be as many as there are components
  if( m_referenceTemperature.empty())
  {
    m_referenceTemperature.resize( numComps );
    m_referenceTemperature.zero();
  }
  else
  {
    GEOS_THROW_IF_NE_MSG( m_referenceTemperature.size(), numComps,
                          GEOS_FMT( "{}: '{}' there must be as many reference temperatures provided as there are components",
                                    fluid->getFullName(),
                                    viewKeyStruct::enthalpyReferenceTemperatureString() ),
                          InputError );
  }

  // If the reference enthalpies are given, then there must be as many as there are components
  if( m_referenceEnthalpy.empty())
  {
    m_referenceEnthalpy.resize( numComps );
    m_referenceEnthalpy.zero();
  }
  else
  {
    GEOS_THROW_IF_NE_MSG( m_referenceEnthalpy.size(), numComps,
                          GEOS_FMT( "{}: '{}' there must be as many reference enthalpy values provided as there are components",
                                    fluid->getFullName(),
                                    viewKeyStruct::referenceEnthalpyString() ),
                          InputError );
  }

  integer const dim0 = m_coefficients.size( 0 );
  integer const dim1 = m_coefficients.size( 1 );
  integer const dim2 = m_coefficients.size( 2 );

  // First dimension must be equal to number of phases
  GEOS_THROW_IF_NE_MSG( dim0, numPhases,
                        GEOS_FMT( "{}: '{}' the first dimension must be equal to the number of phases {}",
                                  fluid->getFullName(),
                                  viewKeyStruct::componentHeatCapacityCoefficientsString(),
                                  numPhases ),
                        InputError );
  // Second dimension must be equal to number of components
  GEOS_THROW_IF_NE_MSG( dim1, numComps,
                        GEOS_FMT( "{}: '{}' the second dimension must be equal to the number of components {}",
                                  fluid->getFullName(),
                                  viewKeyStruct::componentHeatCapacityCoefficientsString(),
                                  numComps ),
                        InputError );
  // Third dimension must be equal to 4
  GEOS_THROW_IF_NE_MSG( dim2, 5,
                        GEOS_FMT( "{}: '{}' the third dimension must be equal 5",
                                  fluid->getFullName(),
                                  viewKeyStruct::componentHeatCapacityCoefficientsString() ),
                        InputError );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
