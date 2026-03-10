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
 * @file BrineSalinity.cpp
 */

#include "BrineSalinity.hpp"

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

BrineSalinity::BrineSalinity( std::unique_ptr< ModelParameters > parameters ):
  ModelParameters( std::move( parameters ) )
{}

std::unique_ptr< ModelParameters > BrineSalinity::create( std::unique_ptr< ModelParameters > parameters )
{
  if( parameters && parameters->get< BrineSalinity >() != nullptr )
  {
    return parameters;
  }
  return std::make_unique< BrineSalinity >( std::move( parameters ) );
}

void BrineSalinity::registerParametersImpl( MultiFluidBase * fluid )
{
  fluid->registerWrapper( viewKeyStruct::salinityString(), &m_salinity ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_salinity ).
    setDescription( "Brine salinity" );

  fluid->registerWrapper( viewKeyStruct::waterCompressibilityString(), &m_waterCompressibility ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_waterCompressibility ).
    setDescription( "The compressibility of pure water." );

  fluid->registerWrapper( viewKeyStruct::saltMolarWeightString(), &m_saltMolarWeight ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_saltMolarWeight ).
    setDescription( "The molar weight for the salt component" );
}

void BrineSalinity::postInputInitializationImpl( MultiFluidBase const * fluid,
                                                 ComponentProperties const & componentProperties )
{
  GEOS_UNUSED_VAR( componentProperties );

  real64 constexpr epsilon = MultiFluidConstants::epsilon;

  GEOS_THROW_IF_LT_MSG( m_salinity, -epsilon,
                        GEOS_FMT( "{}: invalid salinity value provided in '{}'. Salinity should not be negative.",
                                  fluid->getFullName(),
                                  viewKeyStruct::salinityString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_waterCompressibility, epsilon,
                        GEOS_FMT( "{}: invalid salinity value provided in '{}'. Compressibility should be positive.",
                                  fluid->getFullName(),
                                  viewKeyStruct::waterCompressibilityString() ),
                        InputError );

  GEOS_THROW_IF_LT_MSG( m_saltMolarWeight, epsilon,
                        GEOS_FMT( "{}: invalid salt molar weight value provided in '{}'. Molar weight should be positive.",
                                  fluid->getFullName(),
                                  viewKeyStruct::saltMolarWeightString() ),
                        InputError );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos
