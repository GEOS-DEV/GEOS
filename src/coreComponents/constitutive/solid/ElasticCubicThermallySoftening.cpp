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
 *  @file ElasticCubicThermallySoftening.cpp
 */

#include "ElasticCubicThermallySoftening.hpp"

namespace geos
{
using namespace dataRepository;

namespace constitutive
{

ElasticCubicThermallySoftening::ElasticCubicThermallySoftening( string const & name, Group * const parent ):
  ElasticCubic( name, parent ),
//   m_referenceC11(),
//   m_referenceC12(),
//   m_referenceC44(),
  m_firstOrderC11ThermalCoefficient( 0.0 ),
  m_firstOrderC12ThermalCoefficient( 0.0 ),
  m_firstOrderC44ThermalCoefficient( 0.0 ),
  m_secondOrderC11ThermalCoefficient( 0.0 ),
  m_secondOrderC12ThermalCoefficient( 0.0 ),
  m_secondOrderC44ThermalCoefficient( 0.0 ),
  m_referenceTemperature( 300.0 ),
  m_temperature()
{
//   registerWrapper< real64 >( viewKeyStruct::referenceC11String() ).
//     setApplyDefaultValue( m_defaultC11 ).
//     setInputFlag( InputFlags::FALSE ).
//     setDescription( "Reference C11" );

//   registerWrapper< real64 >( viewKeyStruct::referenceC12String() ).
//     setApplyDefaultValue( m_defaultC12 ).
//     setInputFlag( InputFlags::FALSE ).
//     setDescription( "Reference C12" );

//   registerWrapper< real64 >( viewKeyStruct::referenceC44String() ).
//     setApplyDefaultValue( m_defaultC44 ).
//     setInputFlag( InputFlags::FALSE ).
//     setDescription( "Reference C44" );

  registerWrapper< real64 >( viewKeyStruct::firstOrderC11ThermalCoefficientString(), &m_firstOrderC11ThermalCoefficient ).
    setApplyDefaultValue( m_firstOrderC11ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "First order C11 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::firstOrderC12ThermalCoefficientString(), &m_firstOrderC12ThermalCoefficient ).
    setApplyDefaultValue( m_firstOrderC12ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "First order C12 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::firstOrderC44ThermalCoefficientString(), &m_firstOrderC44ThermalCoefficient ).
    setApplyDefaultValue( m_firstOrderC44ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "First order C44 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::secondOrderC11ThermalCoefficientString(), &m_secondOrderC11ThermalCoefficient ).
    setApplyDefaultValue( m_secondOrderC11ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Second order C11 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::secondOrderC12ThermalCoefficientString(), &m_secondOrderC12ThermalCoefficient ).
    setApplyDefaultValue( m_secondOrderC12ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Second order C12 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::secondOrderC44ThermalCoefficientString(), &m_secondOrderC44ThermalCoefficient ).
    setApplyDefaultValue( m_secondOrderC44ThermalCoefficient ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Second order C44 thermal coefficient" );

  registerWrapper< real64 >( viewKeyStruct::referenceTemperatureString(), &m_referenceTemperature ).
    setApplyDefaultValue( m_referenceTemperature ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference temperature" );

  registerWrapper< array1d< real64 > >( viewKeyStruct::temperatureString(), &m_temperature ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Temperature" );
}

ElasticCubicThermallySoftening::~ElasticCubicThermallySoftening()
{}

void ElasticCubicThermallySoftening::allocateConstitutiveData( dataRepository::Group & parent, 
                                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticCubic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_temperature.resize( 0 );
}

void ElasticCubicThermallySoftening::postProcessInput()
{
  ElasticCubic::postProcessInput();
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticCubicThermallySoftening, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
