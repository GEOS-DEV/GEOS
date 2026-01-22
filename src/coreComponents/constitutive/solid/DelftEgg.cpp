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
 *  @file DelftEgg.cpp
 */

#include "DelftEgg.hpp"
#include "SolidFields.hpp"

namespace geos
{
using namespace dataRepository;
namespace constitutive
{

DelftEgg::DelftEgg( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent )
{
  // register default values

  registerWrapper( viewKeyStruct::defaultRecompressionIndexString(), &m_defaultRecompressionIndex ).
    setApplyDefaultValue( 2e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Recompresion Index" );

  registerWrapper( viewKeyStruct::defaultVirginCompressionIndexString(), &m_defaultVirginCompressionIndex ).
    setApplyDefaultValue( 5e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Virgin compression index" );

  registerWrapper( viewKeyStruct::defaultCslSlopeString(), &m_defaultCslSlope ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::defaultShapeParameterString(), &m_defaultShapeParameter ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::defaultPreConsolidationPressureString(), &m_defaultPreConsolidationPressure ).
    setApplyDefaultValue( -1.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial preconsolidation pressure" );

  // register fields

  registerField< fields::solid::recompressionIndex >( &m_recompressionIndex );

  registerField< fields::solid::virginCompressionIndex >( &m_virginCompressionIndex );

  registerField< fields::solid::cslSlope >( &m_cslSlope );

  registerField< fields::solid::shapeParameter >( &m_shapeParameter );

  registerField< fields::solid::preConsolidationPressure >( &m_newPreConsolidationPressure );

  registerField< fields::solid::oldPreConsolidationPressure >( &m_oldPreConsolidationPressure );
}


void DelftEgg::allocateConstitutiveData( Group & parent, localIndex const numPts )
{
  m_newPreConsolidationPressure.resize( 0, numPts );
  m_oldPreConsolidationPressure.resize( 0, numPts );

  ElasticIsotropic::allocateConstitutiveData( parent, numPts );
}


void DelftEgg::postInputInitialization()
{
  ElasticIsotropic::postInputInitialization();

  GEOS_THROW_IF( m_defaultCslSlope <= 0,
                 getFullName() << ": Non-positive slope of critical state line detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultShapeParameter < 1.,
                 getFullName() << ": Shape parameter for yield surface must be greater than or equal to one",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultVirginCompressionIndex <= 0,
                 getFullName() << ": Non-positive virgin compression index detected",
                 InputError, getDataContext() );
  GEOS_THROW_IF( m_defaultVirginCompressionIndex <= m_defaultRecompressionIndex,
                 getFullName() << ": Recompression index should exceed virgin recompression index",
                 InputError, getDataContext() );

  // set results as array default values

  getField< fields::solid::oldPreConsolidationPressure >().
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getField< fields::solid::preConsolidationPressure >().
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getField< fields::solid::recompressionIndex >().
    setApplyDefaultValue( m_defaultRecompressionIndex );

  getField< fields::solid::virginCompressionIndex >().
    setApplyDefaultValue( m_defaultVirginCompressionIndex );

  getField< fields::solid::cslSlope >().
    setApplyDefaultValue( m_defaultCslSlope );

  getField< fields::solid::shapeParameter >().
    setApplyDefaultValue( m_defaultShapeParameter );
}


void DelftEgg::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPreConsolidationPressure = m_newPreConsolidationPressure;
  arrayView2d< real64 > oldPreConsolidationPressure = m_oldPreConsolidationPressure;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldPreConsolidationPressure( k, q ) = newPreConsolidationPressure( k, q );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DelftEgg, string const &, Group * const )
}
} /* namespace geos */
