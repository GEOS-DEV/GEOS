/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file CamClay.cpp
 */

#include "CamClay.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

CamClay::CamClay( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_defaultRefPressure(),
  m_defaultRefStrainVol(),
  m_defaultRecompressionIndex(),
  m_defaultVirginCompressionIndex(),
  m_defaultCslSlope(),
  m_defaultShapeParameter(),
  m_defaultPreConsolidationPressure(),
  m_refPressure(),
  m_refStrainVol(),
  m_recompressionIndex(),
  m_virginCompressionIndex(),
  m_cslSlope(),
  m_shapeParameter(),
  m_newPreConsolidationPressure(),
  m_oldPreConsolidationPressure()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultRefPressureString(), &m_defaultRefPressure ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::defaultRefStrainVolString(), &m_defaultRefStrainVol ).
    setApplyDefaultValue( 1e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Reference volumetric strain" );

  registerWrapper( viewKeyStruct::defaultRecompressionIndexString(), &m_defaultRecompressionIndex ).
    setApplyDefaultValue( 2e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Recompresion index" );

  registerWrapper( viewKeyStruct::defaultVirginCompressionIndexString(), &m_defaultVirginCompressionIndex ).
    setApplyDefaultValue( 5e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Virgin Compression index" );

  registerWrapper( viewKeyStruct::defaultCslSlopeString(), &m_defaultCslSlope ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::defaultShapeParameterString(), &m_defaultShapeParameter ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::defaultPreConsolidationPressureString(), &m_defaultPreConsolidationPressure ).
    setApplyDefaultValue( 5e6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial preconsolidation pressure" );

  // register fields

  registerWrapper( viewKeyStruct::refPressureString(), &m_refPressure ).
    setApplyDefaultValue( -1 ).
    setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::refStrainVolString(), &m_refStrainVol ).
    setApplyDefaultValue( -1 ).
    setDescription( "eference volumetric strain" );

  registerWrapper( viewKeyStruct::recompressionIndexString(), &m_recompressionIndex ).
    setApplyDefaultValue( -1 ).
    setDescription( "Recompression index" );

  registerWrapper( viewKeyStruct::virginCompressionIndexString(), &m_virginCompressionIndex ).
    setApplyDefaultValue( -1 ).
    setDescription( "Virgin compression index" );

  registerWrapper( viewKeyStruct::cslSlopeString(), &m_cslSlope ).
    setApplyDefaultValue( -1 ).
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::shapeParameterString(), &m_shapeParameter ).
    setApplyDefaultValue( -1 ).
    setDescription( "Shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::newPreConsolidationPressureString(), &m_newPreConsolidationPressure ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_3 ).
    setDescription( "New preconsolidation pressure" );

  registerWrapper( viewKeyStruct::oldPreConsolidationPressureString(), &m_oldPreConsolidationPressure ).
    setApplyDefaultValue( -1 ).
    setDescription( "Old preconsolidation pressure" );
}


CamClay::~CamClay()
{}


void CamClay::allocateConstitutiveData( dataRepository::Group & parent,
                                        localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPreConsolidationPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldPreConsolidationPressure.resize( 0, numConstitutivePointsPerParentIndex );

  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void CamClay::postProcessInput()
{
  ElasticIsotropic::postProcessInput();

  GEOSX_THROW_IF( m_defaultCslSlope <= 0, "Negative slope of critical state line detected", InputError );
  GEOSX_THROW_IF( m_defaultRecompressionIndex <= 0, "Negative recompresion index detected", InputError );
  GEOSX_THROW_IF( m_defaultVirginCompressionIndex <= 0, "Negative virgin compression index detected", InputError );
  GEOSX_THROW_IF( m_defaultVirginCompressionIndex <= m_defaultRecompressionIndex, "Recompression index should exceed virgin recompression index", InputError );


  // set results as array default values

  getWrapper< array2d< real64 > >( viewKeyStruct::oldPreConsolidationPressureString() ).
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getWrapper< array2d< real64 > >( viewKeyStruct::newPreConsolidationPressureString() ).
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getWrapper< array1d< real64 > >( viewKeyStruct::refPressureString() ).
    setApplyDefaultValue( m_defaultRefPressure );

  getWrapper< array1d< real64 > >( viewKeyStruct::refStrainVolString() ).
    setApplyDefaultValue( m_defaultRefStrainVol );

  getWrapper< array1d< real64 > >( viewKeyStruct::recompressionIndexString() ).
    setApplyDefaultValue( m_defaultRecompressionIndex );

  getWrapper< array1d< real64 > >( viewKeyStruct::virginCompressionIndexString() ).
    setApplyDefaultValue( m_defaultVirginCompressionIndex );

  getWrapper< array1d< real64 > >( viewKeyStruct::cslSlopeString() ).
    setApplyDefaultValue( m_defaultCslSlope );

  getWrapper< array1d< real64 > >( viewKeyStruct::shapeParameterString() ).
    setApplyDefaultValue( m_defaultShapeParameter );
}


void CamClay::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPreConsolidationPressure = m_newPreConsolidationPressure;
  arrayView2d< real64 > oldPreConsolidationPressure = m_oldPreConsolidationPressure;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldPreConsolidationPressure( k, q ) = newPreConsolidationPressure( k, q );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CamClay, std::string const &, Group * const )
}
} /* namespace geosx */
