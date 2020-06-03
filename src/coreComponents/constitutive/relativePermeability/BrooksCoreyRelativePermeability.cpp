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
 * @file BrooksCoreyRelativePermeability.cpp
 */

#include "BrooksCoreyRelativePermeability.hpp"

#include <cmath>

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

//START_SPHINX_INCLUDE_00

BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( std::string const & name,
                                                                  Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString, &m_phaseMinVolumeFraction )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseRelPermExponentString, &m_phaseRelPermExponent )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "MinimumRel perm power law exponent for each phase" );


  registerWrapper( viewKeyStruct::phaseRelPermMaxValueString, &m_phaseRelPermMaxValue )->
    setApplyDefaultValue( 0.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Maximum rel perm value for each phase" );

}

BrooksCoreyRelativePermeability::~BrooksCoreyRelativePermeability()
{}

void
BrooksCoreyRelativePermeability::DeliverClone( string const & name,
                                               Group * const parent,
                                               std::unique_ptr< ConstitutiveBase > & clone ) const
{
  std::unique_ptr< BrooksCoreyRelativePermeability > newModel = std::make_unique< BrooksCoreyRelativePermeability >( name, parent );

  newModel->m_phaseNames = this->m_phaseNames;
  newModel->m_phaseTypes = this->m_phaseTypes;
  newModel->m_phaseOrder = this->m_phaseOrder;

  newModel->m_phaseMinVolumeFraction = this->m_phaseMinVolumeFraction;
  newModel->m_phaseRelPermExponent   = this->m_phaseRelPermExponent;
  newModel->m_phaseRelPermMaxValue   = this->m_phaseRelPermMaxValue;

  newModel->m_satScale = this->m_satScale;

  clone = std::move( newModel );
}


void BrooksCoreyRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  #define COREY_CHECK_INPUT_LENGTH( data, expected, attr ) \
    if( LvArray::integerConversion< localIndex >((data).size()) != LvArray::integerConversion< localIndex >( expected )) \
    { \
      GEOSX_ERROR( "BrooksCoreyRelativePermeability: invalid number of entries in " \
                   << (attr) << " attribute (" \
                   << (data).size() << "given, " \
                   << (expected) << " expected)" ); \
    }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction, NP, viewKeyStruct::phaseMinVolumeFractionString )
  COREY_CHECK_INPUT_LENGTH( m_phaseRelPermExponent, NP, viewKeyStruct::phaseRelPermExponentString )
  COREY_CHECK_INPUT_LENGTH( m_phaseRelPermMaxValue, NP, viewKeyStruct::phaseRelPermMaxValueString )

#undef COREY_CHECK_INPUT_LENGTH

  m_satScale = 1.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    GEOSX_ERROR_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip] );
    m_satScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_ERROR_IF( m_phaseRelPermExponent[ip] < 0.0,
                    "BrooksCoreyRelativePermeability: invalid exponent value: " << m_phaseRelPermExponent[ip] );

    GEOSX_ERROR_IF( m_phaseRelPermMaxValue[ip] < 0.0 || m_phaseRelPermMaxValue[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid maximum value: " << m_phaseRelPermMaxValue[ip] );
  }

  GEOSX_ERROR_IF( m_satScale < 0.0, "BrooksCoreyRelativePermeability: sum of min volume fractions exceeds 1.0" );
}


void BrooksCoreyRelativePermeability::BatchUpdate( arrayView2d< real64 const > const & phaseVolumeFraction )
{

  arrayView1d< real64 const > const & phaseMinVolumeFraction = m_phaseMinVolumeFraction;
  arrayView1d< real64 const > const & phaseRelPermExponent = m_phaseRelPermExponent;
  arrayView1d< real64 const > const & phaseRelPermMaxValue = m_phaseRelPermMaxValue;


  RelativePermeabilityBase::BatchUpdateKernel< BrooksCoreyRelativePermeability >( phaseVolumeFraction,
                                                                                  phaseMinVolumeFraction,
                                                                                  phaseRelPermExponent,
                                                                                  phaseRelPermMaxValue,
                                                                                  m_satScale );
}


void BrooksCoreyRelativePermeability::PointUpdate( arraySlice1d< real64 const > const & phaseVolFraction,
                                                   localIndex const k,
                                                   localIndex const q )
{
  arraySlice1d< real64 > const relPerm           = m_phaseRelPerm[k][q];
  arraySlice2d< real64 > const dRelPerm_dVolFrac = m_dPhaseRelPerm_dPhaseVolFrac[k][q];

  localIndex const NP = numFluidPhases();

  Compute( NP,
           phaseVolFraction,
           relPerm,
           dRelPerm_dVolFrac,
           m_phaseMinVolumeFraction,
           m_phaseRelPermExponent,
           m_phaseRelPermMaxValue,
           m_satScale );
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyRelativePermeability, std::string const &, Group * const )
} // namespace constitutive

} // namespace geosx
