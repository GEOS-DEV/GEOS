/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

BrooksCoreyRelativePermeability::BrooksCoreyRelativePermeability( string const & name,
                                                                  Group * const parent )
  : RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString(), &m_phaseMinVolumeFraction ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseRelPermExponentString(), &m_phaseRelPermExponent ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "MinimumRel perm power law exponent for each phase" );


  registerWrapper( viewKeyStruct::phaseRelPermMaxValueString(), &m_phaseRelPermMaxValue ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum rel perm value for each phase" );

  registerWrapper( viewKeyStruct::volFracScaleString(), &m_volFracScale ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Factor used to scale the phase capillary pressure, defined as: one minus the sum of the phase minimum volume fractions." );

}

BrooksCoreyRelativePermeability::~BrooksCoreyRelativePermeability()
{}

namespace
{

template< typename ARRAY >
void checkInputSize( ARRAY const & array, localIndex const expected, string const & attr )
{
  GEOSX_THROW_IF_NE_MSG( array.size(), expected,
                         "BrooksCoreyRelativePermeability: invalid number of entries in " << attr << " attribute",
                         InputError );

}

}

void BrooksCoreyRelativePermeability::postProcessInput()
{
  RelativePermeabilityBase::postProcessInput();

  localIndex const numPhases = numFluidPhases();

  checkInputSize( m_phaseMinVolumeFraction, numPhases, viewKeyStruct::phaseMinVolumeFractionString() );
  checkInputSize( m_phaseRelPermExponent, numPhases, viewKeyStruct::phaseRelPermExponentString() );
  checkInputSize( m_phaseRelPermMaxValue, numPhases, viewKeyStruct::phaseRelPermMaxValueString() );

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    GEOSX_THROW_IF( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid min volume fraction value: " << m_phaseMinVolumeFraction[ip],
                    InputError );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_THROW_IF( m_phaseRelPermExponent[ip] < 0.0,
                    "BrooksCoreyRelativePermeability: invalid exponent value: " << m_phaseRelPermExponent[ip],
                    InputError );

    GEOSX_THROW_IF( m_phaseRelPermMaxValue[ip] < 0.0 || m_phaseRelPermMaxValue[ip] > 1.0,
                    "BrooksCoreyRelativePermeability: invalid maximum value: " << m_phaseRelPermMaxValue[ip],
                    InputError );
  }

  GEOSX_THROW_IF( m_volFracScale < 0.0, "BrooksCoreyRelativePermeability: sum of min volume fractions exceeds 1.0", InputError );
}

BrooksCoreyRelativePermeability::KernelWrapper BrooksCoreyRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseRelPermExponent,
                        m_phaseRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( ConstitutiveBase, BrooksCoreyRelativePermeability, string const &, Group * const )

} // namespace constitutive

} // namespace geosx
