/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file CeramicDamage.cpp
 */

#include "CeramicDamage.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

CeramicDamage::CeramicDamage( string const & name, Group * const parent ):
  ElasticIsotropic( name, parent ),
  m_damage(),
  m_defaultFailureStress()
{
  // register default values
  registerWrapper( viewKeyStruct::defaultFailureStressString(), &m_defaultFailureStress ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default failure stress" );

  // register fields
  registerWrapper( viewKeyStruct::damageString(), &m_damage ).
    setApplyDefaultValue( 0.0 ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( "Array of quadrature point damage values" );
}


CeramicDamage::~CeramicDamage()
{}


void CeramicDamage::allocateConstitutiveData( dataRepository::Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  ElasticIsotropic::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_damage.resize( 0, numConstitutivePointsPerParentIndex );
}


void CeramicDamage::postProcessInput()
{
  ElasticIsotropic::postProcessInput();

  GEOSX_THROW_IF( m_defaultFailureStress < 0.0, "Invalid failure stress, please specify a non-negative value.", InputError );
}


void CeramicDamage::saveConvergedState() const
{
  SolidBase::saveConvergedState();
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, CeramicDamage, std::string const &, Group * const )
}
} /* namespace geosx */
