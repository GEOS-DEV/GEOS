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
 * @file PressureDependentPorosity.cpp
 */

#include "PressureDependentPorosity.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

PressureDependentPorosity::PressureDependentPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Solid compressibility" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Reference pressure for fluid compressibility" );
}

PressureDependentPorosity::~PressureDependentPorosity() = default;

std::unique_ptr< ConstitutiveBase >
PressureDependentPorosity::deliverClone( string const & name,
                                         Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PressureDependentPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                                          localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PressureDependentPorosity::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PressureDependentPorosity, string const &, Group * const )
}
} /* namespace geosx */
