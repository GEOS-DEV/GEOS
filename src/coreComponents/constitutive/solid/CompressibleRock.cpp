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
 * @file CompressibleRock.cpp
 */

#include "CompressibleRock.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

CompressibleRock::CompressibleRock( string const & name, Group * const parent ):
  RockBase( name, parent ),
  m_referencePorosity()
{
  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure for solid compressibility" );

  registerWrapper( viewKeyStruct::defaultRefererencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

CompressibleRock::~CompressibleRock() = default;

std::unique_ptr< ConstitutiveBase >
CompressibleRock::deliverClone( string const & name,
                                         Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void CompressibleRock::allocateConstitutiveData( dataRepository::Group & parent,
                                                          localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void CompressibleRock::postProcessInput()
{
  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRock, string const &, Group * const )
}
} /* namespace geosx */
