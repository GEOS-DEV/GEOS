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
 * @file PermeabilityBase.cpp
 */

#include "../permeability/PermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PermeabilityBase::PermeabilityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{}

PermeabilityBase::~PermeabilityBase() = default;

std::unique_ptr< ConstitutiveBase >
PermeabilityBase::deliverClone( string const & name,
                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PermeabilityBase::allocateConstitutiveData( dataRepository::Group * const parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_permeability.resize( parent->size(), numConstitutivePointsPerParentIndex, 3 );
}

void PermeabilityBase::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PermeabilityBase, string const &, Group * const )
}
} /* namespace geosx */
