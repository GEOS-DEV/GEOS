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
 * @file ParallelPlatesPermeability.cpp
 */

#include "ParallelPlatesPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ParallelPlatesPermeability::ParallelPlatesPermeability( string const & name, Group * const parent ):
  FracturePermeabilityBase( name, parent )
{}

ParallelPlatesPermeability::~ParallelPlatesPermeability() = default;

std::unique_ptr< ConstitutiveBase >
ParallelPlatesPermeability::deliverClone( string const & name,
                                          Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void ParallelPlatesPermeability::allocateConstitutiveData( dataRepository::Group * const parent,
                                                           localIndex const numConstitutivePointsPerParentIndex )
{
  FracturePermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, ParallelPlatesPermeability, string const &, Group * const )

}
} /* namespace geosx */
