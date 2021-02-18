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
 * @file FracturePermeability.cpp
 */

#include "../permeability/CarmanKozenyPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


CarmanKozenyPermeability::CarmanKozenyPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{}

CarmanKozenyPermeability::~CarmanKozenyPermeability() = default;

std::unique_ptr< ConstitutiveBase >
CarmanKozenyPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void CarmanKozenyPermeability::allocateConstitutiveData( dataRepository::Group * const parent,
                                                         localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  m_dPerm_dPorosity.resize( parent->size(), numConstitutivePointsPerParentIndex, 3 );
}

void CarmanKozenyPermeability::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CarmanKozenyPermeability, string const &, Group * const )

}
} /* namespace geosx */
