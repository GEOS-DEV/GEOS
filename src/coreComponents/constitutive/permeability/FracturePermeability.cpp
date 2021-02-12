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

#include "../permeability/FracturePermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


FracturePermeability::FracturePermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{}

FracturePermeability::~FracturePermeability() = default;

std::unique_ptr< ConstitutiveBase >
FracturePermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void FracturePermeability::allocateConstitutiveData( dataRepository::Group * const parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_permeability.resize( parent->size(), numConstitutivePointsPerParentIndex, 1 );
  m_dPerm_dAperture.resize( parent->size(), numConstitutivePointsPerParentIndex, 1 );
}

void FracturePermeability::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, FracturePermeability, string const &, Group * const )

void FracturePermeabilityUpdate::compute( real64 const & effectiveAperture,
                                          arraySlice1d< real64 > const & permeability,
                                          arraySlice1d< real64 > const & dPerm_dAperture )
{
  permeability[0] = effectiveAperture*effectiveAperture*effectiveAperture /12;
  dPerm_dAperture[0]  = effectiveAperture*effectiveAperture / 12;
}

}
} /* namespace geosx */
