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
 * @file PorosityModel.cpp
 */

#include "PorosityModel.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PorosityModel::PorosityModel( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
}

PorosityModel::~PorosityModel() = default;

std::unique_ptr< ConstitutiveBase >
PorosityModel::deliverClone( string const & name,
                        Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PorosityModel::allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );


}

void PorosityModel::postProcessInput()
{
  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorosityModel, string const &, Group * const )
}
} /* namespace geosx */
