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
 * @file ProppantPorosityModel.cpp
 */

#include "ProppantPorosity.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ProppantPorosity::ProppantPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent ),
  m_maxProppantConcentration()
{
  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Maximum proppant concentration " );
}

ProppantPorosity::~ProppantPorosity() = default;

std::unique_ptr< ConstitutiveBase >
ProppantPorosity::deliverClone( string const & name,
                                Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void ProppantPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void ProppantPorosity::postProcessInput()
{
  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ProppantPorosity, string const &, Group * const )
}
} /* namespace geosx */
