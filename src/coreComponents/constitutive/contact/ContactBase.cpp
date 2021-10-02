/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ContactBase.cpp
 */

#include "ContactBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

ContactBase::ContactBase( string const & name,
                          Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_penaltyStiffness( 0.0 )
{
  registerWrapper( viewKeyStruct::penaltyStiffnessString(), &m_penaltyStiffness ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of the penetration penalty stiffness. Units of Pressure/length" );
}

ContactBase::~ContactBase()
{}

} /* end namespace constitutive */

} /* end namespace geosx */
