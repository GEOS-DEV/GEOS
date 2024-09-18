/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 *  @file RateAndStateFriction.cpp
 */

#include "RateAndStateFriction.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

RateAndStateFriction::RateAndStateFriction( string const & name, Group * const parent ):
  FrictionBase( name, parent ),
  m_frictionCoefficient()
{}

RateAndStateFriction::~RateAndStateFriction()
{}

void RateAndStateFriction::postInputInitialization()
{}

void RateAndStateFriction::allocateConstitutiveData( Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  FrictionBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

using RateAndStateFrictionUpdates = RateAndStateFriction::KernelWrapper;
RateAndStateFrictionUpdates RateAndStateFriction::createKernelWrapper() const
{
  return RateAndStateFrictionUpdates( m_displacementJumpThreshold,
                                      m_frictionCoefficient, m_a, m_b,
                                      m_Dc, m_V0, m_mu0 );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, RateAndStateFriction, string const &, Group * const )

} /* namespace constitutive */

} /* namespace geos */
