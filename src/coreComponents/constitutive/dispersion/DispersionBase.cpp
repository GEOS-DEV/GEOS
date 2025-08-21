/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file DispersionBase.cpp
 */

#include "constitutive/dispersion/DispersionBase.hpp"
#include "constitutive/dispersion/DispersionFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

DispersionBase::DispersionBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{}

void DispersionBase::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  auto subdomain = dynamic_cast< ElementSubRegionBase * >(&parent); // TODO remove

  // NOTE: enforcing 1 quadrature point
  subdomain->registerField< fields::dispersion::dispersivity >( getName(), &m_dispersivity ).
    reference().resizeDimension< 1, 2 >( 1, 3 );
}

} // namespace constitutive

} // namespace geos
