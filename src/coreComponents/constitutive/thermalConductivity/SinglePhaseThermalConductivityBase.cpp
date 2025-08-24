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
 * @file SinglePhaseThermalConductivityBase.cpp
 */

#include "SinglePhaseThermalConductivityBase.hpp"
#include "ThermalConductivityFields.hpp"
#include "mesh/ElementSubRegionBase.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseThermalConductivityBase::SinglePhaseThermalConductivityBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent )
{}

void SinglePhaseThermalConductivityBase::allocateConstitutiveData( Group & parent,
                                                                   localIndex const numPts )
{
  auto subregion = dynamic_cast< ElementSubRegionBase * >( &parent ); // TODO remove

  // NOTE: enforcing 1 quadrature point
  subregion->registerField< fields::thermalconductivity::effectiveConductivity >( getName(), &m_effectiveConductivity ).
    reference().resizeDimension< 1, 2 >( 1, 3 );
  subregion->registerField< fields::thermalconductivity::dEffectiveConductivity_dT >( getName(), &m_dEffectiveConductivity_dT ).
    reference().resizeDimension< 1, 2 >( 1, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numPts );
}

} // namespace constitutive

} // namespace geos
