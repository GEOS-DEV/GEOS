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
 * @file SinglePhaseConstantThermalConductivity.cpp
 */

#include "SinglePhaseConstantThermalConductivity.hpp"

namespace geos
{

using namespace dataRepository;

namespace constitutive
{

SinglePhaseConstantThermalConductivity::SinglePhaseConstantThermalConductivity( string const & name, Group * const parent ):
  SinglePhaseThermalConductivityBase( name, parent )
{
  registerWrapper( viewKeyStruct::thermalConductivityComponentsString(), &m_thermalConductivityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy, and zz components of a diagonal thermal conductivity tensor [J/(s.m.K)]" );
}

std::unique_ptr< ConstitutiveBase >
SinglePhaseConstantThermalConductivity::deliverClone( string const & name,
                                                      Group * const parent ) const
{
  return SinglePhaseThermalConductivityBase::deliverClone( name, parent );
}

void SinglePhaseConstantThermalConductivity::allocateConstitutiveData( dataRepository::Group & parent,
                                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SinglePhaseThermalConductivityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void SinglePhaseConstantThermalConductivity::postProcessInput()
{
  GEOS_THROW_IF( m_thermalConductivityComponents[0] <= 0 ||
                 m_thermalConductivityComponents[1] <= 0 ||
                 m_thermalConductivityComponents[2] <= 0,
                 GEOS_FMT( "{}: the components of the thermal conductivity tensor must be strictly positive",
                           getFullName() ),
                 InputError );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, SinglePhaseConstantThermalConductivity, string const &, Group * const )

} // namespace constitutive

} // namespace geos
