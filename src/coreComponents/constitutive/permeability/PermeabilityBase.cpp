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
 * @file PermeabilityBase.cpp
 */

#include "../permeability/PermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PermeabilityBase::PermeabilityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::permeabilityString(), &m_permeability ).setPlotLevel( PlotLevel::LEVEL_0 );
}

PermeabilityBase::~PermeabilityBase() = default;

std::unique_ptr< ConstitutiveBase >
PermeabilityBase::deliverClone( string const & name,
                                Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PermeabilityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  m_permeability.resize( 0, numConstitutivePointsPerParentIndex, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PermeabilityBase::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PermeabilityBase, string const &, Group * const )
}
} /* namespace geosx */
