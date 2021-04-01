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
 * @file PorosityBase.cpp
 */

#include "PorosityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PorosityBase::PorosityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::porosityString(), &m_porosity ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::porosityOldString(), &m_porosityOld ).
    setPlotLevel( PlotLevel::LEVEL_3 );
  registerWrapper( viewKeyStruct::dPorosity_dPressureString(), &m_dPorosity_dPressure );
}

PorosityBase::~PorosityBase() = default;

std::unique_ptr< ConstitutiveBase >
PorosityBase::deliverClone( string const & name,
                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void PorosityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  m_porosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_porosityOld.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PorosityBase::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorosityBase, string const &, Group * const )
}
} /* namespace geosx */
