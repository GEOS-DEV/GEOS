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
 * @file CoupledSolid.cpp
 */

#include "CoupledSolid.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::CoupledSolid( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_solidModel( nullptr ),
  m_porosityModel( nullptr ),
  m_solidModelName(),
  m_porosityModelName()
{
  registerWrapper( viewKeyStruct::solidModelNameString(), &m_solidModelName ).
      setInputFlag( InputFlags::REQUIRED ).
      setDescription( "Name of the solid model." );

    registerWrapper( viewKeyStruct::porosityModelNameString(), &m_porosityModelName ).
        setInputFlag( InputFlags::REQUIRED ).
        setDescription( "Name of the porosity model." );
}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
CoupledSolid< SOLID_TYPE, PORO_TYPE >::~CoupledSolid()
{}

template< typename SOLID_TYPE,
          typename PORO_TYPE >
void CoupledSolid< SOLID_TYPE, PORO_TYPE >::postProcessInput()
{
  m_solidModel = &this->getParent().template getGroup< SOLID_TYPE >( m_solidModelName );
  m_porosityModel = &this->getParent().template getGroup< PORO_TYPE >( m_porosityModelName );
}

}
} /* namespace geosx */
