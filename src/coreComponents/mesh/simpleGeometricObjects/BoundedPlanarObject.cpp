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
 * @file BoundedPlanarObject.cpp
 */

#include "BoundedPlanarObject.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

BoundedPlanarObject::BoundedPlanarObject( const string & name, Group * const parent ):
  SimpleGeometricObjectBase( name, parent ),
  m_normal{ 0.0, 0.0, 1.0 },
  m_lengthVector{ 0.0, 0.0, 0.0 },
  m_widthVector{ 0.0, 0.0, 0.0 }
{
  registerWrapper( viewKeyStruct::normalString(), &m_normal ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)" );

  registerWrapper( viewKeyStruct::mLengthVectorString(), &m_lengthVector ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent vector defining the orthonormal basis along with the normal." );

  registerWrapper( viewKeyStruct::mWidthVectorString(), &m_widthVector ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Tangent vector defining the orthonormal basis along with the normal." );  
}

BoundedPlanarObject::~BoundedPlanarObject()
{}

//REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase, BoundedPlanarObject, string const &, Group * const )

} /* namespace geosx */